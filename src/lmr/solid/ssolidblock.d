/**
 * ssolidblock.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module ssolidblock;

import std.stdio;
import std.string : toStringz;
import std.algorithm;
import std.parallelism;
import std.conv;
import std.format;
import std.array;
import std.math;
import ntypes.complex;
import nm.number;
import std.json;
import util.lua;
import util.lua_service;
import nm.rsla;
import nm.smla;
import geom.luawrap.luasgrid : checkStructuredGrid;
import json_helper;
import lua_helper;
import gzip;
import geom;
import globalconfig;
import globaldata : dedicatedConfig;
import fluidblock;
import sfluidblock;
import solidblock;
import solidfvcell;
import solidfvinterface;
import solidfvvertex;
import solidbc;
import block;
import jacobian;

import lmr.solid.solidthermalmodel;

enum ghost_cell_start_id = 1_000_000_000;

class SSolidBlock : SolidBlock {
public:
    size_t n_ghost_cell_layers;
    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax;
    size_t jmin, jmax;
    size_t kmin, kmax;
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
private:
    StructuredGrid grid; // for reading and writing

    size_t _nidim;
    size_t _njdim;
    size_t _nkdim;

    SolidFVCell[] _ctr;
    SolidFVInterface[] _ifi;
    SolidFVInterface[] _ifj;
    SolidFVInterface[] _ifk;
    SolidFVVertex[] _vtx;
    SolidFVInterface[] _sifi;
    SolidFVInterface[] _sifj;
    SolidFVInterface[] _sifk;

public:
    // Constructors
    this(int id, size_t nicell, size_t njcell, size_t nkcell, string label)
    {
        super(id, label);
        myConfig = dedicatedConfig[id];
        this.n_ghost_cell_layers = GlobalConfig.n_ghost_cell_layers;
        this.nicell = nicell;
        this.njcell = njcell;
        this.nkcell = nkcell;
        fillInOtherSizeData();
    }

    this(int id, JSONValue jsonData)
    {
        this.id = id;
        nicell = getJSONint(jsonData, "nic", 0);
        njcell = getJSONint(jsonData, "njc", 0);
        nkcell = getJSONint(jsonData, "nkc", 0);
        label = getJSONstring(jsonData, "label", "");
        this(id, nicell, njcell, nkcell, label);
        active = getJSONbool(jsonData, "active", true);
    }

    this(lua_State* L)
    {
        auto grid = checkStructuredGrid(L, 2);
        nicell = grid.niv - 1;
        njcell = grid.njv - 1;
        nkcell = grid.nkv - 1;
        if (GlobalConfig.dimensions == 2) nkcell = 1;
        const size_t ncells = nicell*njcell*nkcell;
        cells.length = ncells;
        super(-1, "nil");
        // check where the initial temperature is coming from
        double T_init;
        bool lua_fn = false;
        if (lua_isnumber(L, 3)) {
            T_init = luaL_checknumber(L, 3);
        }
        else if (lua_isfunction(L, 3)) {
            lua_fn = true;
        }
        auto stm = lmr.solid.solidthermalmodel.initSolidThermalModel(L, 4);
        Vector3 pos;
        number volume, iLen, jLen, kLen;
        size_t[3] ijk;
        size_t i, j, k;
        foreach(cell_idx, ref cell; cells) {
            // get the cell index in terms of structured grid indices
            ijk = cell_index_to_logical_coordinates(cell_idx, nicell, njcell);
            i = ijk[0]; j = ijk[1]; k = ijk[2];
            // get the cell data
            Vector3 p000 = *grid[i,j,k];
            Vector3 p100 = *grid[i+1,j,k];
            Vector3 p110 = *grid[i+1,j+1,k];
            Vector3 p010 = *grid[i,j+1,k];
            if (GlobalConfig.dimensions == 2) {
                number xyplane_area;
                xyplane_quad_cell_properties(p000, p100, p110, p010, pos, xyplane_area, iLen, jLen, kLen);
                volume = xyplane_area * ((GlobalConfig.axisymmetric) ? pos.y : to!number(1.0) );
            }
            else if (GlobalConfig.dimensions == 3) {
                Vector3 p001 = *grid[i,j,k+1];
                Vector3 p101 = *grid[i+1,j,k+1];
                Vector3 p111 = *grid[i+1,j+1,k+1];
                Vector3 p011 = *grid[i,j+1,k+1];
                hex_cell_properties(p000, p100, p110, p010, p001, p101, p111, p011,
                                    GlobalConfig.true_centroids, pos, volume, iLen, jLen, kLen);
            }
            else {
                throw new Exception("GlobalConfig.dimensions not 2 or 3.");
            }
            if (lua_fn) {
                // Now grab temperature via Lua function call.
                lua_pushvalue(L, 3);
                lua_pushnumber(L, pos.x);
                lua_pushnumber(L, pos.y);
                lua_pushnumber(L, pos.z);
                if (lua_pcall(L, 3, 1, 0) != 0) {
                    string errMsg = "Error in Lua function call for setting initial temperature in solid\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    luaL_error(L, errMsg.toStringz);
                }
                if (lua_isnumber(L, -1)) {
                    T_init = luaL_checknumber(L, -1);
                }
                else {
                    string errMsg = "Error in Lua function call for setting initial temperature in solid.\n";
                    errMsg ~= "A single number value should be returned.";
                    luaL_error(L, errMsg.toStringz);
                }
            }
            // make the cell
            cell = new SolidFVCell(pos, volume, T_init, stm, to!int(cell_idx));
        }
        lua_settop(L, 0);
    }

    pragma(inline, true)
    @nogc size_t ijk_0n_indices_to_cell_id(size_t i, size_t j, size_t k=0) const
    // ijk indices into the hypothetical block of active cells.
    // where 0<k<nkcell, 0<j<njcell, 0<i<nicell are the indices
    // into the hypothetical block of active cells.
    // This cell_id also the index into the single-dimensional cells array,
    // that is held in the FluidBlock base class.
    // Note that the hypothetical block of active cells is embedded in
    // a larger array that includes surrounding layers of ghost cells.
    {
        return (k*njcell + j)*nicell + i;
    }

    override void initLuaGlobals()
    {
        lua_pushinteger(myL, n_ghost_cell_layers); lua_setglobal(myL, "n_ghost_cell_layers");
        lua_pushinteger(myL, nicell); lua_setglobal(myL, "nicell");
        lua_pushinteger(myL, njcell); lua_setglobal(myL, "njcell");
        lua_pushinteger(myL, nkcell); lua_setglobal(myL, "nkcell");
        // Although we make the helper functions available within
        // the block-specific Lua interpreter, we should use
        // those functions only in the context of the master thread.
        setSampleHelperFunctions(myL);
        // Generally, the sampleFluidCell function can be expected to work only in serial mode.
        // Once it is called from a thread, other than the main thread, it may not
        // have access to properly initialized data for any other block.
    }

    void copy_values_from(SolidFVCell other) {

    }

    override void initBoundaryConditions(JSONValue jsonData)
    {
        foreach (boundary; 0 .. (myConfig.dimensions == 3 ? 6 : 4)) {
            string jsonKey = "face_" ~ face_name[boundary];
            auto bcJsonData = jsonData[jsonKey];
            bc ~= makeSolidBCFromJson(bcJsonData, id, boundary,
                                      nicell, njcell, nkcell);
        }
        foreach (bci; bc) bci.postBCconstruction();
    }

    void fillInOtherSizeData()
    // Helper function for the constructor
    {
        _nidim = nicell + 2 * n_ghost_cell_layers;
        _njdim = njcell + 2 * n_ghost_cell_layers;
        // Indices, in each grid direction for the active cells.
        // These limits are inclusive. The mincell and max cell
        // are both within the active set of cells.
        imin = n_ghost_cell_layers; imax = imin + nicell - 1;
        jmin = n_ghost_cell_layers; jmax = jmin + njcell - 1;
        if ( myConfig.dimensions == 2 ) {
            // In 2D simulations, the k range is from 0 to 0 for the
            // storage arrays of cells and relevant faces.
            if ( nkcell != 1 ) {
                writeln("Warning: inconsistent dimensions nkcell set to 1 for 2D");
                nkcell = 1;
            }
            _nkdim = 1;
            kmin = 0; kmax = 0;
        } else {
            // In 3D simulations the k index is just like the i and j indices.
            _nkdim = nkcell + 2 * n_ghost_cell_layers;
            kmin = n_ghost_cell_layers; kmax = kmin + nkcell - 1;
        }
    } // end fill_in_other_size_data()

    // -- Service methods
    @nogc
    size_t toGlobalIndex(size_t i, size_t j, size_t k) const
    in {
        assert(i < _nidim && j < _njdim && k < _nkdim, "Index out of bounds.");
    }
    do {
        return k * (_njdim * _nidim) + j * _nidim + i;
    }

    size_t[] toIJKIndices(size_t gid) const
    {
        size_t k = gid / (njcell * nicell);
        size_t j = (gid - k * (njcell * nicell)) / nicell;
        size_t i = gid - k * (njcell * nicell) - j * njcell;
        return [i, j, k];
    }

    @nogc ref SolidFVCell getCell(size_t i, size_t j, size_t k=0) { return _ctr[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfi(size_t i, size_t j, size_t k=0) { return _ifi[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfj(size_t i, size_t j, size_t k=0) { return _ifj[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfk(size_t i, size_t j, size_t k=0) { return _ifk[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVVertex getVtx(size_t i, size_t j, size_t k=0) { return _vtx[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifi(size_t i, size_t j, size_t k=0) { return _sifi[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifj(size_t i, size_t j, size_t k=0) { return _sifj[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifk(size_t i, size_t j, size_t k=0) { return _sifk[toGlobalIndex(i,j,k)]; }



    override void assembleArrays()
    {
        if ( myConfig.verbosity_level >= 2 )
            writefln("SSolidBlock.assembleArrays(): Begin for solid_block %d", id);
        // Check for obvious errors.
        if ( _nidim <= 0 || _njdim <= 0 || _nkdim <= 0 ) {
            throw new Error(text("SSolidBlock.assembleArrays(): invalid dimensions nidim=",
                                 _nidim, " njdim=", _njdim, " nkdim=", _nkdim));
        }
        size_t ntot = _nidim * _njdim * _nkdim;
        try {
            // Create the cell and interface objects for the entire block.
            foreach (gid; 0 .. ntot) {
                // we want the ghost cells to have an id that can distinguish them from interior cells
                // note that we will loop through interior cells later and correct their ids
                _ctr ~= new SolidFVCell(myConfig); _ctr[gid].id = to!int(ghost_cell_start_id+gid);
                _ifi ~= new SolidFVInterface(); _ifi[gid].id = gid;
                _ifj ~= new SolidFVInterface(); _ifj[gid].id = gid;
                if ( myConfig.dimensions == 3 ) {
                    _ifk ~= new SolidFVInterface(); _ifk[gid].id = gid;
                }
                _vtx ~= new SolidFVVertex(); _vtx[gid].id = gid;
                _sifi ~= new SolidFVInterface(); _sifi[gid].id = gid;
                _sifj ~= new SolidFVInterface(); _sifj[gid].id = gid;
                if ( myConfig.dimensions == 3 ) {
                    _sifk ~= new SolidFVInterface(); _sifk[gid].id = gid;
                }
            } // gid loop
            // Now, assemble the lists of references to the cells and faces in standard order for a structured grid.
            // These arrays are held by the SolidBlock base class and allow us to handle
            // a structured-grid block much as we would an unstructured-grid block.
            // We also make the cell and face id value consistent with the index in the array.
            // We will depend on this equality in other parts of the solver.
            // We also note that these cells are interior to the block (i.e. not ghost cells)
            //
            // Gather cells
            if (myConfig.dimensions == 2) {
                foreach (j; jmin .. jmax+1) {
                    foreach (i; imin .. imax+1) {
                        cells ~= getCell(i, j);
                        cells[$-1].id = cells.length - 1;
                        cells[$-1].is_ghost = false;
                    }
                }
            } else { // assume 3D
                foreach (k; kmin .. kmax+1) {
                    foreach (j; jmin .. jmax+1) {
                        foreach (i; imin .. imax+1) {
                            cells ~= getCell(i, j, k);
                            cells[$-1].id = cells.length - 1;
                            cells[$-1].is_ghost = false;
                        }
                    }
                }
            } // end if dimensions
            //
            // Gather interfaces
            if ( myConfig.dimensions == 2 ) {
                // East-facing interfaces
                for ( size_t j = jmin; j <= jmax; ++j ) {
                    for ( size_t i = imin; i <= imax + 1; ++i ) {
                        faces ~= getIfi(i,j);
                        faces[$-1].id = faces.length-1;
                    }
                }
                // North-facing interfaces
                for ( size_t j = jmin; j <= jmax + 1; ++j ) {
                    for ( size_t i = imin; i <= imax; ++i ) {
                        faces ~= getIfj(i,j);
                        faces[$-1].id = faces.length-1;
                    }
                }
            } else { // assume myConfig.dimensions == 3 (3D)
                // East-facing interfaces
                for ( size_t k = kmin; k <= kmax; ++k) {
                    for ( size_t j = jmin; j <= jmax; ++j) {
                        for ( size_t i = imin; i <= imax+1; ++i) {
                            faces ~= getIfi(i,j,k);
                            faces[$-1].id = faces.length-1;
                        }
                    }
                }
                // North-facing interfaces
                for ( size_t k = kmin; k <= kmax; ++k) {
                    for ( size_t i = imin; i <= imax; ++i) {
                        for ( size_t j = jmin; j <= jmax+1; ++j) {
                            faces ~= getIfj(i,j,k);
                            faces[$-1].id = faces.length-1;
                        }
                    }
                }
                // Top-facing interfaces
                for ( size_t i = imin; i <= imax; ++i) {
                    for ( size_t j = jmin; j <= jmax; ++j) {
                        for ( size_t k = kmin; k <= kmax+1; ++k) {
                            faces ~= getIfk(i,j,k);
                            faces[$-1].id = faces.length-1;
                        }
                    }
                }
            } // end if dimensions
            // Set references to boundary faces in bc objects.
            SolidFVInterface face;
            // North boundary
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t i = imin; i <= imax; ++i ) {
                    face = getIfj(i, jmax + 1, k);
                    face.is_on_boundary = true;
                    face.bc_id = Face.north;
                    face.bc_idx = bc[Face.north].faces.length;
                    bc[Face.north].faces ~= face;
                    bc[Face.north].outsigns ~= 1;
                }
            }
            // East boundary
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for (size_t j = jmin; j <= jmax; ++j ) {
                    face = getIfi(imax + 1, j, k);
                    face.is_on_boundary = true;
                    face.bc_id = Face.east;
                    face.bc_idx = bc[Face.east].faces.length;
                    bc[Face.east].faces ~= face;
                    bc[Face.east].outsigns ~= 1;
                }
            }
            // South boundary
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t i = imin; i <= imax; ++i ) {
                    face = getIfj(i, jmin, k);
                    face.is_on_boundary = true;
                    face.bc_id = Face.south;
                    face.bc_idx = bc[Face.south].faces.length;
                    bc[Face.south].faces ~= face;
                    bc[Face.south].outsigns ~= -1;
                }
            }
            // West boundary
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t j = jmin; j <= jmax; ++j ) {
                    face = getIfi(imin, j, k);
                    face.is_on_boundary = true;
                    face.bc_id = Face.west;
                    face.bc_idx = bc[Face.west].faces.length;
                    bc[Face.west].faces ~= face;
                    bc[Face.west].outsigns ~= -1;
                }
            }
            if ( myConfig.dimensions == 3 ) {
                // Top boundary
                for ( size_t i = imin; i <= imax; ++i ) {
                    for ( size_t j = jmin; j <= jmax; ++j ) {
                        face = getIfk(i, j, kmax + 1);
                        face.is_on_boundary = true;
                        face.bc_id = Face.top;
                        face.bc_idx = bc[Face.top].faces.length;
                        bc[Face.top].faces ~= face;
                        bc[Face.top].outsigns ~= 1;
                    }
                }
                // Bottom boundary
                for ( size_t i = imin; i <= imax; ++i ) {
                    for ( size_t j = jmin; j <= jmax; ++j ) {
                        face = getIfk(i, j, kmin);
                        face.is_on_boundary = true;
                        face.bc_id = Face.bottom;
                        face.bc_idx = bc[Face.bottom].faces.length;
                        bc[Face.bottom].faces ~= face;
                        bc[Face.bottom].outsigns ~= -1;
                    }
                }
            } // end if dimensions
        } catch (Error err) {
            writeln("Crapped out while assembling solid_block arrays.");
            writefln("nicell=%d njcell=%d nkcell=%d", nicell, njcell, nkcell);
            writefln("nidim=%d njdim=%d nkdim=%d", _nidim, _njdim, _nkdim);
            writeln("Probably ran out of memory.");
            writeln("Be a little less ambitious and try a smaller grid next time.");
            writefln("System message: %s", err.msg);
            throw new Error("SolidBlock.assembleArrays() failed.");
        }
        if ( myConfig.verbosity_level >= 2 ) {
            writefln("Done assembling arrays for %d solid cells.", ntot);
        }
    }
    override void bindFacesAndVerticesToCells()
    {
        size_t kstart, kend;
        if ( myConfig.dimensions == 3 ) {
            kstart = kmin - 1;
            kend = kmax + 1;
        } else {
            kstart = 0;
            kend = 0;
        }
        // With these ranges, we also do the first layer of ghost cells.
        for ( size_t k = kstart; k <= kend; ++k ) {
            for ( size_t j = jmin-1; j <= jmax+1; ++j ) {
                for ( size_t i = imin-1; i <= imax+1; ++i ) {
                    SolidFVCell cell = getCell(i,j,k);
                    cell.iface.length = (myConfig.dimensions == 3) ? 6 : 4;
                    cell.iface[Face.north] = getIfj(i,j+1,k);
                    cell.iface[Face.east] = getIfi(i+1,j,k);
                    cell.iface[Face.south] = getIfj(i,j,k);
                    cell.iface[Face.west] = getIfi(i,j,k);
                    cell.vtx ~= getVtx(i,j,k);
                    cell.vtx ~= getVtx(i+1,j,k);
                    cell.vtx ~= getVtx(i+1,j+1,k);
                    cell.vtx ~= getVtx(i,j+1,k);
                    if ( myConfig.dimensions == 3 ) {
                        cell.iface[Face.top] = getIfk(i,j,k+1);
                        cell.iface[Face.bottom] = getIfk(i,j,k);
                        cell.vtx ~= getVtx(i,j,k+1);
                        cell.vtx ~= getVtx(i+1,j,k+1);
                        cell.vtx ~= getVtx(i+1,j+1,k+1);
                        cell.vtx ~= getVtx(i,j+1,k+1);
                    } // end if
                } // for i
            } // for j
        } // for k
    }

    override void bindCellsToFaces() {
        // When computing the fluxes through an interface it is convenient to have
        // easy access to the cells that are either side of it, this routine
        // gathers that information for each of the interfaces.
        size_t i, j, k;
        SolidFVInterface IFace;
        if ( myConfig.dimensions == 2 ) {
            // East-facing interfaces
            for ( j = jmin; j <= jmax; ++j ) {
                for ( i = imin; i <= imax + 1; ++i ) {
                    IFace = getIfi(i, j);
                    IFace.cellLeft = getCell(i-1, j);
                    IFace.cellRight = getCell(i, j);
                }
            }
            // North-facing interfaces
            for ( j = jmin; j <= jmax + 1; ++j ) {
                for ( i = imin; i <= imax; ++i ) {
                    IFace = getIfj(i, j);
                    IFace.cellLeft = getCell(i, j-1);
                    IFace.cellRight = getCell(i, j);
                }
            }
        } else { // assume myConfig.dimensions == 3 (3D)
            // East-facing interfaces
            for (k = kmin; k <= kmax; ++k) {
                for (j = jmin; j <= jmax; ++j) {
                    for (i = imin; i <= imax+1; ++i) {
                        IFace = getIfi(i, j, k);
                        IFace.cellLeft = getCell(i-1, j, k);
                        IFace.cellRight = getCell(i, j, k);
                    }
                }
            }
            // North-facing interfaces
            for (k = kmin; k <= kmax; ++k) {
                for (i = imin; i <= imax; ++i) {
                    for (j = jmin; j <= jmax+1; ++j) {
                        IFace = getIfj(i, j, k);
                        IFace.cellLeft = getCell(i, j-1, k);
                        IFace.cellRight = getCell(i, j, k);
                      }
                }
            }
            // Top-facing interfaces
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    for (k = kmin; k <= kmax+1; ++k) {
                        IFace = getIfk(i, j, k);
                        IFace.cellLeft = getCell(i, j, k-1);
                        IFace.cellRight = getCell(i, j, k);
                    }
                }
            }
        } // end if/else
    } // end bindFacesToCells()

    override void readGrid(string filename)
    {
        size_t nivtx, njvtx, nkvtx;
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
            writeln("SSolidBlock.readGrid(): Start block ", id);
        }
        grid = new StructuredGrid(filename, "gziptext");
        nivtx = grid.niv; njvtx = grid.njv; nkvtx = grid.nkv;
        if ( myConfig.dimensions == 3 ) {
            if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx-1 != nkcell ) {
                throw new Error(text("For solid_block[", id, "] we have a mismatch in 3D grid size.",
                                     " Have read nivtx=", nivtx, " njvtx=", njvtx,
                                     " nkvtx=", nkvtx));
            }
            for ( size_t k = kmin; k <= kmax+1; ++k ) {
                for ( size_t j = jmin; j <= jmax+1; ++j ) {
                    for ( size_t i = imin; i <= imax+1; ++i ) {
                        auto src_vtx = grid[i-imin,j-jmin,k-kmin];
                        auto vtx = getVtx(i,j,k);
                        vtx.pos.x = src_vtx.x;
                        vtx.pos.y = src_vtx.y;
                        vtx.pos.z = src_vtx.z;
                    } // for i
                } // for j
            } // for k
        } else { // 2D case
            if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx != 1 ) {
                throw new Error(text("For solid_block[", id, "] we have a mismatch in 2D grid size.",
                                     " Have read nivtx=", nivtx, " njvtx=", njvtx,
                                     " nkvtx=", nkvtx));
            }
            for ( size_t j = jmin; j <= jmax+1; ++j ) {
                for ( size_t i = imin; i <= imax+1; ++i ) {
                    auto vtx = getVtx(i,j);
                    auto src_vtx = grid[i-imin,j-jmin];
                    vtx.pos.x = src_vtx.x;
                    vtx.pos.y = src_vtx.y;
                    vtx.pos.z = 0.0;
                } // for i
            } // for j
        }
    } // end read_grid()

    override void writeGrid(string filename, double sim_time)
    {
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
            writeln("SSolidBlock.writeGrid(): Start block ", id);
        }
        size_t kmaxrange;
        if ( myConfig.dimensions == 3 ) {
            kmaxrange = kmax + 1;
        } else { // 2D case
            kmaxrange = kmax;
        }
        for ( size_t k = kmin; k <= kmaxrange; ++k ) {
            for ( size_t j = jmin; j <= jmax+1; ++j ) {
                for ( size_t i = imin; i <= imax+1; ++i ) {
                    auto vtx = getVtx(i,j,k);
                    auto dest_vtx = grid[i-imin,j-jmin,k-kmin];
                    dest_vtx.x = vtx.pos.x;
                    dest_vtx.y = vtx.pos.y;
                    dest_vtx.z = vtx.pos.z;
                } // for i
            } // for j
        } // for k
        grid.write_to_gzip_file(filename);
    } // end write_grid()

    override void readSolution(string filename)
    {
        size_t ni, nj, nk;
        double sim_time;
        if ( myConfig.verbosity_level >= 1 && id == 0 ) {
            writeln("read_solution(): Start solid_block ", id);
        }
        auto byLine = new GzipByLine(filename);
        auto line = byLine.front; byLine.popFront();
        formattedRead(line, " %g", &sim_time);
        line = byLine.front; byLine.popFront();
        // ignore second line; it should be just the names of the variables
        // [TODO] We should test the incoming strings against the current variable names.
        line = byLine.front; byLine.popFront();
        formattedRead(line, "%d %d %d", &ni, &nj, &nk);
        if ( ni != nicell || nj != njcell ||
             nk != ((myConfig.dimensions == 3) ? nkcell : 1) ) {
            throw new Error(text("For solid_block[", id, "] we have a mismatch in solution size.",
                                 " Have read ni=", ni, " nj=", nj, " nk=", nk));
        }
        for ( size_t k = kmin; k <= kmax; ++k ) {
            for ( size_t j = jmin; j <= jmax; ++j ) {
                for ( size_t i = imin; i <= imax; ++i ) {
                    line = byLine.front; byLine.popFront();
                    getCell(i,j,k).scanValuesFromString(line);
                } // for i
            } // for j
        } // for k
    } // end read_solution()

    override void computePrimaryCellGeometricData()
    {
        if ( myConfig.dimensions == 2 ) {
            calcVolumes2D();
            calcFaces2D();
            return;
        } else { // 3D
            calcVolumes3D();
            calcFaces3D();
            return;
            //throw new Error("SSolidBlock.computePrimaryCellGeometryData() not implemented for 3D.");
        }
    } // end compute_primary_cell_geometric_data()

    void calcVolumes2D()
    {
        size_t i, j;
        number xA, yA, xB, yB, xC, yC, xD, yD;
        number xN, yN, xS, yS, xE, yE, xW, yW;
        number vol, max_vol, min_vol, xyarea;
        number dx, dy, dxN, dyN, dxE, dyE;
        number lengthN, lengthE, length_max, length_min, length_cross;
        number max_aspect, aspect_ratio;
        SolidFVCell cell, source_cell, target_cell;

        // Cell layout
        // C-----B     3-----2
        // |     |     |     |
        // |  c  |     |  c  |
        // |     |     |     |
        // D-----A     0-----1

        max_vol = 0.0;
        min_vol = 1.0e30;    /* arbitrarily large */
        max_aspect = 0.0;
        for ( i = imin; i <= imax; ++i ) {
            for ( j = jmin; j <= jmax; ++j ) {
                cell = getCell(i,j);
                // These are the corners.
                xA = cell.vtx[1].pos.x;
                yA = cell.vtx[1].pos.y;
                xB = cell.vtx[2].pos.x;
                yB = cell.vtx[2].pos.y;
                xC = cell.vtx[3].pos.x;
                yC = cell.vtx[3].pos.y;
                xD = cell.vtx[0].pos.x;
                yD = cell.vtx[0].pos.y;
                // Cell area in the (x,y)-plane.
                xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                                (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
                // Cell Centroid.
                cell.pos.x = 1.0 / (xyarea * 6.0) *
                    ((yB - yA) * (xA * xA + xA * xB + xB * xB) +
                     (yC - yB) * (xB * xB + xB * xC + xC * xC) +
                     (yD - yC) * (xC * xC + xC * xD + xD * xD) +
                     (yA - yD) * (xD * xD + xD * xA + xA * xA));
                cell.pos.y = -1.0 / (xyarea * 6.0) *
                    ((xB - xA) * (yA * yA + yA * yB + yB * yB) +
                     (xC - xB) * (yB * yB + yB * yC + yC * yC) +
                     (xD - xC) * (yC * yC + yC * yD + yD * yD) +
                     (xA - xD) * (yD * yD + yD * yA + yA * yA));
                cell.pos.z = 0.0;
                // Cell Volume.
                if ( myConfig.axisymmetric ) {
                    // Volume per radian = centroid y-ordinate * cell area
                    vol = xyarea * cell.pos.y;
                } else {
                    // Assume unit depth in the z-direction.
                    vol = xyarea;
                }
                if (vol < 0.0) {
                    throw new Error(text("Negative cell volume: solid_block ", id,
                                         " vol[", i, " ,", j, "]= ", vol));
                }
                if (vol > max_vol) max_vol = vol;
                if (vol < min_vol) min_vol = vol;
                cell.volume = vol;
                cell.areaxy = xyarea;
            } // j loop
        } // i loop
    } // end calc_volumes_2D()

    void calcFaces2D()
    {
        SolidFVInterface iface;
        size_t i, j;
        number xA, xB, yA, yB, xC, yC;
        number LAB, LBC;

        // East-facing interfaces.
        for (i = imin; i <= imax+1; ++i) {
            for (j = jmin; j <= jmax; ++j) {
                iface = getIfi(i,j);
                // These are the corners.
                xA = getVtx(i,j).pos.x;
                yA = getVtx(i,j).pos.y;
                xB = getVtx(i,j+1).pos.x;
                yB = getVtx(i,j+1).pos.y;
                LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
                if (LAB < 1.0e-9) {
                    writefln("Zero length ifi[%d,%d]: %e", i, j, LAB);
                }
                // Direction cosines for the unit normal.
                iface.n.x = (yB - yA) / LAB;
                iface.n.y = -(xB - xA) / LAB;
                iface.n.z = 0.0;  // 2D plane
                iface.t2 = Vector3(0.0, 0.0, 1.0);
                iface.t1 = cross(iface.n, iface.t2);
                // Length in the XY-plane.
                iface.length = LAB;
                // Mid-point and area.
                iface.Ybar = 0.5 * (yA + yB);
                if ( myConfig.axisymmetric ) {
                    // Interface area per radian.
                    iface.area = LAB * iface.Ybar;
                } else {
                    // Assume unit depth in the Z-direction.
                    iface.area = LAB;
                }
                iface.pos = (getVtx(i,j).pos + getVtx(i,j+1).pos)/2.0;
            } // j loop
        } // i loop

        // North-facing interfaces.
        for (i = imin; i <= imax; ++i) {
            for (j = jmin; j <= jmax+1; ++j) {
                iface = getIfj(i,j);
                // These are the corners.
                xB = getVtx(i+1,j).pos.x;
                yB = getVtx(i+1,j).pos.y;
                xC = getVtx(i,j).pos.x;
                yC = getVtx(i,j).pos.y;
                LBC = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB));
                if (LBC < 1.0e-9) {
                    writefln("Zero length ifj[%d,%d]: %e", i, j, LBC);
                }
                // Direction cosines for the unit normal.
                iface.n.x = (yC - yB) / LBC;
                iface.n.y = -(xC - xB) / LBC;
                iface.n.z = 0.0;  // 2D plane
                iface.t2 = Vector3(0.0, 0.0, 1.0);
                iface.t1 = cross(iface.n, iface.t2);
                // Length in the XY-plane.
                iface.length = LBC;
                // Mid-point and area.
                iface.Ybar = 0.5 * (yC + yB);
                if ( myConfig.axisymmetric ) {
                    // Interface area per radian.
                    iface.area = LBC * iface.Ybar;
                } else {
                    // Assume unit depth in the Z-direction.
                    iface.area = LBC;
                }
                iface.pos = (getVtx(i+1,j).pos + getVtx(i,j).pos)/2.0;
            } // j loop
        } // i loop
    } // end calc_faces_2D()

    void calcVolumes3D()
    {
        size_t i, j, k;
        SolidFVCell cell;
        number volume;
        Vector3 centroid, p0, p1, p2, p3, p4, p5, p6, p7;
        //number iLength, jLength, kLength;
        for ( i = imin; i <= imax; ++i ) {
            for ( j = jmin; j <= jmax; ++j ) {
                for ( k = kmin; k <= kmax; ++k ) {
                    cell = getCell(i,j,k);
                    p0.set(cell.vtx[0].pos.x, cell.vtx[0].pos.y, cell.vtx[0].pos.z);
                    p1.set(cell.vtx[1].pos.x, cell.vtx[1].pos.y, cell.vtx[1].pos.z);
                    p2.set(cell.vtx[2].pos.x, cell.vtx[2].pos.y, cell.vtx[2].pos.z);
                    p3.set(cell.vtx[3].pos.x, cell.vtx[3].pos.y, cell.vtx[3].pos.z);
                    p4.set(cell.vtx[4].pos.x, cell.vtx[4].pos.y, cell.vtx[4].pos.z);
                    p5.set(cell.vtx[5].pos.x, cell.vtx[5].pos.y, cell.vtx[5].pos.z);
                    p6.set(cell.vtx[6].pos.x, cell.vtx[6].pos.y, cell.vtx[6].pos.z);
                    p7.set(cell.vtx[7].pos.x, cell.vtx[7].pos.y, cell.vtx[7].pos.z);
                    // Estimate the centroid so that we can use it as the peak
                    // of each of the pyramid sub-volumes.
                    centroid.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    // Mid-points of faces.
                    Vector3 pmN;
                    pmN.set(0.25*(p3.x+p2.x+p6.x+p7.x),
                            0.25*(p3.y+p2.y+p6.y+p7.y),
                            0.25*(p3.z+p2.z+p6.z+p7.z));
                    Vector3 pmE;
                    pmE.set(0.25*(p1.x+p2.x+p6.x+p5.x),
                            0.25*(p1.y+p2.y+p6.y+p5.y),
                            0.25*(p1.z+p2.z+p6.z+p5.z));
                    Vector3 pmS;
                    pmS.set(0.25*(p0.x+p1.x+p5.x+p4.x),
                            0.25*(p0.y+p1.y+p5.y+p4.y),
                            0.25*(p0.z+p1.z+p5.z+p4.z));
                    Vector3 pmW;
                    pmW.set(0.25*(p0.x+p3.x+p7.x+p4.x),
                            0.25*(p0.y+p3.y+p7.y+p4.y),
                            0.25*(p0.z+p3.z+p7.z+p4.z));
                    Vector3 pmT;
                    pmT.set(0.25*(p4.x+p5.x+p6.x+p7.x),
                            0.25*(p4.y+p5.y+p6.y+p7.y),
                            0.25*(p4.z+p5.z+p6.z+p7.z));
                    Vector3 pmB;
                    pmB.set(0.25*(p0.x+p1.x+p2.x+p3.x),
                            0.25*(p0.y+p1.y+p2.y+p3.y),
                            0.25*(p0.z+p1.z+p2.z+p3.z));
                    // Lengths between mid-points of faces.
                    // Note that we are assuming that the hexahedron is not very skewed
                    // when we later use these values as the widths of the hex cell.
                    //number dx, dy, dz;
                    //dx = pmE.x - pmW.x; dy = pmE.y - pmW.y; dz = pmE.z - pmW.z;
                    //iLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                    //dx = pmN.x - pmS.x; dy = pmN.y - pmS.y; dz = pmN.z - pmS.z;
                    //jLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                    //dx = pmT.x - pmB.x; dy = pmT.y - pmB.y; dz = pmT.z - pmB.z;
                    //kLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                    // writeln("Single hexahedron divided into six tetragonal dipyramids.");
                    // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
                    // Base of each dipyramid is specified clockwise from the outside.
                    number sub_volume; Vector3 sub_centroid;
                    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
                    pyramid_properties(p6, p7, p3, p2, centroid, false, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p5, p6, p2, p1, centroid, false, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p4, p5, p1, p0, centroid, false, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p7, p4, p0, p3, centroid, false, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p7, p6, p5, p4, centroid, false, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    pyramid_properties(p0, p1, p2, p3, centroid, false, sub_centroid, sub_volume);
                    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
                    //
                    if (volume < 0.0) {
                        throw new Error(text("Negative cell volume: solid_block ", id,
                                             " vol[", i, " ,", j, "]= ", volume));
                    }
                    //if (vol > max_vol) max_vol = vol;
                    //if (vol < min_vol) min_vol = vol;
                    cell.volume = volume;
                    //cell.areaxy = xyarea;
                } // j loop
            } // i loop
        } // k loop
    } // end calc_volumes_3D()

    void calcFaces3D()
    {
        SolidFVInterface iface;
        size_t i, j, k;
        Vector3 centroid, n, t1, t2, p0, p1, p2, p3;
        number area;
        // ifi interfaces are west interfaces, with their unit normal pointing east.
        for (k = kmin; k <= kmax; ++k) {
            for (j = jmin; j <= jmax; ++j) {
                for (i = imin; i <= imax+1; ++i) {
                    iface = getIfi(i,j,k);
                    // These are the corners.
                    p0.set(getVtx(i,j,k).pos.x, getVtx(i,j,k).pos.y, getVtx(i,j,k).pos.z);
                    p1.set(getVtx(i,j+1,k).pos.x, getVtx(i,j+1,k).pos.y, getVtx(i,j+1,k).pos.z);
                    p2.set(getVtx(i,j+1,k+1).pos.x, getVtx(i,j+1,k+1).pos.y, getVtx(i,j+1,k+1).pos.z);
                    p3.set(getVtx(i,j,k+1).pos.x, getVtx(i,j,k+1).pos.y, getVtx(i,j,k+1).pos.z);
                    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

                    // set properties.
                    iface.n.x = n.x;
                    iface.n.y = n.y;
                    iface.n.z = n.z;
                    iface.t2.x = t2.x;
                    iface.t2.y = t2.y;
                    iface.t2.z = t2.z;
                    iface.t1.x = t1.x;
                    iface.t1.y = t1.y;
                    iface.t1.z = t1.z;
                    iface.area = area;
                    iface.pos = (p0 + p1 + p2 + p3)/4.0;
                } // j loop
            } // j loop
        } // k loop
        // ifj interfaces are south interfaces, with their unit normal pointing north.
        for (k = kmin; k <= kmax; ++k) {
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax+1; ++j) {
                    iface = getIfj(i,j,k);
                    // These are the corners.
                    p0.set(getVtx(i,j,k).pos.x, getVtx(i,j,k).pos.y, getVtx(i,j,k).pos.z);
                    p1.set(getVtx(i,j,k+1).pos.x, getVtx(i,j,k+1).pos.y, getVtx(i,j,k+1).pos.z);
                    p2.set(getVtx(i+1,j,k+1).pos.x, getVtx(i+1,j,k+1).pos.y, getVtx(i+1,j,k+1).pos.z);
                    p3.set(getVtx(i+1,j,k).pos.x, getVtx(i+1,j,k).pos.y, getVtx(i+1,j,k).pos.z);
                    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

                    // set properties.
                    iface.n.x = n.x;
                    iface.n.y = n.y;
                    iface.n.z = n.z;
                    iface.t2.x = t2.x;
                    iface.t2.y = t2.y;
                    iface.t2.z = t2.z;
                    iface.t1.x = t1.x;
                    iface.t1.y = t1.y;
                    iface.t1.z = t1.z;
                    iface.area = area;
                    iface.pos = (p0 + p1 + p2 + p3)/4.0;
                } // j loop
            } // i loop
        } // k loop
        // ifk interfaces are bottom interfaces, with unit normal pointing to top.
        for (i = imin; i <= imax; ++i) {
            for (j = jmin; j <= jmax; ++j) {
                for (k = kmin; k <= kmax+1; ++k) {
                    iface = getIfk(i,j,k);
                    // These are the corners.
                    p0.set(getVtx(i,j,k).pos.x, getVtx(i,j,k).pos.y, getVtx(i,j,k).pos.z);
                    p1.set(getVtx(i+1,j,k).pos.x, getVtx(i+1,j,k).pos.y, getVtx(i+1,j,k).pos.z);
                    p2.set(getVtx(i+1,j+1,k).pos.x, getVtx(i+1,j+1,k).pos.y, getVtx(i+1,j+1,k).pos.z);
                    p3.set(getVtx(i,j+1,k).pos.x, getVtx(i,j+1,k).pos.y, getVtx(i,j+1,k).pos.z);
                    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

                    // set properties.
                    iface.n.x = n.x;
                    iface.n.y = n.y;
                    iface.n.z = n.z;
                    iface.t2.x = t2.x;
                    iface.t2.y = t2.y;
                    iface.t2.z = t2.z;
                    iface.t1.x = t1.x;
                    iface.t1.y = t1.y;
                    iface.t1.z = t1.z;
                    iface.area = area;
                    iface.pos = (p0 + p1 + p2 + p3)/4.0;
                } // j loop
            } // i loop
        } // k loop
    } // end calc_faces_2D()

    void setupSpatialDerivativeCalc()
    {
        int dim = myConfig.dimensions;
        foreach (c; cells) {
            gradients_T_lsq_setup(c, dim);
        }
    }

    override void assignCellLocationsForDerivCalc()
    {
        // This locations is only valid for the weighted least squares calculation.
        foreach (c; cells) {
            // First cell in the cloud is the cell itself.  Differences are taken about it.
            c.cloud_pos ~= c.pos;
            c.cloud_T ~= &(c.T);
            // Subsequent cells are the surrounding interfaces.
            foreach (i, f; c.iface) {
                c.cloud_pos ~= f.pos;
                c.cloud_T ~= &(f.T);
            } // end foreach face
        }
    }

    override void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel)
    {
        foreach(boundary; bc) { boundary.applyPreSpatialDerivActionAtBndryFaces(t, tLevel); }
    }

    override void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel, SolidFVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.whichBoundary == f.bc_id) { boundary.applyPreSpatialDerivActionAtBndryFaces(t, tLevel, f); }
        }
    }

    override void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel)
    {
        foreach(boundary; bc) { boundary.applyPreSpatialDerivActionAtBndryCells(t, tLevel); }
    }

    override void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel, SolidFVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.whichBoundary == f.bc_id) { boundary.applyPreSpatialDerivActionAtBndryCells(t, tLevel, f); }
        }
    }

    override void applyPostFluxAction(double t, int tLevel)
    {
        foreach(boundary; bc) { boundary.applyPostFluxAction(t, tLevel); }
    }

    override void applyPostFluxAction(double t, int tLevel, SolidFVInterface f)
    {
        foreach(boundary; bc) {
            if (boundary.whichBoundary == f.bc_id) { boundary.applyPostFluxAction(t, tLevel, f); }
        }
    }

    @nogc
    override void computeSpatialDerivatives(int ftl)
    {
        int dim = myConfig.dimensions;
        foreach (c; cells) {
            gradients_T_lsq(c, dim);
        }
    }

    @nogc
    override void averageTemperatures() {
        foreach (f; faces) {
            f.averageTemperature();
        }
    }

    @nogc
    override void averageProperties() {
        foreach (f; faces) {
            f.averageProperties();
        }
    }

    @nogc
    override void averageTGradients() {
        foreach (f; faces) {
            f.averageTGradient();
        }
        if (myConfig.solid_domain_augmented_deriv_avg) {
            foreach (f; faces) {
                f.augmentTGradient();
            }
        }
    }

    @nogc
    override void computeFluxes() {
        foreach (f; faces) {
            bool flux_already_set = f.is_on_boundary && bc[f.bc_id].setsFluxDirectly;
            if (!flux_already_set) {
                f.computeFlux(myConfig.dimensions, myConfig.solid_has_isotropic_properties);
            }
        }
    }

    override void clearSources()
    {
        foreach ( cell; cells ) cell.Q = 0.0;
    }

    @nogc
    override double determine_time_step_size(double cfl_value)
    {
        number alpha, L_min, dt_local, dt_allow;
        bool first = true;
        foreach(cell; cells) {

            if (GlobalConfig.dimensions == 2) {
                L_min = double.max;
                foreach (iface; cell.iface) { L_min = fmin(L_min, iface.length.re); }
            } else { // dimensions == 3
                number dx, dy, dz;
                auto nface = cell.iface[Face.north]; auto eface = cell.iface[Face.east];
                auto sface = cell.iface[Face.south]; auto wface = cell.iface[Face.west];
                auto tface = cell.iface[Face.top]; auto bface = cell.iface[Face.bottom];
                dx = eface.pos.x - wface.pos.x; dy = eface.pos.y - wface.pos.y; dz = eface.pos.z - wface.pos.z;
                auto iLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                dx = nface.pos.x - sface.pos.x; dy = nface.pos.y - sface.pos.y; dz = nface.pos.z - sface.pos.z;
                auto jLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                dx = tface.pos.x - bface.pos.x; dy = tface.pos.y - bface.pos.y; dz = tface.pos.z - bface.pos.z;
                auto kLen = sqrt(dx^^2 + dy^^2 + dz^^2);
                L_min = min(iLen, jLen, kLen);
            }
            stm.updateProperties(cell.ss);
            alpha = cell.ss.k / (cell.ss.rho*cell.ss.Cp);
            // Computational Fluid Mechanics and Heat Transfer
            // J. C. Tannehill, D. A. Anderson, and R. H Pletcher
            // Second Edition, Taylor & Francis, 1997
            //
            // The above reference defines the stable time-step for the heat equation (on pg.126 eq. 4.75) as
            // dt = r * (dx^2/alpha) where 0 <= r <= 0.5 for an explicit one-step scheme
            // since the CFL value set by the user may be used in both the fluid and solid domain, we absorb the
            // factor of 1/2 into the time-step calculation below. In effect it means that if a user sets
            // the cfl_value to be 1, to satisfy a stable explicit one-step scheme applied to the wave equation,
            // then that will equate to a cfl of 0.5 for the time-step used here in the heat equation
            dt_local = cfl_value * (L_min^^2) / (2.0*alpha);
            if (first) {
                dt_allow = dt_local;
                first = false;
            } else {
                dt_allow = fmin(dt_allow, dt_local);
            }
        }
        return dt_allow.re;
    } // end determine_time_step_size()

    void initialize_jacobian(int spatial_order_of_jacobian, double eps, int fill_in)
    {
        // This routine initializes the Jacobian matrix attached to this SolidBlock object.
        if (spatial_order_of_jacobian < 0 || spatial_order_of_jacobian > 2) {
            throw new FlowSolverException("ERROR: invalid preconditioner approximation for solid domain, select either 0 1 2");
        }

        size_t nentry = 0;
        // gather the expected number of non-zero entries in the flow Jacobian
        foreach (cell; cells) {
            cell.gather_residual_stencil_lists(spatial_order_of_jacobian);
            nentry += cell.cell_list.length;
        }

        // TODO: gather ghost-cell stencils

        // pre-size the sparse matrix arrays associated with the Jacobian matrix
        jacobian_nonzero_pattern(spatial_order_of_jacobian, nentry, cells.length, eps, fill_in);
    } // end initialize_jacobian()

    void jacobian_nonzero_pattern(int spatial_order_of_jacobian, size_t nentry, size_t ncells, double sigma, int iluFill) {
        // This routine pre-sizes the sparse matrix with the appropriate non-zero structure

        jacobian = new FlowJacobian(sigma, myConfig.dimensions, 1, spatial_order_of_jacobian, nentry, ncells);
        jacobian.local.ia[jacobian.ia_idx] = 0; // the first entry will always be filled
        jacobian.ia_idx += 1;
        foreach (pcell; cells) {
            // loop through cells that will have non-zero entries
            foreach(cell; pcell.cell_list) {
                assert(!cell.is_ghost, "Oops, we expect to not find a ghost cell at this point.");
                size_t jidx = cell.id; // column index
                // populate entry with a place holder value
                jacobian.local.aa[jacobian.aa_idx] = 1.0;
                // fill out the sparse matrix ready for the next entry in the row
                jacobian.aa_idx += 1;
                jacobian.local.ja[jacobian.ja_idx] = jidx;
                jacobian.ja_idx += 1;
            }
            // prepare the sparse matrix for a new row
            jacobian.local.ia[jacobian.ia_idx] = jacobian.aa_idx;
            jacobian.ia_idx += 1;
        }

        // For ILU(p>0) we need to modify the sparsity pattern
        if (iluFill > 0) {

            // [TODO] think about making the lev matrix sparse as well KAD 2022-03-31
            // construct a level matrix
            int p = iluFill;
            int n = to!int(jacobian.local.ia.length)-1;
            int[][] lev; // fill levels
            lev.length = n;
            foreach ( i; 0..n) lev[i].length = n;
            // assign initial fill levels
            foreach ( i; 0 .. n ) {
                foreach ( j; 0 .. n ) {
                    if (jacobian.local[i,j] < 1.0) lev[i][j] = n-1;
                    else lev[i][j] = 0;
                }
            }

            // apply symbolic phase algorithm
            foreach ( i; 1 .. n ) { // Begin from 2nd row
                foreach ( k; 0 .. i ) {
                    if (lev[i][k] <= p) {
                        foreach ( j ; k+1..n) {
                            lev[i][j] = min(lev[i][j], lev[i][k]+lev[k][j]+1);
                        }
                    }
                }
            }

            // modify the sparse matrix non-zero pattern based on the level matrix
            foreach ( i; 0..n) {
                foreach ( j; 0..n) {
                    // add new entry
                    if (lev[i][j] <= p) { jacobian.local[i,j] = 1.0; }
                }
            }
        }

    } // end jacobian_nonzero_pattern()

    void evaluate_jacobian()
    {
        // Higher level routine used to evaluate the Jacobian attached to this SolidBlock object.

        // zero entries
        jacobian.local.aa[] = 0.0;
        version(complex_numbers) {
            // do nothing
        } else {
            // the real-valued finite difference needs a base residual (R0)
            foreach(cell; cells) { evalRHS(myConfig.dimensions, 0, 0, cell.cell_list, cell.face_list, cell); }
        }

        // construct the Jacobian
        foreach (cell; cells) { evaluate_cell_contribution_to_jacobian(cell); }

        // TODO: add boundary condition corrections to boundary cells

    } // end evaluate_jacobian()

    void evaluate_cell_contribution_to_jacobian(SolidFVCell pcell)
    {
        int dim = myConfig.dimensions;
        int ftl = to!int(myConfig.n_flow_time_levels-1);
        int gtl = 0;
        auto eps0 = jacobian.eps;
        number eps;

        // save a copy of the original cell state
        number T_save = pcell.T;
        number e_save = pcell.e[0];

        // perturb the current cell conserved quantities and then evaluate
        // the residuals for each cell in the local domain of influence
        version(complex_numbers) { eps = complex(0.0, eps0.re); }
        else { eps = eps0*fabs(pcell.e[0]) + eps0; }
        pcell.e[ftl] = pcell.e[0];
        pcell.e[ftl] += eps;
        pcell.ss.e = pcell.e[ftl];
        stm.updateTemperature(pcell.ss);
        pcell.T = pcell.ss.T;

        // evaluate perturbed residuals in local stencil
        evalRHS(dim, gtl, ftl, pcell.cell_list, pcell.face_list, pcell);

        // fill local Jacobian entries
        foreach (cell; pcell.cell_list) {
            version(complex_numbers) {
                cell.dRde = cell.dedt[ftl].im/eps.im;
            } else {
                cell.dRde = (cell.dedt[ftl]-cell.dedt[0])/eps;
            }
        }

        // clear imaginary components
        version (complex_numbers) {
            foreach (cell; pcell.cell_list) {
                cell.clear_imaginary_components();
            }
        } else {
            pcell.T = T_save;
            pcell.e[0] = e_save;
            evalRHS(dim, gtl, ftl, pcell.cell_list, pcell.face_list, pcell);
        }

        // we now populate the pre-sized sparse matrix representation of the Jacobian
        foreach(cell; pcell.cell_list) {
            assert(!cell.is_ghost, "Oops, we expect to not find a ghost cell at this point.");
            size_t I = cell.id;
            size_t J = pcell.id;
            jacobian.local[I,J] = cell.dRde.re;
        }
    } // end evaluate_cell_contribution_to_jacobian()

    void evalRHS(int dim, int gtl, int ftl, ref SolidFVCell[] cell_list, SolidFVInterface[] iface_list, SolidFVCell pcell)
    {
        // This routine evaluates the RHS residual on a subset of cells in the SolidBlock.
        // It is used when constructing the numerical Jacobian.
        // Its effect should replicate evalRHS() in solid_loose_coupling_update.d.

        // average T at cell interfaces
        foreach (f; iface_list) {
            f.averageTemperature();
        }

        // apply BCs
        foreach(f; iface_list) {
            if (f.is_on_boundary) { applyPreSpatialDerivActionAtBndryFaces(0.0, 0, f); }
        }

        // compute lsq gradient
        foreach(cell; cell_list) {
            gradients_T_lsq(cell, dim);
        }

        // evaluate fluxes
        foreach (f; iface_list) {
            f.averageTGradient();
            if (myConfig.solid_domain_augmented_deriv_avg) {
                f.augmentTGradient();
            }
            bool flux_already_set = f.is_on_boundary && bc[f.bc_id].setsFluxDirectly;
            if (!flux_already_set) {
                f.computeFlux(myConfig.dimensions, myConfig.solid_has_isotropic_properties);
            }
        }

        // apply BCs
        foreach(f; iface_list) {
            if (f.is_on_boundary) { applyPostFluxAction(0.0, 0, f); }
        }

        // sum fluxes
        foreach(cell; cell_list) {
            cell.timeDerivatives(ftl, dim);
        }

    } // end evalRHS()
}

@nogc
void gradients_T_lsq_setup(SolidFVCell c, int dimensions)
{
    size_t n = c.cloud_pos.length;
    number[12] weights2;
    size_t loop_init;
    loop_init = 1; // All points count.
    weights2[0] = 0.0; // and doesn't enter into the sum itself.
    //
    // Calculate weights used in the least-squares gradient calculation.
    // These are the square of the weights on the original linear constraint eqns.
    // These weights are calculated with the current interface/vertex
    // as the reference point (pos).
    // For the "faces" spatial location we are expecting the primary point
    // (i.e. the face at which we are calculating the gradients) to be in
    // the first cloud position.
    number x0 = c.pos.x; number y0 = c.pos.y; number z0 = c.pos.z;
    if (dimensions == 2) {
        foreach (i; loop_init .. n) {
            number dx = c.cloud_pos[i].x - x0;
            number dy = c.cloud_pos[i].y - y0;
            weights2[i] = 1.0/(dx*dx+dy*dy);
        }
    } else { //3D
        foreach (i; loop_init .. n) {
            number dx = c.cloud_pos[i].x - x0;
            number dy = c.cloud_pos[i].y - y0;
            number dz = c.cloud_pos[i].z - z0;
            weights2[i] = 1.0/(dx*dx+dy*dy+dz*dz);
        }
    }

    number[12] dx, dy, dz;
    //
    // Assemble and invert the normal matrix.
    // We'll reuse the resulting inverse for each flow-field quantity.
    number[12] wx, wy, wz;
    if (dimensions == 3) {
        number[6][3] xTx; // normal matrix, augmented to give 6 entries per row
        number xx = 0.0; number xy = 0.0; number xz = 0.0;
        number yy = 0.0; number yz = 0.0; number zz = 0.0;
        foreach (i; loop_init .. n) {
            dx[i] = c.cloud_pos[i].x - x0;
            dy[i] = c.cloud_pos[i].y - y0;
            dz[i] = c.cloud_pos[i].z - z0;
            xx += weights2[i]*dx[i]*dx[i];
            xy += weights2[i]*dx[i]*dy[i];
            xz += weights2[i]*dx[i]*dz[i];
            yy += weights2[i]*dy[i]*dy[i];
            yz += weights2[i]*dy[i]*dz[i];
            zz += weights2[i]*dz[i]*dz[i];
        }
        xTx[0][0] = xx; xTx[0][1] = xy; xTx[0][2] = xz;
        xTx[1][0] = xy; xTx[1][1] = yy; xTx[1][2] = yz;
        xTx[2][0] = xz; xTx[2][1] = yz; xTx[2][2] = zz;
        xTx[0][3] = 1.0; xTx[0][4] = 0.0; xTx[0][5] = 0.0;
        xTx[1][3] = 0.0; xTx[1][4] = 1.0; xTx[1][5] = 0.0;
        xTx[2][3] = 0.0; xTx[2][4] = 0.0; xTx[2][5] = 1.0;
        double very_small_value = 1.0e-16*(normInf!(3,3,6,number)(xTx).re)^^3;
        if (0 != computeInverse!(3,3,6,number)(xTx, very_small_value)) {
            throw new FlowSolverException("Failed to invert LSQ normal matrix");
            // Assume that the rows are linearly dependent
            // because the sample points are coplanar or colinear.
        }
        // Prepare final weights for later use in the reconstruction phase.
        foreach (i; loop_init .. n) {
            c.wx[i] = xTx[0][3]*dx[i] + xTx[0][4]*dy[i] + xTx[0][5]*dz[i];
            c.wx[i] *= weights2[i];
            c.wy[i] = xTx[1][3]*dx[i] + xTx[1][4]*dy[i] + xTx[1][5]*dz[i];
            c.wy[i] *= weights2[i];
            c.wz[i] = xTx[2][3]*dx[i] + xTx[2][4]*dy[i] + xTx[2][5]*dz[i];
            c.wz[i] *= weights2[i];
        }
    } else {
        // dimensions == 2
        number[4][2] xTx; // normal matrix, augmented to give 4 entries per row
        number xx = 0.0; number xy = 0.0; number yy = 0.0;
        foreach (i; loop_init .. n) {
            dx[i] = c.cloud_pos[i].x - x0;
            dy[i] = c.cloud_pos[i].y - y0;
            xx += weights2[i]*dx[i]*dx[i];
            xy += weights2[i]*dx[i]*dy[i];
            yy += weights2[i]*dy[i]*dy[i];
        }
        xTx[0][0] = xx; xTx[0][1] = xy;
        xTx[1][0] = xy; xTx[1][1] = yy;
        xTx[0][2] = 1.0; xTx[0][3] = 0.0;
        xTx[1][2] = 0.0; xTx[1][3] = 1.0;
        double very_small_value = 1.0e-16*(normInf!(2,2,4,number)(xTx).re)^^2;
        if (0 != computeInverse!(2,2,4,number)(xTx, very_small_value)) {
            throw new FlowSolverException("Failed to invert LSQ normal matrix");
            // Assume that the rows are linearly dependent
            // because the sample points are colinear.
        }
        //number[12] wx, wy, wz;
        // Prepare final weights for later use in the reconstruction phase.
        foreach (i; loop_init .. n) {
            c.wx[i] = xTx[0][2]*dx[i] + xTx[0][3]*dy[i];
            c.wx[i] *= weights2[i];
            c.wy[i] = xTx[1][2]*dx[i] + xTx[1][3]*dy[i];
            c.wy[i] *= weights2[i];
            c.wz[i] = 0.0;
        }
    }
}

@nogc
void gradients_T_lsq(SolidFVCell c, int dimensions)
{
    size_t n = c.cloud_pos.length;
    size_t loop_init;
    loop_init = 1; // All points count.

    number T0;
    number[3] gradT;
    T0 = *(c.cloud_T[0]);
    gradT[0] = 0.0; gradT[1] = 0.0; gradT[2] = 0.0;
    foreach (i; loop_init .. n) {
        number dT = *(c.cloud_T[i]) - T0;
        gradT[0] += c.wx[i] * dT;
        gradT[1] += c.wy[i] * dT;
        if (dimensions == 3) { gradT[2] += c.wz[i] * dT; }
    }

    c.dTdx = gradT[0];
    c.dTdy = gradT[1];
    if (dimensions == 3) c.dTdz = gradT[2];
    else { c.dTdz = 0.0; }
}
