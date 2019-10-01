// sfluidblock.d
// Class for structured blocks of cells, for use within Eilmer4.
// This is the "classic" block within the mbcns/Eilmer series 
// of flow simulation codes.

// Peter J. 2014-07-20 first cut.

module sfluidblock;

import std.conv;
import std.file;
import std.json;
import std.stdio;
import std.format;
import std.string;
import std.array;
import std.algorithm;
import std.math;
import nm.complex;
import nm.number;

import util.lua;
import json_helper;
import lua_helper;
import gzip;
import geom;
import gas;
import kinetics;
import globalconfig;
import globaldata;
import flowstate;
import fluxcalc;
import flowgradients;
import fvcore;
import fvvertex;
import fvinterface;
import fvcell;
import onedinterp;
import fluidblock;
import bc;
import grid_motion;

// Ghost-cell layers surround the active cells of a block.
// For the high-order reconstruction right to the edge on
// structured-grids, we will need a minimum number of ghost cells.
version(nghost3) {
    // For Lachlan's higher-order work, 2019.
    enum n_ghost_cell_layers = 3;
} else {
    // "Traditional" value (since 1991).
    enum n_ghost_cell_layers = 2;
}

// EPSILON parameter for numerical differentiation of flux jacobian
// Value used based on Vanden and Orkwis (1996), AIAA J. 34:6 pp. 1125-1129
immutable double EPSILON = 1.0e-8;
immutable double ESSENTIALLY_ZERO = 1.0e-15;

class SFluidBlock: FluidBlock {
public:
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
    size_t[] micell, mjcell, mkcell; // locations of monitor cells

    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax; 
    size_t jmin, jmax;
    size_t kmin, kmax;
    
    // A place to store coordinates of the corner vertices.
    // For a moving-grid simulation these will be kept up to date
    // and communicated to user-defined Lua functions via infoFluidBlock.
    double[24] corner_coords;
    int[] inflow_partners;
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    OneDInterpolator one_d;
    
private:
    StructuredGrid grid; // for reading and writing
    // Total number of cells in each direction for this block.
    // these will be used in the array allocation routines.
    size_t _nidim;
    size_t _njdim;
    size_t _nkdim;
    // Most of the data is stored in the following arrays.
    // ctr = cell center values
    // ifi = east-facing face properties and fluxes (unit normal in the i-index direction)
    // ifj = north-facing face properties and fluxes (normal in the j-index direction)
    // ifk = top-facing 
    // vtx = cell vertex values (used for the viscous terms, mostly)
    // sifi, sifj and sifk are secondary-cell faces (whose corner nodes are the
    //                     the primary-cell centres.
    FVCell[] _ctr;
    FVInterface[] _ifi;
    FVInterface[] _ifj;
    FVInterface[] _ifk;
    FVVertex[] _vtx;
    FVInterface[] _sifi;
    FVInterface[] _sifj;
    FVInterface[] _sifk;

public:
    this(int blk_id, size_t nicell, size_t njcell, size_t nkcell, string label)
    {
        super(blk_id, Grid_t.structured_grid, nicell*njcell*nkcell, label);
        this.nicell = nicell;
        this.njcell = njcell;
        this.nkcell = nkcell;
        // Fill in other data sizes.
        _nidim = nicell + 2 * n_ghost_cell_layers;
        _njdim = njcell + 2 * n_ghost_cell_layers;
        // Indices, in each grid direction for the active cells.
        // These limits are inclusive. The mincell and max cell
        // are both within the active set of cells.
        imin = n_ghost_cell_layers; imax = imin + nicell - 1;
        jmin = n_ghost_cell_layers; jmax = jmin + njcell - 1;
        if ( GlobalConfig.dimensions == 2 ) {
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
        // Workspace for flux_calc method.
        one_d = new OneDInterpolator(dedicatedConfig[id]);

    } // end constructor

    this(int blk_id, JSONValue json_data)
    {
        nicell = getJSONint(json_data, "nic", 0);
        njcell = getJSONint(json_data, "njc", 0);
        nkcell = getJSONint(json_data, "nkc", 0);
        label = getJSONstring(json_data, "label", "");
        this(blk_id, nicell, njcell, nkcell, label);
        active = getJSONbool(json_data, "active", true);
        omegaz = getJSONdouble(json_data, "omegaz", 0.0);
    } // end constructor from json

    @nogc override int get_interpolation_order()
    {
        return one_d.get_interpolation_order();
    }

    @nogc override void set_interpolation_order(int order)
    {
        one_d.set_interpolation_order(order);
    }

    override void init_lua_globals()
    {
        lua_pushinteger(myL, nicell); lua_setglobal(myL, "nicell");
        lua_pushinteger(myL, njcell); lua_setglobal(myL, "njcell");
        lua_pushinteger(myL, nkcell); lua_setglobal(myL, "nkcell");
        lua_pushinteger(myL, Face.north); lua_setglobal(myL, "north");
        lua_pushinteger(myL, Face.east); lua_setglobal(myL, "east");
        lua_pushinteger(myL, Face.south); lua_setglobal(myL, "south");
        lua_pushinteger(myL, Face.west); lua_setglobal(myL, "west");
        lua_pushinteger(myL, Face.top); lua_setglobal(myL, "top");
        lua_pushinteger(myL, Face.bottom); lua_setglobal(myL, "bottom");
        lua_pushinteger(myL, n_ghost_cell_layers); lua_setglobal(myL, "n_ghost_cell_layers");
        setSampleHelperFunctions(myL);
        // Note that the sampleFluidCell function can be expected to work only in serial mode.
        // Once it is called from a thread, other than the main thread, it may not
        // have access to properly initialized data for any other block.
        setGridMotionHelperFunctions(myL);
    } // end init_lua_globals()

    override void init_boundary_conditions(JSONValue json_data)
    // Initialize boundary conditions after the blocks are fully constructed,
    // because we want access to the full collection of valid block references.
    {
        foreach (boundary; 0 .. (myConfig.dimensions == 3 ? 6 : 4)) {
            string json_key = "boundary_" ~ face_name[boundary];
            auto bc_json_data = json_data[json_key];
            bc ~= make_BC_from_json(bc_json_data, id, boundary);
        }
        foreach (bci; bc) bci.post_bc_construction();
    } // end init_boundary_conditions()

    override string toString() const
    {
        char[] repr;
        repr ~= "SFluidBlock(";
        repr ~= "id=" ~ to!string(id);
        repr ~= ", label=\"" ~ label ~ "\"";
        repr ~= ", active=" ~ to!string(active);
        repr ~= ", grid_type=\"" ~ gridTypeName(grid_type) ~ "\"";
        repr ~= ", omegaz=" ~ to!string(omegaz);
        repr ~= ", nicell=" ~ to!string(nicell);
        repr ~= ", njcell=" ~ to!string(njcell);
        repr ~= ", nkcell=" ~ to!string(nkcell);
        repr ~= ", n_ghost_cell_layers=" ~ to!string(n_ghost_cell_layers);
        repr ~= ", \n    bc=["~ face_name[0] ~ "=" ~ to!string(bc[0]);
        foreach (i; 1 .. (myConfig.dimensions == 3 ? 6 : 4)) {
            repr ~= ",\n        " ~ face_name[i] ~ "=" ~ to!string(bc[i]);
        }
        repr ~= "\n       ]"; // end bc list
        repr ~= ")";
        return to!string(repr);
    } // end toString()

    @nogc size_t to_single_index()(size_t i, size_t j, size_t k) const
    {
        assert(i < _nidim && j < _njdim && k < _nkdim, "Index out of bounds.");
        return (k*_njdim + j)*_nidim + i; 
    }

    @nogc size_t[3] to_ijk_indices(size_t gid) const
    {
        size_t[3] ijk;
        size_t slabDim = _njdim * _nidim;
        size_t k = gid / slabDim;
        size_t j = (gid - k*slabDim) / _nidim;
        size_t i = gid - k*slabDim - j*_nidim;
        ijk[0] = i; ijk[1] = j; ijk[2] = k;
        return ijk;
    }

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

    @nogc size_t ijk_indices_to_cell_id(size_t i, size_t j, size_t k=0) const
    // ijk indices into the underlying block array.
    {
        i -= n_ghost_cell_layers;
        j -= n_ghost_cell_layers;
        k = (myConfig.dimensions == 2) ? 0 : k - n_ghost_cell_layers;
        return (k*njcell + j)*nicell + i;
    }
    
    @nogc size_t[3] cell_id_to_ijk_indices(size_t id) const
    {
        size_t[3] ijk;
        // Invert the cell_id calculation given above.
        size_t k = (myConfig.dimensions == 2) ? 0 : id/(njcell*nicell);
        size_t j = (id - k*(njcell*nicell))/nicell;
        size_t i = id - k*(njcell*nicell) - j*nicell;
        // Add on the ghost-cell offset.
        ijk[0] = i + n_ghost_cell_layers;
        ijk[1] = j + n_ghost_cell_layers;
        ijk[2] = (myConfig.dimensions == 2) ? 0 : k + n_ghost_cell_layers;
        return ijk;
    } // end cell_id_to_ijk_indices()

    @nogc ref FVCell get_cell()(size_t i, size_t j, size_t k=0) 
    {
        return _ctr[to_single_index!()(i,j,k)];
    }
    @nogc ref FVInterface get_ifi()(size_t i, size_t j, size_t k=0) 
    {
        return _ifi[to_single_index!()(i,j,k)];
    }
    @nogc ref FVInterface get_ifj()(size_t i, size_t j, size_t k=0)
    {
        return _ifj[to_single_index!()(i,j,k)];
    }
    @nogc ref FVInterface get_ifk()(size_t i, size_t j, size_t k=0)
    {
        return _ifk[to_single_index!()(i,j,k)];
    }
    @nogc ref FVVertex get_vtx()(size_t i, size_t j, size_t k=0)
    {
        return _vtx[to_single_index!()(i,j,k)];
    }

    @nogc
    override void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
        grid.find_enclosing_cell(p, indx, found); // delegate to the grid object
    }

    override void init_grid_and_flow_arrays(string gridFileName)
    {
        assemble_arrays();
        bind_interfaces_vertices_and_cells();
        store_references_for_derivative_calc(0);
        if (myConfig.verbosity_level > 1) { writeln("init_grid_and_flow_arrays(): Start block ", id); }
        grid = new StructuredGrid(gridFileName, myConfig.grid_format);
        grid.sort_cells_into_bins();
        sync_vertices_from_underlying_grid(0);
        // Set references to boundary faces in bc objects.
        // north boundary
        foreach (k; kmin .. kmax+1) {
            foreach (i; imin .. imax+1) {
                bc[Face.north].faces ~= get_ifj!()(i, jmax+1, k);
                bc[Face.north].outsigns ~= 1;
            }
        }
        foreach (k; kmin .. kmax+1) {
            foreach (j; jmin .. jmax+1) {
                bc[Face.east].faces ~= get_ifi!()(imax+1, j, k);
                bc[Face.east].outsigns ~= 1;
            }
        }
        foreach (k; kmin .. kmax+1) {
            foreach (i; imin .. imax+1) {
                bc[Face.south].faces ~= get_ifj!()(i, jmin, k);
                bc[Face.south].outsigns ~= -1;
            }
        }
        foreach (k; kmin .. kmax+1) {
            foreach (j; jmin .. jmax+1) {
                bc[Face.west].faces ~= get_ifi!()(imin, j, k);
                bc[Face.west].outsigns ~= -1;
            }
        }
        if (myConfig.dimensions == 3) {
            foreach (j; jmin .. jmax+1) {
                foreach (i; imin .. imax+1) {
                    bc[Face.top].faces ~= get_ifk!()(i, j, kmax+1);
                    bc[Face.top].outsigns ~= 1;
                }
            }
            foreach (j; jmin .. jmax+1) {
                foreach (i; imin .. imax+1) {
                    bc[Face.bottom].faces ~= get_ifk!()(i, j, kmin);
                    bc[Face.bottom].outsigns ~= -1;
                }
            }
        } // end if dimensions == 3
        //
        // Set up the lists of indices for look-up of cells and faces
        // from a given vertex.
        cellIndexListPerVertex.length = vertices.length;
        foreach (i, c; cells) {
            foreach (vtx; c.vtx) { cellIndexListPerVertex[vtx.id] ~= i; }
        }
        faceIndexListPerVertex.length = vertices.length;
        foreach (i, f; faces) {
            foreach (vtx; f.vtx) { faceIndexListPerVertex[vtx.id] ~= i; }
        }
    } // end init_grid_and_flow_arrays()

    void assemble_arrays()
    // We shouldn't be calling this until the essential bits of the GlobalConfig
    // and the local myConfig instances have been set up.
    {
        if (myConfig.verbosity_level >= 2) {
            writefln("assemble_arrays(): Begin for block %d", id);
        }
        // Check for obvious errors.
        if ( _nidim <= 0 || _njdim <= 0 || _nkdim <= 0 ) {
            string msg = text("resize_arrays(): invalid dimensions nidim=",
                              _nidim, " njdim=", _njdim, " nkdim=", _nkdim);
            throw new FlowSolverException(msg);
        }
        size_t ntot = _nidim * _njdim * _nkdim;
        bool lsq_workspace_at_faces = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces);
        bool lsq_workspace_at_vertices = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.vertices);
        try {
            // Create the cell and interface objects for the entire structured block.
            // This includes the layer of surrounding ghost cells.
            // The for each cell, face and vertex, the global index (gid)
            // will be the index into the array held privately by this class.
            //
            // In the FluidBlock base class, we will keep an array of "active" cells
            // that may be accessed directly by other parts of the code.
            // Providing such access brings the structured-grid code a little closer
            // to the flavour of the unstructured-grid code.
            foreach (gid; 0 .. ntot) {
                _ctr ~= new FVCell(myConfig);
                _ifi ~= new FVInterface(myConfig, lsq_workspace_at_faces);
                _ifj ~= new FVInterface(myConfig, lsq_workspace_at_faces);
                if ( myConfig.dimensions == 3 ) {
                    _ifk ~= new FVInterface(myConfig, lsq_workspace_at_faces);
                }
                _vtx ~= new FVVertex(myConfig, lsq_workspace_at_vertices);
            }
            // Now, assemble the lists of references to the cells, vertices and faces
            // in standard order for a structured grid.
            // These arrays are held by the FluidBlock base class and allow us to handle
            // a structured-grid block much as we would an unstructured-grid block.
            if (myConfig.dimensions == 2) {
                foreach (j; jmin .. jmax+1) {
                    foreach (i; imin .. imax+1) { cells ~= get_cell!()(i, j); }
                }
                foreach (j; jmin .. jmax+2) {
                    foreach (i; imin .. imax+2) { vertices ~= get_vtx!()(i, j); }
                }
                foreach (j; jmin .. jmax+1) {
                    foreach (i; imin .. imax+2) { faces ~= get_ifi!()(i, j); }
                }
                foreach (j; jmin .. jmax+2) {
                    foreach (i; imin .. imax+1) { faces ~= get_ifj!()(i, j); }
                }
            } else { // assume 3D
                foreach (k; kmin .. kmax+1) {
                    foreach (j; jmin .. jmax+1) {
                        foreach (i; imin .. imax+1) { cells ~= get_cell!()(i, j, k); }
                    }
                }
                foreach (k; kmin .. kmax+2) {
                    foreach (j; jmin .. jmax+2) {
                        foreach (i; imin .. imax+2) { vertices ~= get_vtx!()(i, j, k); }
                    }
                }
                foreach (k; kmin .. kmax+1) {
                    foreach (j; jmin .. jmax+1) {
                        foreach (i; imin .. imax+2) { faces ~= get_ifi!()(i, j, k); }
                    }
                }
                foreach (k; kmin .. kmax+1) {
                    foreach (j; jmin .. jmax+2) {
                        foreach (i; imin .. imax+1) { faces ~= get_ifj!()(i, j, k); }
                    }
                }
                foreach (k; kmin .. kmax+2) {
                    foreach (j; jmin .. jmax+1) {
                        foreach (i; imin .. imax+1) { faces ~= get_ifk!()(i, j, k); }
                    }
                }
            } // end if dimensions
        } catch (Exception e) {
            writeln("Failed while assembling block arrays.");
            writefln("nicell=%d njcell=%d nkcell=%d", nicell, njcell, nkcell);
            writefln("nidim=%d njdim=%d nkdim=%d", _nidim, _njdim, _nkdim);
            writeln("Probably ran out of memory.");
            writeln("Be a little less ambitious and try a smaller grid next time.");
            writefln("System message: %s", e.msg);
            throw new FlowSolverException("Block.assemble_arrays() failed.");
        }
        //
        // Make the cell. vertex. and face.id consistent with the index in the array.
        // We will depend on this equality in other parts of the flow solver.
        // We also note that these cells are interior to the block (i.e. not ghost cells)
        foreach (i, c; cells) { c.id = to!int(i); c.contains_flow_data = true; c.is_interior_to_domain = true; }
        foreach (i, v; vertices) { v.id = to!int(i); }
        foreach (i, f; faces) { f.id = to!int(i); }
        // Alter the id values of the ghost cells to be a bit like those in the
        // unstructured-grid block.
        int cell_id = ghost_cell_start_id;
        int face_id = ghost_cell_start_id;
        int vtx_id = ghost_cell_start_id;
        if (myConfig.dimensions == 2) {
            foreach (j; 0 .. _njdim) {
                foreach (i; 0 .. _nidim) {
                    auto c = get_cell!()(i, j); if (c.id == -1) { c.id = cell_id; ++cell_id; }
                    auto f = get_ifi!()(i,j); if (f.id == -1) { f.id = face_id; ++face_id; }
                    f = get_ifj!()(i,j); if (f.id == -1) { f.id = face_id; ++face_id; }
                    auto v = get_vtx!()(i,j); if (v.id == -1) { v.id = vtx_id; ++vtx_id; }
                }
            }
        } else { // assume 3D
            foreach (k; 0 .. _nkdim) {
                foreach (j; 0 .. _njdim) {
                    foreach (i; 0 .. _nidim) {
                        auto c = get_cell!()(i, j, k); if (c.id == -1) { c.id = cell_id; ++cell_id; }
                        auto f = get_ifi!()(i,j,k); if (f.id == -1) { f.id = face_id; ++face_id; }
                        f = get_ifj!()(i,j,k); if (f.id == -1) { f.id = face_id; ++face_id; }
                        f = get_ifk!()(i,j,k); if (f.id == -1) { f.id = face_id; ++face_id; }
                        auto v = get_vtx!()(i,j,k); if (v.id == -1) { v.id = vtx_id; ++vtx_id; }
                    }
                }
            }
        } // end if dimensions
        if (myConfig.verbosity_level >= 2) {
            writefln("Done assembling arrays for %d cells in block %d.", ntot, id);
        }
    } // end of assemble_arrays()

    void bind_interfaces_vertices_and_cells()
    {
        // There is a fixed order of faces and vertices for each cell.
        // Refer to fvcore.d
        size_t kstart, kend;
        if (myConfig.dimensions == 3) {
            kstart = kmin - 1;
            kend = kmax + 1;
        } else {
            kstart = 0;
            kend = 0;
        }
        // With the ranges above and in the following nested loops,
        // we make connections for the first layer of ghost cells, also.
        for ( size_t k = kstart; k <= kend; ++k ) {
            for ( size_t j = jmin-1; j <= jmax+1; ++j ) {
                for ( size_t i = imin-1; i <= imax+1; ++i ) {
                    FVCell cell = get_cell!()(i,j,k);
                    cell.iface.length = 0; cell.outsign.length = 0;
                    cell.iface ~= get_ifj!()(i,j+1,k); cell.outsign ~= 1.0; // north
                    cell.iface ~= get_ifi!()(i+1,j,k); cell.outsign ~= 1.0; // east
                    cell.iface ~= get_ifj!()(i,j,k); cell.outsign ~= -1.0; // south
                    cell.iface ~= get_ifi!()(i,j,k); cell.outsign ~= -1.0; // west
                    cell.vtx.length = 0;
                    cell.vtx ~= get_vtx!()(i,j,k);
                    cell.vtx ~= get_vtx!()(i+1,j,k);
                    cell.vtx ~= get_vtx!()(i+1,j+1,k);
                    cell.vtx ~= get_vtx!()(i,j+1,k);
                    if (myConfig.dimensions == 3) {
                        cell.iface ~= get_ifk!()(i,j,k+1); cell.outsign ~= 1.0; // top
                        cell.iface ~= get_ifk!()(i,j,k); cell.outsign ~= -1.0; // bottom
                        cell.vtx ~= get_vtx!()(i,j,k+1);
                        cell.vtx ~= get_vtx!()(i+1,j,k+1);
                        cell.vtx ~= get_vtx!()(i+1,j+1,k+1);
                        cell.vtx ~= get_vtx!()(i,j+1,k+1);
                    } // end if
                } // for i
            } // for j
        } // for k
        //
        // Sometimes it is convenient for an interface to come complete 
        // with information about the vertices that define it and also
        // the cells that adjoin it.
        //
        // ifi interfaces are west interfaces, with their unit normal pointing east.
        // In 2D, vtx0==p00, vtx1==p01.
        // In 3D, the cycle [vtx0,vtx1,vtx2,vtx3] progresses counter-clockwise around 
        // the periphery of the face when the normal unit vector is pointing toward you.
        // t1 vector aligned with j-index direction
        // t2 vector aligned with k-index direction
        for (size_t k = kmin; k <= kmax; ++k) {
            for (size_t j = jmin; j <= jmax; ++j) {
                for (size_t i = imin; i <= imax+1; ++i) {
                    auto IFace = get_ifi!()(i,j,k);
                    IFace.vtx.length = 0;
                    if (myConfig.dimensions == 3) {
                        IFace.vtx ~= get_vtx!()(i,j,k);
                        IFace.vtx ~= get_vtx!()(i,j+1,k);
                        IFace.vtx ~= get_vtx!()(i,j+1,k+1);
                        IFace.vtx ~= get_vtx!()(i,j,k+1);
                    } else {
                        IFace.vtx ~= get_vtx!()(i,j);
                        IFace.vtx ~= get_vtx!()(i,j+1);
                    }
                    if (i == imin) {
                        IFace.is_on_boundary = true;
                        IFace.bc_id = Face.west;
                        if (bc[Face.west].ghost_cell_data_available) {
                            IFace.left_cell = (myConfig.dimensions == 3) ? get_cell!()(i-1,j,k) : get_cell!()(i-1,j);
                        }
                        IFace.right_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j,k) : get_cell!()(i,j);
                    } else if (i == imax+1) {
                        IFace.is_on_boundary = true;
                        IFace.bc_id = Face.east;
                        IFace.left_cell = (myConfig.dimensions == 3) ? get_cell!()(i-1,j,k) : get_cell!()(i-1,j);
                        if (bc[Face.east].ghost_cell_data_available) {
                            IFace.right_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j,k) : get_cell!()(i,j);
                        }
                    } else {
                        IFace.left_cell = (myConfig.dimensions == 3) ? get_cell!()(i-1,j,k) : get_cell!()(i-1,j);
                        IFace.right_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j,k) : get_cell!()(i,j);
                    }
                } // i loop
            } // j loop
        } // for k
        // ifj interfaces are south interfaces, with their unit normal pointing north.
        // In 2D, vtx0==p10, vtx1==p00.
        // t1 vector aligned with k-index direction
        // t2 vector aligned with i-index direction
        for (size_t k = kmin; k <= kmax; ++k) {
            for (size_t i = imin; i <= imax; ++i) {
                for (size_t j = jmin; j <= jmax+1; ++j) {
                    auto IFace = get_ifj!()(i,j,k);
                    IFace.vtx.length = 0;
                    if (myConfig.dimensions == 3) {
                        IFace.vtx ~= get_vtx!()(i,j,k);
                        IFace.vtx ~= get_vtx!()(i,j,k+1);
                        IFace.vtx ~= get_vtx!()(i+1,j,k+1);
                        IFace.vtx ~= get_vtx!()(i+1,j,k);
                    } else {
                        IFace.vtx ~= get_vtx!()(i+1,j);
                        IFace.vtx ~= get_vtx!()(i,j);
                    }
                    if (j == jmin) {
                        IFace.is_on_boundary = true;
                        IFace.bc_id = Face.south;
                        if (bc[Face.south].ghost_cell_data_available) {
                            IFace.left_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j-1,k) : get_cell!()(i,j-1);
                        }
                        IFace.right_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j,k) : get_cell!()(i,j);
                    } else if (j == jmax+1) {
                        IFace.is_on_boundary = true;
                        IFace.bc_id = Face.north;
                        IFace.left_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j-1,k) : get_cell!()(i,j-1);
                        if (bc[Face.north].ghost_cell_data_available) {
                            IFace.right_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j,k) : get_cell!()(i,j);
                        }
                    } else {
                        IFace.left_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j-1,k) : get_cell!()(i,j-1);
                        IFace.right_cell = (myConfig.dimensions == 3) ? get_cell!()(i,j,k) : get_cell!()(i,j);
                    }
                } // j loop
            } // i loop
        } // for k
        if (myConfig.dimensions == 2) return;
        // ifk interfaces are bottom interfaces, with unit normal pointing to top.
        // t1 vector aligned with i-index direction
        // t2 vector aligned with j-index direction
        for (size_t i = imin; i <= imax; ++i) {
            for (size_t j = jmin; j <= jmax; ++j) {
                for (size_t k = kmin; k <= kmax+1; ++k) {
                    auto IFace = get_ifk!()(i,j,k);
                    IFace.vtx.length = 0;
                    IFace.vtx ~= get_vtx!()(i,j,k);
                    IFace.vtx ~= get_vtx!()(i+1,j,k);
                    IFace.vtx ~= get_vtx!()(i+1,j+1,k);
                    IFace.vtx ~= get_vtx!()(i,j+1,k);
                    if (k == kmin) {
                        IFace.is_on_boundary = true;
                        IFace.bc_id = Face.bottom;
                        if (bc[Face.bottom].ghost_cell_data_available) {
                            IFace.left_cell = get_cell!()(i,j,k-1);
                        }
                        IFace.right_cell = get_cell!()(i,j,k);
                    } else if (k == kmax+1) {
                        IFace.is_on_boundary = true;
                        IFace.bc_id = Face.top;
                        IFace.left_cell = get_cell!()(i,j,k-1);
                        if (bc[Face.top].ghost_cell_data_available) {
                            IFace.right_cell = get_cell!()(i,j,k);
                        }
                    } else {
                        IFace.left_cell = get_cell!()(i,j,k-1);
                        IFace.right_cell = get_cell!()(i,j,k);
                    }
                } // for k
            } // j loop
        } // i loop
        return;
    } // end bind_interfaces_vertices_and_cells()

    @nogc
    override void compute_primary_cell_geometric_data(size_t gtl)
    // Compute cell and interface geometric properties.
    {
        size_t i, j, k;
        Vector3 dummy;
        Vector3 ds;
        if (myConfig.dimensions == 2) {
            foreach (c; cells) { c.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
            foreach (f; faces) { f.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
        } else { // 3D
            foreach (c; cells) { c.update_3D_geometric_data(gtl); }
            foreach (f; faces) { f.update_3D_geometric_data(gtl); }
        }
        //
        // Propagate cross-cell lengths into the ghost cells.
        // *Assuming* that the grid is not changing rapidly toward the wall,
        // we can extrapolate the cell geometry into the halo of ghost cells.
        // 25-Feb-2014
        // Jason Qin and Paul Petrie-Repar have identified the lack of exact symmetry in
        // the reconstruction process at the wall as being a cause of the leaky wall
        // boundary conditions.  Note that the symmetry is not consistent with the 
        // linear extrapolation used for the positions and volumes in the next section.
        // [TODO] -- think about this carefully.
        auto option = CopyDataOption.cell_lengths_only;
        for (j = jmin; j <= jmax; ++j) {
            for (k = kmin; k <= kmax; ++k) {
                i = imin;
                get_cell!()(i-1,j,k).copy_values_from(get_cell!()(i,j,k), option);
                get_cell!()(i-2,j,k).copy_values_from(get_cell!()(i+1,j,k), option);
                version(nghost3) { get_cell!()(i-3,j,k).copy_values_from(get_cell!()(i+2,j,k), option); }
                i = imax;
                get_cell!()(i+1,j,k).copy_values_from(get_cell!()(i,j,k), option);
                get_cell!()(i+2,j,k).copy_values_from(get_cell!()(i-1,j,k), option);
                version(nghost3) { get_cell!()(i+3,j,k).copy_values_from(get_cell!()(i-2,j,k), option); }
            }
        }
        for (i = imin; i <= imax; ++i) {
            for (k = kmin; k <= kmax; ++k) {
                j = jmin;
                get_cell!()(i,j-1,k).copy_values_from(get_cell!()(i,j,k), option);
                get_cell!()(i,j-2,k).copy_values_from(get_cell!()(i,j+1,k), option);
                version(nghost3) { get_cell!()(i,j-3,k).copy_values_from(get_cell!()(i,j+2,k), option); }
                j = jmax;
                get_cell!()(i,j+1,k).copy_values_from(get_cell!()(i,j,k), option);
                get_cell!()(i,j+2,k).copy_values_from(get_cell!()(i,j-1,k), option);
                version(nghost3) { get_cell!()(i,j+3,k).copy_values_from(get_cell!()(i,j-2,k), option); }
            }
        }
        if (myConfig.dimensions == 3) {
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    k = kmin;
                    get_cell!()(i,j,k-1).copy_values_from(get_cell!()(i,j,k), option);
                    get_cell!()(i,j,k-2).copy_values_from(get_cell!()(i,j,k+1), option);
                    version(nghost3) { get_cell!()(i,j,k-3).copy_values_from(get_cell!()(i,j,k+2), option); }
                    k = kmax;
                    get_cell!()(i,j,k+1).copy_values_from(get_cell!()(i,j,k), option);
                    get_cell!()(i,j,k+2).copy_values_from(get_cell!()(i,j,k-1), option);
                    version(nghost3) { get_cell!()(i,j,k+3).copy_values_from(get_cell!()(i,j,k-2), option); }
                }
            }
        } // end if dimensions == 3
        /* Extrapolate (with first-order) cell positions and volumes to ghost cells. */
        // TODO -- think about how to make these things consistent.
        @nogc
        void extrap(ref Vector3 pos, ref const(Vector3) p1, ref const(Vector3) p2)
        {
            pos.set(p1); pos *= 2.0; pos -= p2;
        }
        for (j = jmin; j <= jmax; ++j) {
            for (k = kmin; k <= kmax; ++k) {
                i = imin;
                auto cell_1 = get_cell!()(i,j,k);
                auto cell_2 = get_cell!()(i+1,j,k);
                auto ghost_cell = get_cell!()(i-1,j,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                cell_2 = cell_1;
                cell_1 = ghost_cell;
                ghost_cell = get_cell!()(i-2,j,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_2.volume[gtl] - cell_2.volume[gtl];
                version(nghost3) {
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = get_cell!()(i-3,j,k);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_2.volume[gtl] - cell_2.volume[gtl];
                }
                i = imax;
                cell_1 = get_cell!()(i,j,k);
                cell_2 = get_cell!()(i-1,j,k);
                ghost_cell = get_cell!()(i+1,j,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                cell_2 = cell_1;
                cell_1 = ghost_cell;
                ghost_cell = get_cell!()(i+2,j,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                version(nghost3) {
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = get_cell!()(i+3,j,k);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                }
            }
        }
        for (i = imin; i <= imax; ++i) {
            for (k = kmin; k <= kmax; ++k) {
                j = jmin;
                auto cell_1 = get_cell!()(i,j,k);
                auto cell_2 = get_cell!()(i,j+1,k);
                auto ghost_cell = get_cell!()(i,j-1,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                cell_2 = cell_1;
                cell_1 = ghost_cell;
                ghost_cell = get_cell!()(i,j-2,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                version(nghost3) {
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = get_cell!()(i,j-3,k);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                }
                j = jmax;
                cell_1 = get_cell!()(i,j,k);
                cell_2 = get_cell!()(i,j-1,k);
                ghost_cell = get_cell!()(i,j+1,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                cell_2 = cell_1;
                cell_1 = ghost_cell;
                ghost_cell = get_cell!()(i,j+2,k);
                extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                version(nghost3) {
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = get_cell!()(i,j+3,k);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                }
            }
        }
        if (myConfig.dimensions == 3) {
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    k = kmin;
                    auto cell_1 = get_cell!()(i,j,k);
                    auto cell_2 = get_cell!()(i,j,k+1);
                    auto ghost_cell = get_cell!()(i,j,k-1);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = get_cell!()(i,j,k-2);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    version(nghost3) {
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = get_cell!()(i,j,k-3);
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    }
                    k = kmax;
                    cell_1 = get_cell!()(i,j,k);
                    cell_2 = get_cell!()(i,j,k-1);
                    ghost_cell = get_cell!()(i,j,k+1);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = get_cell!()(i,j,k+2);
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    version(nghost3) {
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = get_cell!()(i,j,k+3);
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    }
                }
            }
        } // end if dimensions == 3
    } // end compute_primary_cell_geometric_data()

    @nogc
    override void compute_least_squares_setup(size_t gtl)
    {
        // Update the least-squares geometric weights and the workspaces, if appropriate.
        // The weights should be calculated when the grid is initialised or moved.
        // They are needed for flow gradient calculations that feed into the viscous fluxes.
        if (myConfig.viscous && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)) {
            if (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces) {
                foreach(iface; faces) {
                    iface.grad.set_up_workspace_leastsq(iface.cloud_pos, iface.pos, false, iface.ws_grad);
                }       
            } else { // myConfig.spatial_deriv_locn == vertices
                foreach(vtx; vertices) {
                    vtx.grad.set_up_workspace_leastsq(vtx.cloud_pos, vtx.pos[gtl], true, vtx.ws_grad);
                }
            }
        }
    } // end compute_least_squares_setup()

    void store_references_for_derivative_calc(size_t gtl)
    {
        final switch (myConfig.spatial_deriv_locn) {
        case SpatialDerivLocn.vertices:
            store_references_for_derivative_calc_at_vertices(gtl);
            break;
        case SpatialDerivLocn.faces:
            store_references_for_derivative_calc_at_faces(gtl);
            break;
        case SpatialDerivLocn.cells:
            assert(0, "spatial_deriv_locn at cells not implemented for the structured solver.");
        }
    } // end store_references_for_derivative_calc()

    void store_references_for_derivative_calc_at_faces(size_t gtl)
    {
    // The weighted least squares calculation is expecting the interface
    // at which the gradient is being calculated to be stored in position [0].
    // However the divergence calculation is expecting a specific ordering of
    // the cloud points, as such we must look up the spatial_deriv_calc type
    // to decide which cloud to use.
        size_t i, j, k;
        if (myConfig.dimensions == 2) {
            // First, i-faces
            for (i = imin; i <= imax+1; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    FVInterface face = get_ifi!()(i,j);
                    // Points nearby.
                    if (i == imin) {
                        // west boundary
                        FVInterface D = get_ifj!()(i,j);
                        FVCell E = get_cell!()(i,j);
                        FVInterface F = get_ifj!()(i,j+1);
                        // Retain locations and references to flow states for later.
                        face.cloud_pos = [&(face.pos), &(D.pos), &(E.pos[gtl]), &(F.pos)];
                        face.cloud_fs = [face.fs, D.fs, E.fs, F.fs];
                    } else if (i == imax+1) {
                        // east boundary
                        FVInterface A = get_ifj!()(i-1,j+1);
                        FVCell B = get_cell!()(i-1,j);
                        FVInterface C = get_ifj!()(i-1,j);
                        // Retain locations and references to flow states for later.
                        if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                            face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos[gtl]), &(C.pos)];
                            face.cloud_fs = [face.fs, A.fs, B.fs, C.fs];
                        } else {
                            face.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos), &(face.pos)];
                            face.cloud_fs = [A.fs, B.fs, C.fs, face.fs];
                        }
                    } else {
                        // interior face
                        FVInterface A = get_ifj!()(i-1,j+1);
                        FVCell B = get_cell!()(i-1,j);
                        FVInterface C = get_ifj!()(i-1,j);
                        FVInterface D = get_ifj!()(i,j);
                        FVCell E = get_cell!()(i,j);
                        FVInterface F = get_ifj!()(i,j+1);
                        // Retain locations and references to flow states for later.
                        if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                            face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos[gtl]), &(C.pos), 
                                              &(D.pos), &(E.pos[gtl]), &(F.pos)];
                            face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs, F.fs];
                        } else {
                            face.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos), 
                                              &(D.pos), &(E.pos[gtl]), &(F.pos)];
                            face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs, F.fs];
                        }
                    }
                } // j loop
            } // i loop
            // Now, j-faces
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax+1; ++j) {
                    FVInterface face = get_ifj!()(i,j);
                    // Points nearby.
                    if (j == jmin) {
                        // south boundary
                        FVInterface D = get_ifi!()(i+1,j);
                        FVCell E = get_cell!()(i,j);
                        FVInterface F = get_ifi!()(i,j);
                        // Retain locations and references to flow states for later.
                        face.cloud_pos = [&(face.pos), &(D.pos), &(E.pos[gtl]), &(F.pos)];
                        face.cloud_fs = [face.fs, D.fs, E.fs, F.fs];
                    } else if (j == jmax+1) {
                        // north boundary
                        FVInterface A = get_ifi!()(i,j-1);
                        FVCell B = get_cell!()(i,j-1);
                        FVInterface C = get_ifi!()(i+1,j-1);
                        // Retain locations and references to flow states for later.
                        if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                            face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos[gtl]), &(C.pos)];
                            face.cloud_fs = [face.fs, A.fs, B.fs, C.fs];
                        } else {
                            face.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos), &(face.pos)];
                            face.cloud_fs = [A.fs, B.fs, C.fs, face.fs];
                        }
                    } else {
                        // interior face
                        FVInterface A = get_ifi!()(i,j-1);
                        FVCell B = get_cell!()(i,j-1);
                        FVInterface C = get_ifi!()(i+1,j-1);
                        FVInterface D = get_ifi!()(i+1,j);
                        FVCell E = get_cell!()(i,j);
                        FVInterface F = get_ifi!()(i,j);
                        // Retain locations and references to flow states for later.
                        face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos[gtl]), &(C.pos), 
                                          &(D.pos), &(E.pos[gtl]), &(F.pos)];
                        face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs, F.fs];
                        if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                            face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos[gtl]), &(C.pos), 
                                              &(D.pos), &(E.pos[gtl]), &(F.pos)];
                            face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs, F.fs];
                        } else {
                            face.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos), 
                                              &(D.pos), &(E.pos[gtl]), &(F.pos)];
                            face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs, F.fs];
                        }
                    }
                } // j loop
            } // i loop
        } else { // for 3D.
            // First, i-faces
            for (i = imin; i <= imax+1; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    for (k = kmin; k <= kmax; ++k) {
                        FVInterface face = get_ifi!()(i,j,k);
                        // Points nearby.
                        if (i == imin) {
                            // west boundary
                            FVInterface F = get_ifj!()(i,j+1,k);
                            FVInterface G = get_ifj!()(i,j,k);
                            FVInterface H = get_ifk!()(i,j,k+1);
                            FVInterface I = get_ifk!()(i,j,k);
                            FVCell J = get_cell!()(i,j,k);
                            // Retain locations and references to flow states for later.
                            face.cloud_pos = [&(face.pos), &(F.pos), &(G.pos), &(H.pos),
                                              &(I.pos), &(J.pos[gtl])];
                            face.cloud_fs = [face.fs, F.fs, G.fs, H.fs, I.fs, J.fs];
                        } else if (i == imax+1) {
                            // east boundary
                            FVInterface A = get_ifj!()(i-1,j+1,k);
                            FVInterface B = get_ifj!()(i-1,j,k);
                            FVInterface C = get_ifk!()(i-1,j,k+1);
                            FVInterface D = get_ifk!()(i-1,j,k);
                            FVCell E = get_cell!()(i-1,j,k);
                            // Retain locations and references to flow states for later.
                            if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                                face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos), &(C.pos), &(D.pos),
                                                  &(E.pos[gtl])];
                                face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs];
                            } else {
                                face.cloud_pos = [&(A.pos), &(B.pos), &(C.pos), &(D.pos),
                                                  &(E.pos[gtl]), &(face.pos)];
                                face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs, face.fs];
                            }
                        } else {
                            // interior face
                            FVInterface A = get_ifj!()(i-1,j+1,k);
                            FVInterface B = get_ifj!()(i-1,j,k);
                            FVInterface C = get_ifk!()(i-1,j,k+1);
                            FVInterface D = get_ifk!()(i-1,j,k);
                            FVCell E = get_cell!()(i-1,j,k);
                            FVInterface F = get_ifj!()(i,j+1,k);
                            FVInterface G = get_ifj!()(i,j,k);
                            FVInterface H = get_ifk!()(i,j,k+1);
                            FVInterface I = get_ifk!()(i,j,k);
                            FVCell J = get_cell!()(i,j,k);
                            // Retain locations and references to flow states for later.
                            if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                                face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos), &(C.pos), &(D.pos), &(E.pos[gtl]), 
                                                  &(F.pos), &(G.pos), &(H.pos), &(I.pos), &(J.pos[gtl])];
                                face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs,
                                                 F.fs, G.fs, H.fs, I.fs, J.fs];
                            } else {
                                face.cloud_pos = [&(A.pos), &(B.pos), &(C.pos), &(D.pos), &(E.pos[gtl]), 
                                                  &(F.pos), &(G.pos), &(H.pos), &(I.pos), &(J.pos[gtl])];
                                face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs,
                                                 F.fs, G.fs, H.fs, I.fs, J.fs];
                            }
                        }
                    } // k loop
                } // j loop
            } // i loop
            // Next, j-faces
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax+1; ++j) {
                    for (k = kmin; k <= kmax; ++k) {
                        FVInterface face = get_ifj!()(i,j,k);
                        // Points nearby.
                        if (j == jmin) {
                            // south boundary
                            FVInterface F = get_ifi!()(i+1,j,k);
                            FVInterface G = get_ifi!()(i,j,k);
                            FVInterface H = get_ifk!()(i,j,k+1);
                            FVInterface I = get_ifk!()(i,j,k);
                            FVCell J = get_cell!()(i,j,k);
                            // Retain locations and references to flow states for later.
                            face.cloud_pos = [&(face.pos), &(F.pos), &(G.pos), &(H.pos),
                                              &(I.pos), &(J.pos[gtl])];
                            face.cloud_fs = [face.fs, F.fs, G.fs, H.fs, I.fs, J.fs];
                        } else if (j == jmax+1) {
                            // north boundary
                            FVInterface A = get_ifi!()(i+1,j-1,k);
                            FVInterface B = get_ifi!()(i,j-1,k);
                            FVInterface C = get_ifk!()(i,j-1,k+1);
                            FVInterface D = get_ifk!()(i,j-1,k);
                            FVCell E = get_cell!()(i,j-1,k);
                            // Retain locations and references to flow states for later.
                            if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                                face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos), &(C.pos), &(D.pos),
                                                  &(E.pos[gtl])];
                                face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs];
                            } else {
                                face.cloud_pos = [&(A.pos), &(B.pos), &(C.pos), &(D.pos),
                                                  &(E.pos[gtl]), &(face.pos)];
                                face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs, face.fs];
                            }
                        } else {
                            // interior face
                            FVInterface A = get_ifi!()(i+1,j-1,k);
                            FVInterface B = get_ifi!()(i,j-1,k);
                            FVInterface C = get_ifk!()(i,j-1,k+1);
                            FVInterface D = get_ifk!()(i,j-1,k);
                            FVCell E = get_cell!()(i,j-1,k);
                            FVInterface F = get_ifi!()(i+1,j,k);
                            FVInterface G = get_ifi!()(i,j,k);
                            FVInterface H = get_ifk!()(i,j,k+1);
                            FVInterface I = get_ifk!()(i,j,k);
                            FVCell J = get_cell!()(i,j,k);
                            // Retain locations and references to flow states for later.
                            if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                                face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos), &(C.pos), &(D.pos), &(E.pos[gtl]), 
                                                  &(F.pos), &(G.pos), &(H.pos), &(I.pos), &(J.pos[gtl])];
                                face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs,
                                                 F.fs, G.fs, H.fs, I.fs, J.fs];
                            } else {
                                face.cloud_pos = [&(A.pos), &(B.pos), &(C.pos), &(D.pos), &(E.pos[gtl]), 
                                                  &(F.pos), &(G.pos), &(H.pos), &(I.pos), &(J.pos[gtl])];
                                face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs,
                                                 F.fs, G.fs, H.fs, I.fs, J.fs];
                            }
                        }
                    } // k loop
                } // j loop
            } // i loop
            // Finally, k-faces
            for (i = imin; i <= imax; ++i) {
                for (j = jmin; j <= jmax; ++j) {
                    for (k = kmin; k <= kmax+1; ++k) {
                        FVInterface face = get_ifk!()(i,j,k);
                        // Points nearby.
                        if (k == kmin) {
                            // bottom boundary
                            FVInterface F = get_ifj!()(i,j+1,k);
                            FVInterface G = get_ifj!()(i,j,k);
                            FVInterface H = get_ifi!()(i+1,j,k);
                            FVInterface I = get_ifi!()(i,j,k);
                            FVCell J = get_cell!()(i,j,k);
                            // Retain locations and references to flow states for later.
                            face.cloud_pos = [&(face.pos), &(F.pos), &(G.pos), &(H.pos),
                                              &(I.pos), &(J.pos[gtl])];
                            face.cloud_fs = [face.fs, F.fs, G.fs, H.fs, I.fs, J.fs];
                        } else if (k == kmax+1) {
                            // top boundary
                            FVInterface A = get_ifj!()(i,j+1,k-1);
                            FVInterface B = get_ifj!()(i,j,k-1);
                            FVInterface C = get_ifi!()(i+1,j,k-1);
                            FVInterface D = get_ifi!()(i,j,k-1);
                            FVCell E = get_cell!()(i,j,k-1);
                            // Retain locations and references to flow states for later.
                            if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                                face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos), &(C.pos), &(D.pos),
                                                  &(E.pos[gtl])];
                                face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs];
                            } else {
                                face.cloud_pos = [&(A.pos), &(B.pos), &(C.pos), &(D.pos),
                                                  &(E.pos[gtl]), &(face.pos)];
                                face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs, face.fs];
                            }
                        } else {
                            // interior face
                            FVInterface A = get_ifj!()(i,j+1,k-1);
                            FVInterface B = get_ifj!()(i,j,k-1);
                            FVInterface C = get_ifi!()(i+1,j,k-1);
                            FVInterface D = get_ifi!()(i,j,k-1);
                            FVCell E = get_cell!()(i,j,k-1);
                            FVInterface F = get_ifj!()(i,j+1,k);
                            FVInterface G = get_ifj!()(i,j,k);
                            FVInterface H = get_ifi!()(i+1,j,k);
                            FVInterface I = get_ifi!()(i,j,k);
                            FVCell J = get_cell!()(i,j,k);
                            // Retain locations and references to flow states for later.
                            if (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares) {
                                face.cloud_pos = [&(face.pos), &(A.pos), &(B.pos), &(C.pos), &(D.pos), &(E.pos[gtl]), 
                                                  &(F.pos), &(G.pos), &(H.pos), &(I.pos), &(J.pos[gtl])];
                                face.cloud_fs = [face.fs, A.fs, B.fs, C.fs, D.fs, E.fs,
                                             F.fs, G.fs, H.fs, I.fs, J.fs];
                            } else {
                                face.cloud_pos = [&(A.pos), &(B.pos), &(C.pos), &(D.pos), &(E.pos[gtl]), 
                                                  &(F.pos), &(G.pos), &(H.pos), &(I.pos), &(J.pos[gtl])];
                                face.cloud_fs = [A.fs, B.fs, C.fs, D.fs, E.fs,
                                                 F.fs, G.fs, H.fs, I.fs, J.fs];
                            }
                        }
                    } // k loop
                } // j loop
            } // i loop
        } // end if (myConfig.dimensions
    } // end store_references_for_derivative_calc_at_faces()

    void store_references_for_derivative_calc_at_vertices(size_t gtl)
    {
        size_t i, j, k;
        if (myConfig.dimensions == 2) {
            // First, do all of the internal secondary cells.
            // i.e. Those not on a boundary.
            for ( i = imin+1; i <= imax; ++i ) {
                for ( j = jmin+1; j <= jmax; ++j ) {
                    // Secondary-cell centre is a primary-cell vertex.
                    FVVertex vtx = get_vtx!()(i,j);
                    // These are the corners of the secondary cell.
                    FVCell A = get_cell!()(i,j-1);
                    FVCell B = get_cell!()(i,j);
                    FVCell C = get_cell!()(i-1,j);
                    FVCell D = get_cell!()(i-1,j-1);
                    // Retain locations and references to flow states for later.
                    vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos[gtl]), &(C.pos[gtl]), &(D.pos[gtl])];
                    vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
                } // j loop
            } // i loop
            // Half-cells along the edges of the block.
            // East boundary
            i = imax+1;
            for (j = jmin+1; j <= jmax; ++j) {
                FVVertex vtx = get_vtx!()(i,j);
                FVInterface A = get_ifi!()(i,j-1);
                FVInterface B = get_ifi!()(i,j);
                FVCell C = get_cell!()(i-1,j);
                FVCell D = get_cell!()(i-1,j-1);
                vtx.cloud_pos = [&(A.pos), &(B.pos), &(C.pos[gtl]), &(D.pos[gtl])];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // j loop
            // West boundary
            i = imin;
            for (j = jmin+1; j <= jmax; ++j) {
                FVVertex vtx = get_vtx!()(i,j);
                // These are the corners of the secondary cell.
                FVCell A = get_cell!()(i,j-1);
                FVCell B = get_cell!()(i,j);
                FVInterface C = get_ifi!()(i,j);
                FVInterface D = get_ifi!()(i,j-1);
                vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos[gtl]), &(C.pos), &(D.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // j loop
            // North boundary
            j = jmax+1;
            for (i = imin+1; i <= imax; ++i) {
                FVVertex vtx = get_vtx!()(i,j);
                FVCell A = get_cell!()(i,j-1);
                FVInterface B = get_ifj!()(i,j);
                FVInterface C = get_ifj!()(i-1,j);
                FVCell D = get_cell!()(i-1,j-1);
                vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos), &(C.pos), &(D.pos[gtl])];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // i loop
            // South boundary
            j = jmin;
            for (i = imin+1; i <= imax; ++i) {
                FVVertex vtx = get_vtx!()(i,j);
                FVInterface A = get_ifj!()(i,j);
                FVCell B = get_cell!()(i,j);
                FVCell C = get_cell!()(i-1,j);
                FVInterface D = get_ifj!()(i-1,j);
                vtx.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos[gtl]), &(D.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // i loop
            // For the corners, we are going to use the same divergence-theorem-based
            // gradient calculator and let one edge collapse to a point, thus giving
            // it a triangle to compute over.  This should be fine. 
            // North-east corner
            {
                i = imax+1; j = jmax+1;
                FVVertex vtx = get_vtx!()(i,j);
                FVInterface A = get_ifi!()(i,j-1);
                FVInterface B = get_ifj!()(i-1,j);
                FVCell C = get_cell!()(i-1,j-1);
                vtx.cloud_pos = [&(A.pos), &(B.pos), &(C.pos[gtl])];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
            // South-east corner
            {
                i = imax+1; j = jmin;
                FVVertex vtx = get_vtx!()(i,j);
                FVInterface A = get_ifi!()(i,j);
                FVCell B = get_cell!()(i-1,j);
                FVInterface C = get_ifj!()(i-1,j);
                vtx.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
            // South-west corner
            {
                i = imin; j = jmin;
                FVVertex vtx = get_vtx!()(i,j);
                FVInterface A = get_ifj!()(i,j);
                FVCell B = get_cell!()(i,j);
                FVInterface C = get_ifi!()(i,j);
                vtx.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
            // North-west corner
            {
                i = imin; j = jmax+1;
                FVVertex vtx = get_vtx!()(i,j);
                FVCell A = get_cell!()(i,j-1);
                FVInterface B = get_ifj!()(i,j);
                FVInterface C = get_ifi!()(i,j-1);
                vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos), &(C.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
        } else { // Flow quantity derivatives for 3D.
            // Internal secondary cell geometry information
            for ( i = imin; i <= imax-1; ++i ) {
                for ( j = jmin; j <= jmax-1; ++j ) {
                    for ( k = kmin; k <= kmax-1; ++k ) {
                        FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                        FVCell c0 = get_cell!()(i,j,k);
                        FVCell c1 = get_cell!()(i+1,j,k);
                        FVCell c2 = get_cell!()(i+1,j+1,k);
                        FVCell c3 = get_cell!()(i,j+1,k);
                        FVCell c4 = get_cell!()(i,j,k+1);
                        FVCell c5 = get_cell!()(i+1,j,k+1);
                        FVCell c6 = get_cell!()(i+1,j+1,k+1);
                        FVCell c7 = get_cell!()(i,j+1,k+1);
                        vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos[gtl]), &(c3.pos[gtl]),
                                         &(c4.pos[gtl]), &(c5.pos[gtl]), &(c6.pos[gtl]), &(c7.pos[gtl])];
                        vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                    }
                }
            }
            // East boundary secondary cell geometry information
            i = imax;
            for ( j = jmin; j <= jmax-1; ++j ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                    FVCell c0 = get_cell!()(i,j,k);
                    FVInterface c1 = get_ifi!()(i+1,j,k);
                    FVInterface c2 = get_ifi!()(i+1,j+1,k);
                    FVCell c3 = get_cell!()(i,j+1,k);
                    FVCell c4 = get_cell!()(i,j,k+1);
                    FVInterface c5 = get_ifi!()(i+1,j,k+1);
                    FVInterface c6 = get_ifi!()(i+1,j+1,k+1);
                    FVCell c7 = get_cell!()(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos[gtl]),
                                     &(c4.pos[gtl]), &(c5.pos), &(c6.pos), &(c7.pos[gtl])];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // West boundary secondary cell geometry information
            i = imin - 1;
            for ( j = jmin; j <= jmax-1; ++j ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                    FVInterface c0 = get_ifi!()(i+1,j,k);
                    FVCell c1 = get_cell!()(i+1,j,k);
                    FVCell c2 = get_cell!()(i+1,j+1,k);
                    FVInterface c3 = get_ifi!()(i+1,j+1,k);
                    FVInterface c4 = get_ifi!()(i+1,j,k+1);
                    FVCell c5 = get_cell!()(i+1,j,k+1);
                    FVCell c6 = get_cell!()(i+1,j+1,k+1);
                    FVInterface c7 = get_ifi!()(i+1,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos), &(c1.pos[gtl]), &(c2.pos[gtl]), &(c3.pos),
                                     &(c4.pos), &(c5.pos[gtl]), &(c6.pos[gtl]), &(c7.pos)];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // North boundary secondary cell geometry information
            j = jmax;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                    FVCell c0 = get_cell!()(i,j,k);
                    FVCell c1 = get_cell!()(i+1,j,k);
                    FVInterface c2 = get_ifj!()(i+1,j+1,k);
                    FVInterface c3 = get_ifj!()(i,j+1,k);
                    FVCell c4 = get_cell!()(i,j,k+1);
                    FVCell c5 = get_cell!()(i+1,j,k+1);
                    FVInterface c6 = get_ifj!()(i+1,j+1,k+1);
                    FVInterface c7 = get_ifj!()(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos), &(c3.pos),
                                     &(c4.pos[gtl]), &(c5.pos[gtl]), &(c6.pos), &(c7.pos)];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // South boundary secondary cell geometry information
            j = jmin - 1;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( k = kmin; k <= kmax-1; ++k ) {
                    FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                    FVInterface c0 = get_ifj!()(i,j+1,k);
                    FVInterface c1 = get_ifj!()(i+1,j+1,k);
                    FVCell c2 = get_cell!()(i+1,j+1,k);
                    FVCell c3 = get_cell!()(i,j+1,k);
                    FVInterface c4 = get_ifj!()(i,j+1,k+1);
                    FVInterface c5 = get_ifj!()(i+1,j+1,k+1);
                    FVCell c6 = get_cell!()(i+1,j+1,k+1);
                    FVCell c7 = get_cell!()(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos), &(c1.pos), &(c2.pos[gtl]), &(c3.pos[gtl]),
                                     &(c4.pos), &(c5.pos), &(c6.pos[gtl]), &(c7.pos[gtl])];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // Top boundary secondary cell geometry information
            k = kmax;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( j = jmin; j <= jmax-1; ++j ) {
                    FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                    FVCell c0 = get_cell!()(i,j,k);
                    FVCell c1 = get_cell!()(i+1,j,k);
                    FVCell c2 = get_cell!()(i+1,j+1,k);
                    FVCell c3 = get_cell!()(i,j+1,k);
                    FVInterface c4 = get_ifk!()(i,j,k+1);
                    FVInterface c5 = get_ifk!()(i+1,j,k+1);
                    FVInterface c6 = get_ifk!()(i+1,j+1,k+1);
                    FVInterface c7 = get_ifk!()(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos[gtl]), &(c3.pos[gtl]),
                                     &(c4.pos), &(c5.pos), &(c6.pos), &(c7.pos)];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // Bottom boundary secondary cell geometry information
            k = kmin - 1;
            for ( i = imin; i <= imax-1; ++i ) {
                for ( j = jmin; j <= jmax-1; ++j ) {
                    FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                    FVInterface c0 = get_ifk!()(i,j,k+1);
                    FVInterface c1 = get_ifk!()(i+1,j,k+1);
                    FVInterface c2 = get_ifk!()(i+1,j+1,k+1);
                    FVInterface c3 = get_ifk!()(i,j+1,k+1);
                    FVCell c4 = get_cell!()(i,j,k+1);
                    FVCell c5 = get_cell!()(i+1,j,k+1);
                    FVCell c6 = get_cell!()(i+1,j+1,k+1);
                    FVCell c7 = get_cell!()(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos), &(c1.pos), &(c2.pos), &(c3.pos),
                                     &(c4.pos[gtl]), &(c5.pos[gtl]), &(c6.pos[gtl]), &(c7.pos[gtl])];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // Now, do the 4 edges around the bottom face.
            // Bottom-South edge [0]-->[1]
            j = jmin; k = kmin;         
            for ( i = imin+1; i <= imax; ++i ) {
                FVVertex vtx = get_vtx!()(i,j,k);
                FVCell c0 = get_cell!()(i-1,j,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifj!()(i-1,j,k);
                FVInterface c3 = get_ifk!()(i-1,j,k);
                FVInterface c4 = get_ifj!()(i,j,k);
                FVInterface c5 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Bottom-North edge [3]-->[2]
            j = jmax; k = kmin;
            for ( i = imin+1; i <= imax; ++i ) {
                FVVertex vtx = get_vtx!()(i,j+1,k);
                FVCell c0 = get_cell!()(i-1,j,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifj!()(i-1,j+1,k);
                FVInterface c3 = get_ifk!()(i-1,j,k);
                FVInterface c4 = get_ifj!()(i,j+1,k);
                FVInterface c5 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Bottom-West edge [0]-->[3]
            i = imin; k = kmin;
            for ( j = jmin+1; j <= jmax; ++j ) {
                FVVertex vtx = get_vtx!()(i,j,k);
                FVCell c0 = get_cell!()(i,j-1,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i,j-1,k);
                FVInterface c3 = get_ifk!()(i,j-1,k);
                FVInterface c4 = get_ifi!()(i,j,k);
                FVInterface c5 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Bottom-East edge [1]-->[2]
            i = imax; k = kmin;
            for ( j = jmin+1; j <= jmax; ++j ) {
                FVVertex vtx = get_vtx!()(i+1,j,k);
                FVCell c0 = get_cell!()(i,j-1,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i+1,j-1,k);
                FVInterface c3 = get_ifk!()(i,j-1,k);
                FVInterface c4 = get_ifi!()(i+1,j,k);
                FVInterface c5 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // 4 edges around the top face.
            // Top-South edge [4]-->[5]
            j = jmin; k = kmax;
            for ( i = imin+1; i <= imax; ++i ) {
                FVVertex vtx = get_vtx!()(i,j,k+1);
                FVCell c0 = get_cell!()(i-1,j,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifj!()(i-1,j,k);
                FVInterface c3 = get_ifk!()(i-1,j,k+1);
                FVInterface c4 = get_ifj!()(i,j,k);
                FVInterface c5 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Top-North edge [7]-->[6]
            j = jmax; k = kmax;
            for ( i = imin+1; i <= imax; ++i ) {
                FVVertex vtx = get_vtx!()(i,j+1,k+1);
                FVCell c0 = get_cell!()(i-1,j,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifj!()(i-1,j+1,k);
                FVInterface c3 = get_ifk!()(i-1,j,k+1);
                FVInterface c4 = get_ifj!()(i,j+1,k);
                FVInterface c5 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Top-West edge [4]-->[7]
            i = imin; k = kmax;
            for ( j = jmin+1; j <= jmax; ++j ) {
                FVVertex vtx = get_vtx!()(i,j,k+1);
                FVCell c0 = get_cell!()(i,j-1,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i,j-1,k);
                FVInterface c3 = get_ifk!()(i,j-1,k+1);
                FVInterface c4 = get_ifi!()(i,j,k);
                FVInterface c5 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Top-East edge [5]-->[6]
            i = imax; k = kmax;
            for ( j = jmin+1; j <= jmax; ++j ) {
                FVVertex vtx = get_vtx!()(i+1,j,k+1);
                FVCell c0 = get_cell!()(i,j-1,k);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i+1,j-1,k);
                FVInterface c3 = get_ifk!()(i,j-1,k+1);
                FVInterface c4 = get_ifi!()(i+1,j,k);
                FVInterface c5 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // 4 edges running from bottom to top.
            // South-West edge [0]-->[4]
            i = imin; j = jmin;
            for ( k = kmin+1; k <= kmax; ++k ) {
                FVVertex vtx = get_vtx!()(i,j,k);
                FVCell c0 = get_cell!()(i,j,k-1);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i,j,k-1);
                FVInterface c3 = get_ifj!()(i,j,k-1);
                FVInterface c4 = get_ifi!()(i,j,k);
                FVInterface c5 = get_ifj!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // South-East edge [1]-->[5]
            i = imax; j = jmin;
            for ( k = kmin+1; k <= kmax; ++k ) {
                FVVertex vtx = get_vtx!()(i+1,j,k);
                FVCell c0 = get_cell!()(i,j,k-1);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i+1,j,k-1);
                FVInterface c3 = get_ifj!()(i,j,k-1);
                FVInterface c4 = get_ifi!()(i+1,j,k);
                FVInterface c5 = get_ifj!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // North-East edge [2]-->[6]
            i = imax; j = jmax;
            for ( k = kmin+1; k <= kmax; ++k ) {
                FVVertex vtx = get_vtx!()(i+1,j+1,k);
                FVCell c0 = get_cell!()(i,j,k-1);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i+1,j,k-1);
                FVInterface c3 = get_ifj!()(i,j+1,k-1);
                FVInterface c4 = get_ifi!()(i+1,j,k);
                FVInterface c5 = get_ifj!()(i,j+1,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // North-West edge [3]-->[7]
            i = imin; j = jmax;
            for ( k = kmin+1; k <= kmax; ++k ) {
                FVVertex vtx = get_vtx!()(i,j+1,k);
                FVCell c0 = get_cell!()(i,j,k-1);
                FVCell c1 = get_cell!()(i,j,k);
                FVInterface c2 = get_ifi!()(i,j,k-1);
                FVInterface c3 = get_ifj!()(i,j+1,k-1);
                FVInterface c4 = get_ifi!()(i,j,k);
                FVInterface c5 = get_ifj!()(i,j+1,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Finally, the 8 corners.
            // South-West-Bottom corner [0]
            i = imin; j = jmin; k = kmin;
            {
                FVVertex vtx = get_vtx!()(i,j,k);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i,j,k);
                FVInterface c2 = get_ifj!()(i,j,k);
                FVInterface c3 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // South-East-Bottom corner [1]
            i = imax; j = jmin; k = kmin;
            {
                FVVertex vtx = get_vtx!()(i+1,j,k);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i+1,j,k);
                FVInterface c2 = get_ifj!()(i,j,k);
                FVInterface c3 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-East-Bottom corner [2]
            i = imax; j = jmax; k = kmin;
            {
                FVVertex vtx = get_vtx!()(i+1,j+1,k);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i+1,j,k);
                FVInterface c2 = get_ifj!()(i,j+1,k);
                FVInterface c3 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-West-Bottom corner [3]
            i = imin; j = jmax; k = kmin;
            {
                FVVertex vtx = get_vtx!()(i,j+1,k);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i,j,k);
                FVInterface c2 = get_ifj!()(i,j+1,k);
                FVInterface c3 = get_ifk!()(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // South-West-Top corner [4]
            i = imin; j = jmin; k = kmax;
            {
                FVVertex vtx = get_vtx!()(i,j,k+1);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i,j,k);
                FVInterface c2 = get_ifj!()(i,j,k);
                FVInterface c3 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // South-East-Top corner [5]
            i = imax; j = jmin; k = kmax;
            {
                FVVertex vtx = get_vtx!()(i+1,j,k+1);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i+1,j,k);
                FVInterface c2 = get_ifj!()(i,j,k);
                FVInterface c3 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-East-Top corner [6]
            i = imax; j = jmax; k = kmax;
            {
                FVVertex vtx = get_vtx!()(i+1,j+1,k+1);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i+1,j,k);
                FVInterface c2 = get_ifj!()(i,j+1,k);
                FVInterface c3 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-West-Top corner [7]
            i = imin; j = jmax; k = kmax;
            {
                FVVertex vtx = get_vtx!()(i,j+1,k+1);
                FVCell c0 = get_cell!()(i,j,k);
                FVInterface c1 = get_ifi!()(i,j,k);
                FVInterface c2 = get_ifj!()(i,j+1,k);
                FVInterface c3 = get_ifk!()(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
        } // end if (myConfig.dimensions
    } // end store_references_for_derivative_calc_at_vertices()

    @nogc
    override void sync_vertices_from_underlying_grid(size_t gtl=0)
    {
        size_t nivtx, njvtx, nkvtx;
        nivtx = grid.niv; njvtx = grid.njv; nkvtx = grid.nkv;
        if (myConfig.dimensions == 3) {
            if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx-1 != nkcell ) {
                string msg = "Mismatch in 3D grid size";
                debug {
                    msg ~= text("\nFor block[", id, "] we have read nivtx=", nivtx,
                                " njvtx=", njvtx, " nkvtx=", nkvtx);
                }
                throw new FlowSolverException(msg);
            }
            for (size_t k = kmin; k <= kmax+1; ++k) {
                for (size_t j = jmin; j <= jmax+1; ++j) {
                    for (size_t i = imin; i <= imax+1; ++i) {
                        auto vtx = get_vtx!()(i,j,k);
                        auto src_vtx = grid[i-imin,j-jmin,k-kmin];
                        vtx.pos[gtl].set(src_vtx);
                    } // for i
                } // for j
            } // for k
        } else { // 2D case
            if (nivtx-1 != nicell || njvtx-1 != njcell || nkvtx != 1) {
                string msg = "Mismatch in 2D grid size";
                debug {
                    msg ~= text("\nFor block[", id, "] we have read nivtx=", nivtx,
                                " njvtx=", njvtx, " nkvtx=", nkvtx);
                }
                throw new FlowSolverException(msg);
            }
            for (size_t j = jmin; j <= jmax+1; ++j) {
                for (size_t i = imin; i <= imax+1; ++i) {
                    auto vtx = get_vtx!()(i,j);
                    auto src_vtx = grid[i-imin,j-jmin];
                    vtx.pos[gtl].set(src_vtx.x, src_vtx.y, to!number(0.0));
                } // for i
            } // for j
        }
    } // end sync_vertices_from_underlying_grid()

    @nogc
    override void sync_vertices_to_underlying_grid(size_t gtl=0)
    // Note that we reuse the StructuredGrid object that was created on the
    // use of init_grid_and_flow_arrays().
    {
        size_t kmaxrange = (myConfig.dimensions == 3) ? kmax+1 : kmax;
        for (size_t k = kmin; k <= kmaxrange; ++k) {
            for (size_t j = jmin; j <= jmax+1; ++j) {
                for (size_t i = imin; i <= imax+1; ++i) {
                    auto vtx = get_vtx!()(i,j,k);
                    auto dest_vtx = grid[i-imin,j-jmin,k-kmin];
                    dest_vtx.set(vtx.pos[gtl]);
                } // for i
            } // for j
        } // for k
    } // end sync_vertices_to_underlying_grid()

    override void read_new_underlying_grid(string fileName)
    {
        if (myConfig.verbosity_level > 1) { writeln("read_new_underlying_grid() for block ", id); }
        grid = new StructuredGrid(fileName, myConfig.grid_format);
        grid.sort_cells_into_bins();
    }

    override void write_underlying_grid(string fileName)
    {
        if (myConfig.verbosity_level > 1) { writeln("write_underlying_grid() for block ", id); }
        grid.write(fileName, myConfig.grid_format);
    }
    
    override double read_solution(string filename, bool overwrite_geometry_data)
    // Note that the position data is read into grid-time-level 0
    // by scan_values_from_string(). 
    // Returns sim_time from file.
    // Keep in sync with write_initial_flow_file() in flowstate.d
    // and write_solution below.
    {
        if (myConfig.verbosity_level > 1) { writeln("read_solution(): Start block ", id); }
        double sim_time; // to be read from file
        string myLabel;
        int nvariables;
        string[] variableNames;
        switch (myConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File fin = File(filename, "rb");
            string expected_header = "structured_grid_flow 1.0";
            char[] found_header = new char[expected_header.length];
            fin.rawRead(found_header);
            if (found_header != expected_header) {
                throw new FlowSolverException("SFluidBlock.read_solution from raw_binary_file: " ~
                                              "unexpected header: " ~ to!string(found_header)); 
            }
            int[1] int1; fin.rawRead(int1);
            int label_length = int1[0];
            if (label_length > 0) {
                char[] found_label = new char[label_length];
                fin.rawRead(found_label);
                myLabel = to!string(found_label);
            }
            double[1] dbl1; fin.rawRead(dbl1); sim_time = dbl1[0];
            fin.rawRead(int1); nvariables = int1[0];
            foreach(i; 0 .. nvariables) {
                char[] varname; fin.rawRead(int1); varname.length = int1[0]; 
                fin.rawRead(varname);
            }
            int[4] int4; fin.rawRead(int4);
            int my_dimensions = int4[0];
            size_t nic = int4[1]; size_t njc = int4[2]; size_t nkc = int4[3];
            if (my_dimensions != myConfig.dimensions) {
                string msg = text("dimensions found: " ~ to!string(my_dimensions));
                throw new FlowSolverException(msg);
            }
            if (nic != nicell || njc != njcell || 
                nkc != ((myConfig.dimensions == 3) ? nkcell : 1)) {
                string msg = text("For block[", id, "] we have a mismatch in solution size.",
                                  " Have read nic=", nic, " njc=", njc, " nkc=", nkc);
                throw new FlowSolverException(msg);
            }   
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t j = jmin; j <= jmax; ++j ) {
                    for ( size_t i = imin; i <= imax; ++i ) {
                        get_cell!()(i,j,k).read_values_from_raw_binary(fin, overwrite_geometry_data);
                    }
                }
            }
            break;
        default:
            auto byLine = new GzipByLine(filename);
            auto line = byLine.front; byLine.popFront();
            string format_version;
            formattedRead(line, "structured_grid_flow %s", &format_version);
            if (format_version != "1.0") {
                string msg = text("File format version found: " ~ format_version);
                throw new FlowSolverException(msg); 
            }
            line = byLine.front; byLine.popFront();
            formattedRead(line, "label: %s", &myLabel);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "sim_time: %g", &sim_time);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "variables: %d", &nvariables);
            line = byLine.front; byLine.popFront();
            variableNames = line.strip().split();
            foreach (i; 0 .. variableNames.length) { variableNames[i] = variableNames[i].strip("\""); }
            // If all of the variables line up,
            // it is probably faster to use fixed-order, so default to that.
            bool useFixedOrder = true;
            foreach (i, varName; myConfig.flow_variable_list) {
                if (!std.algorithm.canFind(variableNames, varName)) {
                    throw new Exception("Could not find variable: " ~ varName);
                }
                if (!std.algorithm.equal(variableNames[i], varName)) {
                    useFixedOrder = false;
                }
            }
            line = byLine.front; byLine.popFront();
            size_t my_dimensions, nic, njc, nkc;
            formattedRead(line, "dimensions: %d", &my_dimensions);
            if (my_dimensions != myConfig.dimensions) {
                string msg = text("dimensions found: " ~ to!string(my_dimensions));
                throw new FlowSolverException(msg);
            }
            line = byLine.front; byLine.popFront();
            formattedRead(line, "nicell: %d", &nic);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "njcell: %d", &njc);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "nkcell: %d", &nkc);
            if (nic != nicell || njc != njcell || 
                nkc != ((myConfig.dimensions == 3) ? nkcell : 1)) {
                string msg = text("For block[", id, "] we have a mismatch in solution size.",
                                  " Have read nic=", nic, " njc=", njc, " nkc=", nkc);
                throw new FlowSolverException(msg);
            }   
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t j = jmin; j <= jmax; ++j ) {
                    for ( size_t i = imin; i <= imax; ++i ) {
                        line = byLine.front; byLine.popFront();
                        get_cell!()(i,j,k).scan_values_from_string(line, variableNames, useFixedOrder,
                                                                myConfig.gmodel, overwrite_geometry_data);
                    } // for i
                } // for j
            } // for k
        } // end switch flow_format
        if (myConfig.verbosity_level > 1) { writeln("read_solution(): Done block ", id); }
        return sim_time;
    } // end read_solution()

    override void write_solution(string filename, double sim_time)
    // Write the flow solution (i.e. the primary variables at the cell centers)
    // for a single block.
    // The text format is almost Tecplot POINT format; the raw binary format has the same layout.
    // Keep in sync with write_initial_flow_file() in flowstate.d and read_solution above.
    {
        if (myConfig.verbosity_level > 1) { writeln("write_solution(): Start block ", id); }
        switch (myConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File outfile = File(filename, "wb");
            int[1] int1; int[4] int4; double[1] dbl1; // buffer arrays
            string header = "structured_grid_flow 1.0";
            outfile.rawWrite(to!(char[])(header));
            int1[0] = to!int(label.length); outfile.rawWrite(int1);
            if (label.length > 0) { outfile.rawWrite(to!(char[])(label)); }
            dbl1[0] = sim_time; outfile.rawWrite(dbl1);
            int1[0] = to!int(myConfig.flow_variable_list.length); outfile.rawWrite(int1);
            foreach(varname; myConfig.flow_variable_list) {
                int1[0] = to!int(varname.length); outfile.rawWrite(int1);
                outfile.rawWrite(to!(char[])(varname));
            }
            int4[0] = myConfig.dimensions;
            int4[1] = to!int(nicell); int4[2] = to!int(njcell); int4[3] = to!int(nkcell);
            outfile.rawWrite(int4);
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t j = jmin; j <= jmax; ++j ) {
                    for ( size_t i = imin; i <= imax; ++i ) {
                        get_cell!()(i,j,k).write_values_to_raw_binary(outfile);
                    }
                }
            }
            outfile.close();
            break;
        default:
            auto outfile = new GzipOut(filename);
            auto writer = appender!string();
            formattedWrite(writer, "structured_grid_flow 1.0\n");
            formattedWrite(writer, "label: %s\n", label);
            formattedWrite(writer, "sim_time: %.18e\n", sim_time);
            formattedWrite(writer, "variables: %d\n", myConfig.flow_variable_list.length);
            foreach(varname; myConfig.flow_variable_list) { formattedWrite(writer, " \"%s\"", varname); }
            formattedWrite(writer, "\n");
            formattedWrite(writer, "dimensions: %d\n", myConfig.dimensions);
            formattedWrite(writer, "nicell: %d\n", nicell);
            formattedWrite(writer, "njcell: %d\n", njcell);
            formattedWrite(writer, "nkcell: %d\n", nkcell);
            outfile.compress(writer.data);
            for ( size_t k = kmin; k <= kmax; ++k ) {
                for ( size_t j = jmin; j <= jmax; ++j ) {
                    for ( size_t i = imin; i <= imax; ++i ) {
                        outfile.compress(" " ~ get_cell!()(i,j,k).write_values_to_string() ~ "\n");
                    }
                }
            }
            outfile.finish();
        } // end switch flow_format
        if (myConfig.verbosity_level > 1) { writeln("write_solution(): Done block ", id); }
    } // end write_solution()

    @nogc
    override void propagate_inflow_data_west_to_east()
    {
        // Assume that the west-face ghost cells have appropriate data.
        for ( size_t k = kmin; k <= kmax; ++k ) {
            for ( size_t j = jmin; j <= jmax; ++j ) {
                auto src_cell = get_cell!()(imin-1,j,k);
                for ( size_t i = imin; i <= imax; ++i ) {
                    auto dest_cell = get_cell!()(i,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.all_flow);
                }
            }
        }
        foreach (cell; cells) {
            cell.encode_conserved(0, 0, omegaz);
            // Even though the following call appears redundant at this point,
            // fills in some gas properties such as Prandtl number that is
            // needed for both the cfd_check and the BaldwinLomax turbulence model.
            if (0 != cell.decode_conserved(0, 0, omegaz)) {
                throw new FlowSolverException("Bad cell while propagating west to east.");
            }
        }
        set_cell_dt_chem(-1.0);
    } // end propagate_inflow_data_west_to_east()

    @nogc
    override void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0)
    // Compute the flux from flow-field data on either-side of the interface.
    {
        // Barring exceptions at the block boundaries, the general process is:
        // (1) interpolate LEFT and RIGHT interface states from cell-center properties.
        // (2) save u, v, w, T for the viscous flux calculation by making a local average.
        // The values for u, v and T may be updated subsequently by the interface-flux function.
        // (3) Apply the flux calculator to the Lft,Rght flow states.
        //
        bool do_reconstruction = allow_high_order_interpolation && (myConfig.interpolation_order > 1);
        //
        if (myConfig.high_order_flux_calculator && n_ghost_cell_layers == 3) {
            if (!bc[Face.north].ghost_cell_data_available) { throw new Error("north ghost cell data missing"); }
            if (!bc[Face.south].ghost_cell_data_available) { throw new Error("south ghost cell data missing"); }
            if (!bc[Face.west].ghost_cell_data_available) { throw new Error("west ghost cell data missing"); }
            if (!bc[Face.east].ghost_cell_data_available) { throw new Error("east ghost cell data missing"); }
            FVCell cL0, cL1, cL2, cR0, cR1, cR2;
            // ifi interfaces are East-facing interfaces.
            for (size_t k = kmin; k <= kmax; ++k) {
                for (size_t j = jmin; j <= jmax; ++j) {
                    for (size_t i = imin; i <= imax+1; ++i) {
                        auto IFace = get_ifi!()(i,j,k);
                        cL0 = get_cell!()(i-1,j,k); cL1 = get_cell!()(i-2,j,k); cL2 = get_cell!()(i-3,j,k);
                        cR0 = get_cell!()(i,j,k); cR1 = get_cell!()(i+1,j,k); cR2 = get_cell!()(i+2,j,k);
                        // Low-order reconstruction just copies data from adjacent FV_Cell.
                        // Even for high-order reconstruction, we depend upon this copy for
                        // the viscous-transport and diffusion coefficients.
                        Lft.copy_values_from(cL0.fs);
                        Rght.copy_values_from(cR0.fs);
                        if (do_reconstruction && !IFace.in_suppress_reconstruction_zone &&
                            !(myConfig.suppress_reconstruction_at_shocks && IFace.fs.S)) {
                            one_d.interp_l3r3(IFace,
                                              cL2, cL1, cL0, cR0, cR1, cR2,
                                              cL2.iLength, cL1.iLength, cL0.iLength,
                                              cR0.iLength, cR1.iLength, cR2.iLength,
                                              Lft, Rght);
                        }
                        IFace.fs.copy_average_values_from(Lft, Rght);
                        //
                        if ((i == imin) && bc[Face.west].convective_flux_computed_in_bc) continue;
                        if ((i == imax+1) && bc[Face.east].convective_flux_computed_in_bc) continue;
                        compute_interface_flux(Lft, Rght, IFace, myConfig, omegaz);
                    } // i loop
                } // j loop
            } // for k
            // ifj interfaces are North-facing interfaces.
            for (size_t k = kmin; k <= kmax; ++k) {
                for (size_t j = jmin; j <= jmax+1; ++j) {
                    for (size_t i = imin; i <= imax; ++i) {
                        auto IFace = get_ifj!()(i,j,k);
                        cL0 = get_cell!()(i,j-1,k); cL1 = get_cell!()(i,j-2,k); cL2 = get_cell!()(i,j-3,k);
                        cR0 = get_cell!()(i,j,k); cR1 = get_cell!()(i,j+1,k); cR2 = get_cell!()(i,j+2,k);
                        // Low-order reconstruction just copies data from adjacent FV_Cell.
                        // Even for high-order reconstruction, we depend upon this copy for
                        // the viscous-transport and diffusion coefficients.
                        Lft.copy_values_from(cL0.fs);
                        Rght.copy_values_from(cR0.fs);
                        if (do_reconstruction && !IFace.in_suppress_reconstruction_zone &&
                            !(myConfig.suppress_reconstruction_at_shocks && IFace.fs.S)) {
                            one_d.interp_l3r3(IFace,
                                              cL2, cL1, cL0, cR0, cR1, cR2,
                                              cL2.iLength, cL1.iLength, cL0.iLength,
                                              cR0.iLength, cR1.iLength, cR2.iLength,
                                              Lft, Rght);
                        }
                        IFace.fs.copy_average_values_from(Lft, Rght);
                        //
                        if ((j == jmin) && bc[Face.south].convective_flux_computed_in_bc) continue;
                        if ((j == jmax+1) && bc[Face.north].convective_flux_computed_in_bc) continue;
                        compute_interface_flux(Lft, Rght, IFace, myConfig, omegaz);
                    } // j loop
                } // i loop
            } // for k
    
            if (myConfig.dimensions == 2) return;
    
            // ifk interfaces are Top-facing interfaces.
            if (!bc[Face.top].ghost_cell_data_available) { throw new Error("top ghost cell data missing"); }
            if (!bc[Face.bottom].ghost_cell_data_available) { throw new Error("bottom ghost cell data missing"); }
            for (size_t k = kmin; k <= kmax+1; ++k) {
                for (size_t j = jmin; j <= jmax; ++j) {
                    for (size_t i = imin; i <= imax; ++i) {
                        auto IFace = get_ifk!()(i,j,k);
                        cL0 = get_cell!()(i,j,k-1); cL1 = get_cell!()(i,j,k-2); cL2 = get_cell!()(i,j,k-3);
                        cR0 = get_cell!()(i,j,k); cR1 = get_cell!()(i,j,k+1); cR2 = get_cell!()(i,j,k+2);
                        // Low-order reconstruction just copies data from adjacent FV_Cell.
                        // Even for high-order reconstruction, we depend upon this copy for
                        // the viscous-transport and diffusion coefficients.
                        Lft.copy_values_from(cL0.fs);
                        Rght.copy_values_from(cR0.fs);
                        if (do_reconstruction && !IFace.in_suppress_reconstruction_zone &&
                            !(myConfig.suppress_reconstruction_at_shocks && IFace.fs.S)) {
                            one_d.interp_l3r3(IFace,
                                              cL2, cL1, cL0, cR0, cR1, cR2,
                                              cL2.iLength, cL1.iLength, cL0.iLength,
                                              cR0.iLength, cR1.iLength, cR2.iLength,
                                              Lft, Rght);
                        }
                        IFace.fs.copy_average_values_from(Lft, Rght);
                        //
                        if ((k == kmin) && bc[Face.bottom].convective_flux_computed_in_bc) continue;
                        if ((k == kmax+1) && bc[Face.top].convective_flux_computed_in_bc) continue;
                        compute_interface_flux(Lft, Rght, IFace, myConfig, omegaz);
                    } // for k 
                } // j loop
            } // i loop
            return;
        } // end if (high_order_flux_calculator)
        //
        // Continue with the flux calculation being done in the classic
        // (piecewise-parabolic) reconstruction,
        // followed by flux calculation from Left,Right conditions.
        //
        // ifi interfaces are East-facing interfaces.
        for (size_t k = kmin; k <= kmax; ++k) {
            for (size_t j = jmin; j <= jmax; ++j) {
                for (size_t i = imin; i <= imax+1; ++i) {
                    auto IFace = get_ifi!()(i,j,k);
                    auto cL0 = get_cell!()(i-1,j,k); auto cL1 = get_cell!()(i-2,j,k);
                    auto cR0 = get_cell!()(i,j,k); auto cR1 = get_cell!()(i+1,j,k);
                    // Low-order reconstruction just copies data from adjacent FV_Cell.
                    // Even for high-order reconstruction, we depend upon this copy for
                    // the viscous-transport and diffusion coefficients.
                    if ((i == imin) && !(bc[Face.west].ghost_cell_data_available)) {
                        Lft.copy_values_from(cR0.fs);
                    } else {
                        Lft.copy_values_from(cL0.fs);
                    }
                    if ((i == imax+1) && !(bc[Face.east].ghost_cell_data_available)) {
                        Rght.copy_values_from(cL0.fs);
                    } else {
                        Rght.copy_values_from(cR0.fs);
                    }
                    if (do_reconstruction && !IFace.in_suppress_reconstruction_zone &&
                        !(myConfig.suppress_reconstruction_at_shocks && IFace.fs.S)) {
                        if ((i == imin) && !(bc[Face.west].ghost_cell_data_available)) {
                            one_d.interp_l0r2(IFace, cR0, cR1, cR0.iLength, cR1.iLength, Lft, Rght);
                        } else if ((i == imin+1) && !(bc[Face.west].ghost_cell_data_available)) {
                            one_d.interp_l1r2(IFace, cL0, cR0, cR1, cL0.iLength, cR0.iLength, cR1.iLength, Lft, Rght);
                        } else if ((i == imax) && !(bc[Face.east].ghost_cell_data_available)) {
                            one_d.interp_l2r1(IFace, cL1, cL0, cR0, cL1.iLength, cL0.iLength, cR0.iLength, Lft, Rght);
                        } else if ((i == imax+1) && !(bc[Face.east].ghost_cell_data_available)) {
                            one_d.interp_l2r0(IFace, cL1, cL0, cL1.iLength, cL0.iLength, Lft, Rght);
                        } else { // General symmetric reconstruction.
                            one_d.interp_l2r2(IFace, cL1, cL0, cR0, cR1,
                                              cL1.iLength, cL0.iLength, cR0.iLength, cR1.iLength,
                                              Lft, Rght);
                        }
                    }
                    IFace.fs.copy_average_values_from(Lft, Rght);
                    //
                    if ((i == imin) && bc[Face.west].convective_flux_computed_in_bc) continue;
                    if ((i == imax+1) && bc[Face.east].convective_flux_computed_in_bc) continue;
                    //
                    if ((i == imin) && !(bc[Face.west].ghost_cell_data_available)) {
                        compute_flux_at_left_wall(Rght, IFace, myConfig, omegaz);
                    } else if ((i == imax+1) && !(bc[Face.east].ghost_cell_data_available)) {
                        compute_flux_at_right_wall(Lft, IFace, myConfig, omegaz);
                    } else {
                        compute_interface_flux(Lft, Rght, IFace, myConfig, omegaz);
                    }
                } // i loop
            } // j loop
        } // for k
        // ifj interfaces are North-facing interfaces.
        for (size_t k = kmin; k <= kmax; ++k) {
            for (size_t j = jmin; j <= jmax+1; ++j) {
                for (size_t i = imin; i <= imax; ++i) {
                    auto IFace = get_ifj!()(i,j,k);
                    auto cL0 = get_cell!()(i,j-1,k); auto cL1 = get_cell!()(i,j-2,k);
                    auto cR0 = get_cell!()(i,j,k); auto cR1 = get_cell!()(i,j+1,k);
                    // Low-order reconstruction just copies data from adjacent FV_Cell.
                    // Even for high-order reconstruction, we depend upon this copy for
                    // the viscous-transport and diffusion coefficients.
                    if ((j == jmin) && !(bc[Face.south].ghost_cell_data_available)) {
                        Lft.copy_values_from(cR0.fs);
                    } else {
                        Lft.copy_values_from(cL0.fs);
                    }
                    if ((j == jmax+1) && !(bc[Face.north].ghost_cell_data_available)) {
                        Rght.copy_values_from(cL0.fs);
                    } else {
                        Rght.copy_values_from(cR0.fs);
                    }
                    if (do_reconstruction && !IFace.in_suppress_reconstruction_zone &&
                        !(myConfig.suppress_reconstruction_at_shocks && IFace.fs.S)) {
                        if ((j == jmin) && !(bc[Face.south].ghost_cell_data_available)) {
                            one_d.interp_l0r2(IFace, cR0, cR1, cR0.jLength, cR1.jLength, Lft, Rght);
                        } else if ((j == jmin+1) && !(bc[Face.south].ghost_cell_data_available)) {
                            one_d.interp_l1r2(IFace, cL0, cR0, cR1, cL0.jLength, cR0.jLength, cR1.jLength, Lft, Rght);
                        } else if ((j == jmax) && !(bc[Face.north].ghost_cell_data_available)) {
                            one_d.interp_l2r1(IFace, cL1, cL0, cR0, cL1.jLength, cL0.jLength, cR0.jLength, Lft, Rght);
                        } else if ((j == jmax+1) && !(bc[Face.north].ghost_cell_data_available)) {
                            one_d.interp_l2r0(IFace, cL1, cL0, cL1.jLength, cL0.jLength, Lft, Rght);
                        } else { // General symmetric reconstruction.
                            one_d.interp_l2r2(IFace, cL1, cL0, cR0, cR1,
                                              cL1.jLength, cL0.jLength, cR0.jLength, cR1.jLength,
                                              Lft, Rght);
                        }
                    }
                    IFace.fs.copy_average_values_from(Lft, Rght);
                    //
                    if ((j == jmin) && bc[Face.south].convective_flux_computed_in_bc) continue;
                    if ((j == jmax+1) && bc[Face.north].convective_flux_computed_in_bc) continue;
                    //
                    if ((j == jmin) && !(bc[Face.south].ghost_cell_data_available)) {
                        compute_flux_at_left_wall(Rght, IFace, myConfig, omegaz);
                    } else if ((j == jmax+1) && !(bc[Face.north].ghost_cell_data_available)) {
                        compute_flux_at_right_wall(Lft, IFace, myConfig, omegaz);
                    } else {
                        compute_interface_flux(Lft, Rght, IFace, myConfig, omegaz);
                    }
                } // j loop
            } // i loop
        } // for k
    
        if (myConfig.dimensions == 2) return;
    
        // ifk interfaces are Top-facing interfaces.
        for (size_t k = kmin; k <= kmax+1; ++k) {
            for (size_t j = jmin; j <= jmax; ++j) {
                for (size_t i = imin; i <= imax; ++i) {
                    auto IFace = get_ifk!()(i,j,k);
                    auto cL0 = get_cell!()(i,j,k-1); auto cL1 = get_cell!()(i,j,k-2);
                    auto cR0 = get_cell!()(i,j,k); auto cR1 = get_cell!()(i,j,k+1);
                    // Low-order reconstruction just copies data from adjacent FV_Cell.
                    // Even for high-order reconstruction, we depend upon this copy for
                    // the viscous-transport and diffusion coefficients.
                    if ((k == kmin) && !(bc[Face.bottom].ghost_cell_data_available)) {
                        Lft.copy_values_from(cR0.fs);
                    } else {
                        Lft.copy_values_from(cL0.fs);
                    }
                    if ((k == kmax+1) && !(bc[Face.top].ghost_cell_data_available)) {
                        Rght.copy_values_from(cL0.fs);
                    } else {
                        Rght.copy_values_from(cR0.fs);
                    }
                    if (do_reconstruction && !IFace.in_suppress_reconstruction_zone &&
                        !(myConfig.suppress_reconstruction_at_shocks && IFace.fs.S)) {
                        if ((k == kmin) && !(bc[Face.bottom].ghost_cell_data_available)) {
                            one_d.interp_l0r2(IFace, cR0, cR1, cR0.kLength, cR1.kLength, Lft, Rght);
                        } else if ((k == kmin+1) && !(bc[Face.bottom].ghost_cell_data_available)) {
                            one_d.interp_l1r2(IFace, cL0, cR0, cR1, cL0.kLength, cR0.kLength, cR1.kLength, Lft, Rght);
                        } else if ((k == kmax) && !(bc[Face.top].ghost_cell_data_available)) {
                            one_d.interp_l2r1(IFace, cL1, cL0, cR0, cL1.kLength, cL0.kLength, cR0.kLength, Lft, Rght);
                        } else if ((k == kmax+1) && !(bc[Face.top].ghost_cell_data_available)) {
                            one_d.interp_l2r0(IFace, cL1, cL0, cL1.kLength, cL0.kLength, Lft, Rght);
                        } else { // General symmetric reconstruction.
                            one_d.interp_l2r2(IFace, cL1, cL0, cR0, cR1,
                                              cL1.kLength, cL0.kLength, cR0.kLength, cR1.kLength,
                                              Lft, Rght);
                        }
                    }
                    IFace.fs.copy_average_values_from(Lft, Rght);
                    //
                    if ((k == kmin) && bc[Face.bottom].convective_flux_computed_in_bc) continue;
                    if ((k == kmax+1) && bc[Face.top].convective_flux_computed_in_bc) continue;
                    //
                    if ((k == kmin) && !(bc[Face.bottom].ghost_cell_data_available)) {
                        compute_flux_at_left_wall(Rght, IFace, myConfig, omegaz);
                    } else if ((k == kmax+1) && !(bc[Face.top].ghost_cell_data_available)) {
                        compute_flux_at_right_wall(Lft, IFace, myConfig, omegaz);
                    } else {
                        compute_interface_flux(Lft, Rght, IFace, myConfig, omegaz);
                    }
                } // for k 
            } // j loop
        } // i loop
        return;
    } // end convective_flux_phase0()

    @nogc
    override void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0)
    // Compute the flux from data on either-side of the interface.
    // For the structured-grid block, there is nothing to do.
    // The unstructured-grid block needs to work in two phases.
    {
        return;
    }

    @nogc void copy_current_corner_coords()
    {
        if (myConfig.dimensions == 2) {
            FVVertex vtx00 = get_vtx!()(imin,jmin);
            corner_coords[0] = vtx00.pos[0].x.re;
            corner_coords[1] = vtx00.pos[0].y.re;
            corner_coords[2] = vtx00.pos[0].z.re;
            FVVertex vtx10 = get_vtx!()(imax+1,jmin);
            corner_coords[3] = vtx10.pos[0].x.re;
            corner_coords[4] = vtx10.pos[0].y.re;
            corner_coords[5] = vtx10.pos[0].z.re;
            FVVertex vtx11 = get_vtx!()(imax+1,jmax+1);
            corner_coords[6] = vtx11.pos[0].x.re;
            corner_coords[7] = vtx11.pos[0].y.re;
            corner_coords[8] = vtx11.pos[0].z.re;
            FVVertex vtx01 = get_vtx!()(imin,jmax+1);
            corner_coords[9] = vtx01.pos[0].x.re;
            corner_coords[10] = vtx01.pos[0].y.re;
            corner_coords[11] = vtx01.pos[0].z.re;
            // In 2D, the upper layer get the same values.
            corner_coords[12] = vtx00.pos[0].x.re;
            corner_coords[13] = vtx00.pos[0].y.re;
            corner_coords[14] = vtx00.pos[0].z.re;
            //
            corner_coords[15] = vtx10.pos[0].x.re;
            corner_coords[16] = vtx10.pos[0].y.re;
            corner_coords[17] = vtx10.pos[0].z.re;
            //
            corner_coords[18] = vtx11.pos[0].x.re;
            corner_coords[19] = vtx11.pos[0].y.re;
            corner_coords[20] = vtx11.pos[0].z.re;
            //
            corner_coords[21] = vtx01.pos[0].x.re;
            corner_coords[22] = vtx01.pos[0].y.re;
            corner_coords[23] = vtx01.pos[0].z.re;
        } else {
            FVVertex vtx000 = get_vtx!()(imin,jmin,kmin);
            corner_coords[0] = vtx000.pos[0].x.re;
            corner_coords[1] = vtx000.pos[0].y.re;
            corner_coords[2] = vtx000.pos[0].z.re;
            FVVertex vtx100 = get_vtx!()(imax+1,jmin,kmin);
            corner_coords[3] = vtx100.pos[0].x.re;
            corner_coords[4] = vtx100.pos[0].y.re;
            corner_coords[5] = vtx100.pos[0].z.re;
            FVVertex vtx110 = get_vtx!()(imax+1,jmax+1,kmin);
            corner_coords[6] = vtx110.pos[0].x.re;
            corner_coords[7] = vtx110.pos[0].y.re;
            corner_coords[8] = vtx110.pos[0].z.re;
            FVVertex vtx010 = get_vtx!()(imin,jmax+1,kmin);
            corner_coords[9] = vtx010.pos[0].x.re;
            corner_coords[10] = vtx010.pos[0].y.re;
            corner_coords[11] = vtx010.pos[0].z.re;
            FVVertex vtx001 = get_vtx!()(imin,jmin,kmax+1);
            corner_coords[12] = vtx001.pos[0].x.re;
            corner_coords[13] = vtx001.pos[0].y.re;
            corner_coords[14] = vtx001.pos[0].z.re;
            FVVertex vtx101 = get_vtx!()(imax+1,jmin,kmax+1);
            corner_coords[15] = vtx101.pos[0].x.re;
            corner_coords[16] = vtx101.pos[0].y.re;
            corner_coords[17] = vtx101.pos[0].z.re;
            FVVertex vtx111 = get_vtx!()(imax+1,jmax+1,kmax+1);
            corner_coords[18] = vtx111.pos[0].x.re;
            corner_coords[19] = vtx111.pos[0].y.re;
            corner_coords[20] = vtx111.pos[0].z.re;
            FVVertex vtx011 = get_vtx!()(imin,jmax+1,kmax+1);
            corner_coords[21] = vtx011.pos[0].x.re;
            corner_coords[22] = vtx011.pos[0].y.re;
            corner_coords[23] = vtx011.pos[0].z.re;
        }
    } // end copy_current_corner_coords()
    
    @nogc void set_current_corner_coords_to_infinity()
    {
        foreach(ref double cc; corner_coords) { cc = double.infinity; } 
    }
} // end class SFluidBlock


/** Indexing of the data in 2D.
 *
 * \verbatim
 * The following figure shows cell [i,j] and its associated
 * vertices and faces. 
 * (New arrangement, planned August 2006, implemented Nov 2006)
 *
 *
 *
 *     Vertex 3         North face           Vertex 2 
 *   vtx[i,j+1]         ifj[i,j+1]           vtx[i+1,j+1]
 *             +--------------x--------------+
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *   West      |         cell center         |  East 
 *   face      |          ctr[i,j]           |  face
 *   ifi[i,j]  x              o              x  ifi[i+1,j]
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             +--------------x--------------+
 *     Vertex 0           South face         Vertex 1
 *     vtx[i,j]           ifj[i,j]           vtx[i+1,j]
 *
 *
 * Thus...
 * ----
 * Active cells are indexed as ctr[i][i], where
 * imin <= i <= imax, jmin <= j <= jmax.
 *
 * Active east-facing interfaces are indexed as ifi[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax.
 *
 * Active north-facing interfaces are indexed as ifj[i][j], where
 * imin <= i <= imax, jmin <= j <= jmax+1.
 *
 * Active vertices are indexed as vtx[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax+1.
 *
 * Space for ghost cells is available outside these ranges.
 *
 * Indexing for the 3D data -- see page 8 in 3D CFD workbook
 * \endverbatim
 */

