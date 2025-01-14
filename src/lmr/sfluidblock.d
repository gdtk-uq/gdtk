// sfluidblock.d
// Class for structured blocks of cells, for use within Eilmer4.
// This is the "classic" block within the mbcns/Eilmer series
// of flow simulation codes.

// Peter J. 2014-07-20 first cut.

module lmr.sfluidblock;

import std.algorithm;
import std.array;
import std.conv;
import std.file;
import std.format;
import std.json;
import std.math;
import std.range;
import std.stdio;
import std.string;

import gzip;

import gas.luagas_model;
import gas;
import geom;
import geom.luawrap.luasgrid;
import kinetics;
import nm.number;
import ntypes.complex;
import util.json_helper;
import util.lua;
import util.lua_service;

import lmr.bc;
import lmr.block;
import lmr.conservedquantities;
import lmr.coredata;
import lmr.flowgradients;
import lmr.flowstate;
import lmr.fluidblock;
import lmr.fluidfvcell;
import lmr.fluxcalc;
import lmr.fvinterface;
import lmr.fvvertex;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.grid_motion;
import lmr.grid_motion_udf;
import lmr.lua_helper;
import lmr.luawrap.luaflowstate;
import lmr.onedinterp;

// EPSILON parameter for numerical differentiation of flux jacobian
// Value used based on Vanden and Orkwis (1996), AIAA J. 34:6 pp. 1125-1129
immutable double EPSILON = 1.0e-8;
immutable double ESSENTIALLY_ZERO = 1.0e-15;

@nogc
size_t[3] cell_index_to_logical_coordinates(size_t gid, size_t nic, size_t njc)
{
    pragma(inline, true);
    size_t[3] ijk;
    size_t slabDim = njc * nic;
    size_t k = gid / slabDim;
    size_t j = (gid - k*slabDim) / nic;
    size_t i = gid - k*slabDim - j*nic;
    return [i, j, k];
}

class SFluidBlock: FluidBlock {
public:
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
    size_t[] micell, mjcell, mkcell; // locations of monitor cells
    //
    size_t nic, njc, nkc;
    size_t niv, njv, nkv;
    StructuredGrid grid; // for reading and writing
    //
    // A place to store coordinates of the corner vertices.
    // For a moving-grid simulation these will be kept up to date
    // and communicated to user-defined Lua functions via infoFluidBlock.
    double[24] corner_coords;
    //
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    OneDInterpolator one_d;
    ConservedQuantities flux;

public:
    this(int blk_id, size_t nicell, size_t njcell, size_t nkcell, string label)
    {
        super(blk_id, Grid_t.structured_grid, nicell*njcell*nkcell,
              GlobalConfig.n_ghost_cell_layers, label);
        this.n_ghost_cell_layers = GlobalConfig.n_ghost_cell_layers;
        nic = nicell;
        njc = njcell;
        nkc = nkcell;
        // Fill in other data sizes.
        niv = nic + 1;
        njv = njc + 1;
        if (GlobalConfig.dimensions == 2) {
            // In 2D simulations, the k range is from 0 to 0 for the
            // storage of cells, vertices and relevant faces.
            if (nkc != 1) { throw new Error("Inconsistent dimensions and value for nkc."); }
            nkv = 1;
        } else {
            // In 3D simulations the k index is just like the i and j indices.
            nkv = nkc + 1;
        }
    } // end constructor

    this(int blk_id, JSONValue json_data)
    {
        size_t nicell = getJSONint(json_data, "nic", 0);
        size_t njcell = getJSONint(json_data, "njc", 0);
        size_t nkcell = getJSONint(json_data, "nkc", 0);
        label = getJSONstring(json_data, "label", "");
        this(blk_id, nicell, njcell, nkcell, label);
        active = getJSONbool(json_data, "active", true);
        omegaz = getJSONdouble(json_data, "omegaz", 0.0);
        may_be_turbulent = getJSONbool(json_data, "may_be_turbulent", true);
    } // end constructor from json

    this(lua_State* L)
    // Generate a new block and fill it with information from a Lua interpreter state.
    // This particular constructor is only used in the prep stage.
    {
        auto grid = checkStructuredGrid(L, 2);
        double omegaz = luaL_checknumber(L, 5);
        nic = grid.niv - 1;
        njc = grid.njv - 1;
        nkc = grid.nkv - 1;
        if (GlobalConfig.dimensions == 2) nkc = 1;
        ncells = nic*njc*nkc;
        super(-1, "nil");
        grid_type = Grid_t.structured_grid;
        ncells_expected = ncells;
        n_ghost_cell_layers = 0;
        // do we need a lighter weight config instantiation process?
        // where could we get a pre-built instance from, GlobalConfig?
        myConfig = new LocalConfig(-1);
        myConfig.init_gas_model_bits();
        cells.length = ncells; // not defined yet

        // We don't need the full celldata for the prep stage, just the flowstates
        celldata.flowstates.reserve(ncells);

        bool lua_fs = false;
        FlowState* myfs;
        // check where our flowstate is coming from
        if ( isObjType(L, 3, "number") ) {
            myfs = checkFlowState(L, 3);
            lua_settop(L, 0); // clear stack
        } else if (lua_isfunction(L, 3)) {
            lua_fs = true;
        }
        Vector3 pos;
        number volume, iLen, jLen, kLen;
        size_t[3] ijk;
        size_t i, j, k;
        foreach(cell_idx, ref cell; cells) {
            // get the cell index in terms of structured grid indices
            ijk = cell_index_to_logical_coordinates(cell_idx, nic, njc);
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
            } else if (GlobalConfig.dimensions == 3) {
                Vector3 p001 = *grid[i,j,k+1];
                Vector3 p101 = *grid[i+1,j,k+1];
                Vector3 p111 = *grid[i+1,j+1,k+1];
                Vector3 p011 = *grid[i,j+1,k+1];
                hex_cell_properties(p000, p100, p110, p010, p001, p101, p111, p011,
                                    GlobalConfig.true_centroids, pos, volume, iLen, jLen, kLen);
            } else {
                throw new Exception("GlobalConfig.dimensions not 2 or 3.");
            }
            if (omegaz != 0.0) { into_rotating_frame(myfs.vel, pos, omegaz); }
            if (lua_fs) {
                // Now grab flow state via Lua function call.
                // If the block is in a rotating frame with omegaz != 0.0,
                // we presume that the Lua function will provide the velocity
                // components relative to the rotating frame.
                lua_pushvalue(L, 3);
                lua_pushnumber(L, pos.x);
                lua_pushnumber(L, pos.y);
                lua_pushnumber(L, pos.z);
                if (lua_pcall(L, 3, 1, 0) != 0) {
                    string errMsg = "Error in Lua function call for setting FlowState\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    luaL_error(L, errMsg.toStringz);
                }
                if (lua_istable(L, -1)) {
                    myfs = makeFlowStateFromTable(L, lua_gettop(L));
                } else {
                    myfs = checkFlowState(L, -1);
                }
                if (!myfs) {
                    string errMsg = "Error in from Lua function call for setting FlowState\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    errMsg ~= "The returned object is not a proper _FlowState handle or table.";
                    luaL_error(L, errMsg.toStringz);
                }
            }
            // make the cell
            celldata.flowstates ~= *myfs; // Copy the myfs flowstate into the celldata structure
            cell = new FluidFVCell(myConfig, pos, &celldata, volume, to!int(cell_idx));
        }
        if (lua_fs) { lua_settop(L, 0); }
    } // end constructor from Lua state

    override JSONValue get_header()
    // return information in JSON format that describes this block
    {
        JSONValue header = ["structured": true];
        header["nic"] = to!int(nic);
        header["njc"] = to!int(njc);
        header["nkc"] = to!int(nkc);
        return header;
    }

    override void init_workspace()
    {
        super.init_workspace();
        // Workspace for flux_calc method.
        one_d = new OneDInterpolator(dedicatedConfig[id]);
    }

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
        // Lua interpreter for the block.
        // It will be available for computing user-defined source terms.
        lua_pushglobaltable(myL);
        registerGasModel(myL);
        pushObj!(GasModel, GasModelMT)(myL, dedicatedConfig[id].gmodel);
        lua_setglobal(myL, "gmodel");
        lua_pushinteger(myL, dedicatedConfig[id].n_species);
        lua_setglobal(myL, "n_species");
        lua_pushinteger(myL, dedicatedConfig[id].n_modes);
        lua_setglobal(myL, "n_modes");
        lua_pushinteger(myL, nic); lua_setglobal(myL, "nicell");
        lua_pushinteger(myL, njc); lua_setglobal(myL, "njcell");
        lua_pushinteger(myL, nkc); lua_setglobal(myL, "nkcell");
        lua_pushinteger(myL, n_ghost_cell_layers); lua_setglobal(myL, "n_ghost_cell_layers");
        // Although we make the helper functions available within
        // the block-specific Lua interpreter, we should use
        // those functions only in the context of the master thread.
        setSampleHelperFunctions(myL);
        // Generally, the sampleFluidCell function can be expected to work only in serial mode.
        // Once it is called from a thread, other than the main thread, it may not
        // have access to properly initialized data for any other block.
        setGridMotionHelperFunctions(myL);
    } // end init_lua_globals()

    override void init_boundary_conditions(JSONValue json_data)
    // Initialize boundary conditions after the blocks are constructed,
    // because we want access to the full collection of valid block references.
    {
        foreach (boundary; 0 .. (myConfig.dimensions == 3 ? 6 : 4)) {
            string json_key = "boundary_" ~ face_name[boundary];
            auto bc_json_data = json_data[json_key];
            bc ~= make_BC_from_json(bc_json_data, id, boundary);
        }
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
        repr ~= ", may_be_turbulent=" ~ to!string(may_be_turbulent);
        repr ~= ", nic=" ~ to!string(nic);
        repr ~= ", njc=" ~ to!string(njc);
        repr ~= ", nkc=" ~ to!string(nkc);
        repr ~= ", n_ghost_cell_layers=" ~ to!string(n_ghost_cell_layers);
        repr ~= ", \n    bc=["~ face_name[0] ~ "=" ~ to!string(bc[0]);
        foreach (i; 1 .. (myConfig.dimensions == 3 ? 6 : 4)) {
            repr ~= ",\n        " ~ face_name[i] ~ "=" ~ to!string(bc[i]);
        }
        repr ~= "\n       ]"; // end bc list
        repr ~= ")";
        return to!string(repr);
    } // end toString()

    @nogc
    size_t cell_index(size_t i, size_t j, size_t k=0) const
    // The i,j,k indices into the hypothetical block of active cells
    // map to the index into the single-dimensional cells array
    // that is held in the FluidBlock base class.
    {
        pragma(inline, true);
        debug {
            if (!(i < nic && j < njc && k < nkc)) {
                writefln("cell_index[%d,%d,%d] from [%d,%d,%d] in blk.id %d", i, j, k, nic, njc, nkc, id);
            }
        }
        assert(i < nic && j < njc && k < nkc, "Index out of bounds.");
        return (k*njc + j)*nic + i;
    }

    @nogc
    size_t[3] to_ijk_indices_for_cell(size_t gid) const
    {
        pragma(inline, true);
        return cell_index_to_logical_coordinates(gid, nic, njc);
    }

    @nogc
    size_t vertex_index(size_t i, size_t j, size_t k=0) const
    // The i,j,k indices into the hypothetical block of vertices
    // map to the index into the single-dimensional vertices array
    // that is held in the FluidBlock base class.
    {
        pragma(inline, true);
        debug {
            if (!(i < niv && j < njv && k < nkv)) {
                writefln("vertex_index[%d,%d,%d] from [%d,%d,%d] in blk.id %d", i, j, k, niv, njv, nkv, id);
            }
        }
        assert(i < niv && j < njv && k < nkv, "Index out of bounds.");
        return (k*njv + j)*niv + i;
    }

    @nogc
    size_t ifi_index(size_t i, size_t j, size_t k=0) const
    // The i,j,k indices into the hypothetical block of i-faces
    // map to the index into the single-dimensional vertices array
    // that is held in the FluidBlock base class.
    {
        pragma(inline, true);
        debug {
            if (!(i < niv && j < njc && k < nkc)) {
                writefln("ifi_index[%d,%d,%d] from [%d,%d,%d] in blk.id %d", i, j, k, niv, njc, nkc, id);
            }
        }
        assert(i < niv && j < njc && k < nkc, "Index out of bounds.");
        return (k*njc + j)*niv + i;
    }

    @nogc
    size_t ifj_index(size_t i, size_t j, size_t k=0) const
    // The i,j,k indices into the hypothetical block of j-faces
    // map to the index into the single-dimensional vertices array
    // that is held in the FluidBlock base class.
    {
        pragma(inline, true);
        debug {
            if (!(i < nic && j < njv && k < nkc)) {
                writefln("ifj_index[%d,%d,%d] from [%d,%d,%d] in blk.id %d", i, j, k, nic, njv, nkc, id);
            }
        }
        assert(i < nic && j < njv && k < nkc, "Index out of bounds.");
        size_t nifaces = niv*njc*nkc;
        return (k*njv + j)*nic + i + nifaces;
    }

    @nogc
    size_t ifk_index(size_t i, size_t j, size_t k=0) const
    // The i,j,k indices into the hypothetical block of k-faces
    // map to the index into the single-dimensional vertices array
    // that is held in the FluidBlock base class.
    {
        pragma(inline, true);
        debug {
            if (!(i < nic && j < njc && k < nkv)) {
                writefln("ifj_index[%d,%d,%d] from [%d,%d,%d] in blk.id %d", i, j, k, nic, njc, nkv, id);
            }
        }
        assert(i < nic && j < njc && k < nkv, "Index out of bounds.");
        size_t nifaces = niv*njc*nkc;
        size_t njfaces = nic*njv*nkc;
        return (k*njc + j)*nic + i + nifaces + njfaces;
    }

    @nogc ref FluidFVCell get_cell(size_t i, size_t j, size_t k=0)
    {
        pragma(inline, true);
        return cells[cell_index(i,j,k)];
    }
    @nogc ref FVVertex get_vtx(size_t i, size_t j, size_t k=0)
    {
        pragma(inline, true);
        return vertices[vertex_index(i,j,k)];
    }
    @nogc ref FVInterface get_ifi(size_t i, size_t j, size_t k=0)
    {
        pragma(inline, true);
        return faces[ifi_index(i,j,k)];
    }
    @nogc ref FVInterface get_ifj(size_t i, size_t j, size_t k=0)
    {
        pragma(inline, true);
        return faces[ifj_index(i,j,k)];
    }
    @nogc ref FVInterface get_ifk(size_t i, size_t j, size_t k=0)
    {
        pragma(inline, true);
        return faces[ifk_index(i,j,k)];
    }

    @nogc
    override void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
        grid.find_enclosing_cell(p, indx, found); // delegate to the grid object
    }

    override void init_grid_and_flow_arrays(string gridFileName)
    {
        if (myConfig.verbosity_level > 1) { writeln("init_grid_and_flow_arrays(): Start block ", id); }
        //
        bool lsq_workspace_at_faces = (myConfig.viscous)
            && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces);
        bool lsq_workspace_at_vertices = (myConfig.viscous)
            && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.vertices);
        bool lsq_workspace_at_cells = (myConfig.viscous)
            && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.cells);
        try {
            // First allocate densely packed core data structures
            auto gmodel = myConfig.gmodel;
            size_t nturb = myConfig.turb_model.nturb;
            size_t neq = myConfig.cqi.n;
            size_t nftl = myConfig.n_flow_time_levels;

            flux = new_ConservedQuantities(neq);
            ncells = nic*njc*nkc;
            size_t nifaces = niv*njc*nkc;
            size_t njfaces = nic*njv*nkc;
            size_t nkfaces = nic*njc*nkv;
            nfaces = nifaces + njfaces;
            if (myConfig.dimensions == 3) nfaces += nkfaces;
            size_t nbfaces = nkc*nic*2 + nkc*njc*2;
            if (myConfig.dimensions == 3) nbfaces += 2*njc*nic;
            size_t nghost = (2*njc*nkc + 2*nic*nkc)*n_ghost_cell_layers;
            if (myConfig.dimensions == 3) nghost += 2*nic*njc*n_ghost_cell_layers;

            allocate_dense_celldata(ncells, nghost, neq, nftl);
            allocate_dense_facedata(nfaces, nbfaces, neq, nftl);

            // Create the interior cell, vertex and interface objects for the block.
            foreach (n; 0 .. nic*njc*nkc) {
                cells ~= new FluidFVCell(myConfig, &celldata, to!int(n));
            }
            foreach (n; 0 .. niv*njv*nkv) {
                vertices ~= new FVVertex(myConfig, lsq_workspace_at_vertices);
            }
            // First, ifi faces.
            int fid = 0;
            foreach (n; 0 .. niv*njc*nkc) {
                faces ~= new FVInterface(myConfig, IndexDirection.i, &facedata, fid); fid++;
            }
            // Second, ifj faces.
            foreach (n; 0 .. nic*njv*nkc) {
                faces ~= new FVInterface(myConfig, IndexDirection.j, &facedata, fid); fid++;
            }
            // Third, maybe, ifk faces.
            if (myConfig.dimensions == 3) {
                foreach (n; 0 .. nic*njc*nkv) {
                    faces ~= new FVInterface(myConfig, IndexDirection.k, &facedata, fid); fid++;
                }
            }
            // Now, construct the ghost cells, attaching them to the boundary faces.
            int cell_id = to!int(nic*njc*nkc);
            // North and South boundaries.
            if (bc[Face.north].ghost_cell_data_available) {
                bool shock_fitting = (bc[Face.north].type == "inflow_shock_fitting");
                foreach (k; 0 .. nkc) {
                    foreach (i; 0 .. nic) {
                        auto f = get_ifj(i, njc, k);
                        foreach (n; 0 .. n_ghost_cell_layers) {
                            auto c = new FluidFVCell(myConfig, &celldata, cell_id);
                            ++cell_id; c.is_ghost = true;
                            if (shock_fitting) c.is_in_shock_fitting_boundary = true;
                            f.right_cells ~= c;
                        }
                    }
                }
            }
            if (bc[Face.south].ghost_cell_data_available) {
                bool shock_fitting = (bc[Face.south].type == "inflow_shock_fitting");
                foreach (k; 0 .. nkc) {
                    foreach (i; 0 .. nic) {
                        auto f = get_ifj(i, 0, k);
                        foreach (n; 0 .. n_ghost_cell_layers) {
                            auto c = new FluidFVCell(myConfig, &celldata, cell_id);
                            ++cell_id; c.is_ghost = true;
                            if (shock_fitting) c.is_in_shock_fitting_boundary = true;
                            f.left_cells ~= c;
                        }
                    }
                }
            }
            if (bc[Face.east].ghost_cell_data_available) {
                bool shock_fitting = (bc[Face.east].type == "inflow_shock_fitting");
                foreach (k; 0 .. nkc) {
                    foreach (j; 0 .. njc) {
                        auto f = get_ifi(nic, j, k);
                        foreach (n; 0 .. n_ghost_cell_layers) {
                            auto c = new FluidFVCell(myConfig, &celldata, cell_id);
                            ++cell_id; c.is_ghost = true;
                            if (shock_fitting) c.is_in_shock_fitting_boundary = true;
                            f.right_cells ~= c;
                        }
                    }
                }
            }
            if (bc[Face.west].ghost_cell_data_available) {
                bool shock_fitting = (bc[Face.west].type == "inflow_shock_fitting");
                foreach (k; 0 .. nkc) {
                    foreach (j; 0 .. njc) {
                        auto f = get_ifi(0, j, k);
                        foreach (n; 0 .. n_ghost_cell_layers) {
                            auto c = new FluidFVCell(myConfig, &celldata, cell_id);
                            ++cell_id; c.is_ghost = true;
                            if (shock_fitting) c.is_in_shock_fitting_boundary = true;
                            f.left_cells ~= c;
                        }
                    }
                }
            }
            if (myConfig.dimensions == 3) {
                if (bc[Face.top].ghost_cell_data_available) {
                    bool shock_fitting = (bc[Face.top].type == "inflow_shock_fitting");
                    foreach (j; 0 .. njc) {
                        foreach (i; 0 .. nic) {
                            auto f = get_ifk(i, j, nkc);
                            foreach (n; 0 .. n_ghost_cell_layers) {
                                auto c = new FluidFVCell(myConfig, &celldata, cell_id);
                                ++cell_id; c.is_ghost = true;
                                if (shock_fitting) c.is_in_shock_fitting_boundary = true;
                                f.right_cells ~= c;
                            }
                        }
                    }
                }
                if (bc[Face.bottom].ghost_cell_data_available) {
                bool shock_fitting = (bc[Face.bottom].type == "inflow_shock_fitting");
                    foreach (j; 0 .. njc) {
                        foreach (i; 0 .. nic) {
                            auto f = get_ifk(i, j, 0);
                            foreach (n; 0 .. n_ghost_cell_layers) {
                                auto c = new FluidFVCell(myConfig, &celldata, cell_id);
                                ++cell_id; c.is_ghost = true;
                                if (shock_fitting) c.is_in_shock_fitting_boundary = true;
                                f.left_cells ~= c;
                            }
                        }
                    }
                }
            } // end if (myConfig.dimensions == 3)
        } catch (Exception e) {
            writeln("Failed while assembling block arrays.");
            writefln("nic=%d njc=%d nkc=%d", nic, njc, nkc);
            writeln("Probably ran out of memory.");
            writeln("Be a little less ambitious and try a smaller grid next time.");
            writefln("System message: %s", e.msg);
            throw new FlowSolverException("SFluidBlock.init_grid_and_flow_arrays() failed.");
        }
        // Setup dense arrays for fast checks of left or right data available (NNG)
        // North
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                size_t idx = ifj_index(i, njc, k);
                if (!startsWith(bc[Face.north].type, "exchange_"))
                    facedata.left_interior_only[idx] = true;
            }
        }
        // South
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                size_t idx = ifj_index(i, 0, k);
                if (!startsWith(bc[Face.south].type, "exchange_"))
                    facedata.right_interior_only[idx] = true;
            }
        }
        // East
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                size_t idx = ifi_index(nic, j, k);
                if (!startsWith(bc[Face.east].type, "exchange_"))
                    facedata.left_interior_only[idx] = true;
            }
        }
        // West
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                size_t idx = ifi_index(0, j, k);
                if (!startsWith(bc[Face.west].type, "exchange_"))
                    facedata.right_interior_only[idx] = true;
            }
        }
        if (myConfig.dimensions == 3) {
            // Top
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t idx = ifk_index(i, j, nkc);
                    if (!startsWith(bc[Face.top].type, "exchange_"))
                        facedata.left_interior_only[idx] = true;
                }
            }
            // Bottom
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t idx = ifk_index(i, j, 0);
                    if (!startsWith(bc[Face.bottom].type, "exchange_"))
                        facedata.right_interior_only[idx] = true;
                }
            }
        } // end if (myConfig.dimensions == 3)
        //
        // Densified structures for the onedinterp routines in the structured formulation
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. niv) {
                    size_t fid = ifi_index(i,j,k);
                    if (i>2)     facedata.stencil_idxs[fid].L2 = cell_index(i-3, j, k);
                    if (i>1)     facedata.stencil_idxs[fid].L1 = cell_index(i-2, j, k);
                    if (i>0)     facedata.stencil_idxs[fid].L0 = cell_index(i-1, j, k);
                    if (i+1<niv) facedata.stencil_idxs[fid].R0 = cell_index(i+0, j, k);
                    if (i+2<niv) facedata.stencil_idxs[fid].R1 = cell_index(i+1, j, k);
                    if (i+3<niv) facedata.stencil_idxs[fid].R2 = cell_index(i+2, j, k);
                } // i loop
            } // j loop
        } // k loop
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njv) {
                    size_t fid = ifj_index(i,j,k);
                    if (j>2)     facedata.stencil_idxs[fid].L2 = cell_index(i, j-3, k);
                    if (j>1)     facedata.stencil_idxs[fid].L1 = cell_index(i, j-2, k);
                    if (j>0)     facedata.stencil_idxs[fid].L0 = cell_index(i, j-1, k);
                    if (j+1<njv) facedata.stencil_idxs[fid].R0 = cell_index(i, j+0, k);
                    if (j+2<njv) facedata.stencil_idxs[fid].R1 = cell_index(i, j+1, k);
                    if (j+3<njv) facedata.stencil_idxs[fid].R2 = cell_index(i, j+2, k);
                } // j loop
            } // i loop
        } // k loop
        if (myConfig.dimensions == 3) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    foreach (k; 0 .. nkv) {
                        size_t fid = ifk_index(i,j,k);
                        if (k>2)     facedata.stencil_idxs[fid].L2 = cell_index(i, j, k-3);
                        if (k>1)     facedata.stencil_idxs[fid].L1 = cell_index(i, j, k-2);
                        if (k>0)     facedata.stencil_idxs[fid].L0 = cell_index(i, j, k-1);
                        if (k+1<nkv) facedata.stencil_idxs[fid].R0 = cell_index(i, j, k+0);
                        if (k+2<nkv) facedata.stencil_idxs[fid].R1 = cell_index(i, j, k+1);
                        if (k+3<nkv) facedata.stencil_idxs[fid].R2 = cell_index(i, j, k+2);
                    } // k loop
                } // i loop
            } // j loop
        } // end 3d
        //
        // Now, loop through the ghost cells and attach their indices to the face stencil objects
        size_t cid = nic*njc*nkc;
        // North
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                size_t fid0 = ifj_index(i, njc, k);
                size_t fid1 = ifj_index(i, njc-1, k);
                size_t fid2 = ifj_index(i, njc-2, k);
                facedata.stencil_idxs[fid0].R0 = cid;
                facedata.stencil_idxs[fid1].R1 = cid;
                facedata.stencil_idxs[fid2].R2 = cid;
                cid++;
                facedata.stencil_idxs[fid0].R1 = cid;
                facedata.stencil_idxs[fid1].R2 = cid;
                cid++;
                if (n_ghost_cell_layers>2) {
                    facedata.stencil_idxs[fid0].R2 = cid;
                    cid++;
                }
            }
        }
        // South
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                size_t fid0 = ifj_index(i, 0, k);
                size_t fid1 = ifj_index(i, 1, k);
                size_t fid2 = ifj_index(i, 2, k);
                facedata.stencil_idxs[fid0].L0 = cid;
                facedata.stencil_idxs[fid1].L1 = cid;
                facedata.stencil_idxs[fid2].L2 = cid;
                cid++;
                facedata.stencil_idxs[fid0].L1 = cid;
                facedata.stencil_idxs[fid1].L2 = cid;
                cid++;
                if (n_ghost_cell_layers>2) {
                    facedata.stencil_idxs[fid0].L2 = cid;
                    cid++;
                }
            }
        }
        // East
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                size_t fid0 = ifi_index(nic, j, k);
                size_t fid1 = ifi_index(nic-1, j, k);
                size_t fid2 = ifi_index(nic-2, j, k);
                facedata.stencil_idxs[fid0].R0 = cid;
                facedata.stencil_idxs[fid1].R1 = cid;
                facedata.stencil_idxs[fid2].R2 = cid;
                cid++;
                facedata.stencil_idxs[fid0].R1 = cid;
                facedata.stencil_idxs[fid1].R2 = cid;
                cid++;
                if (n_ghost_cell_layers>2) {
                    facedata.stencil_idxs[fid0].R2 = cid;
                    cid++;
                }
            }
        }
        // West
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                size_t fid0 = ifi_index(0, j, k);
                size_t fid1 = ifi_index(1, j, k);
                size_t fid2 = ifi_index(2, j, k);
                facedata.stencil_idxs[fid0].L0 = cid;
                facedata.stencil_idxs[fid1].L1 = cid;
                facedata.stencil_idxs[fid2].L2 = cid;
                cid++;
                facedata.stencil_idxs[fid0].L1 = cid;
                facedata.stencil_idxs[fid1].L2 = cid;
                cid++;
                if (n_ghost_cell_layers>2) {
                    facedata.stencil_idxs[fid0].L2 = cid;
                    cid++;
                }
            }
        }
        if (myConfig.dimensions == 3) {
            // Top
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t fid0 = ifk_index(i, j, nkc);
                    size_t fid1 = ifk_index(i, j, nkc-1);
                    size_t fid2 = ifk_index(i, j, nkc-2);
                    facedata.stencil_idxs[fid0].R0 = cid;
                    facedata.stencil_idxs[fid1].R1 = cid;
                    facedata.stencil_idxs[fid2].R2 = cid;
                    cid++;
                    facedata.stencil_idxs[fid0].R1 = cid;
                    facedata.stencil_idxs[fid1].R2 = cid;
                    cid++;
                    if (n_ghost_cell_layers>2) {
                        facedata.stencil_idxs[fid0].R2 = cid;
                        cid++;
                    }
                }
            }
            // Bottom
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t fid0 = ifk_index(i, j, 0);
                    size_t fid1 = ifk_index(i, j, 1);
                    size_t fid2 = ifk_index(i, j, 2);
                    facedata.stencil_idxs[fid0].L0 = cid;
                    facedata.stencil_idxs[fid1].L1 = cid;
                    facedata.stencil_idxs[fid2].L2 = cid;
                    cid++;
                    facedata.stencil_idxs[fid0].L1 = cid;
                    facedata.stencil_idxs[fid1].L2 = cid;
                    cid++;
                    if (n_ghost_cell_layers>2) {
                        facedata.stencil_idxs[fid0].L2 = cid;
                        cid++;
                    }
                }
            }
        } // end if (myConfig.dimensions == 3)

        // TODO: Redundant storage. But we use these in some areas of the code that are structure agnostic
        foreach(idx; 0 .. nfaces){
            facedata.f2c[idx].left  = facedata.stencil_idxs[idx].L0;
            facedata.f2c[idx].right = facedata.stencil_idxs[idx].R0;
        }
        // Now that all of the cells and faces are constructed,
        // We want to store the local structure of the grid
        // in the array of references stored in the faces.
        //   cL[2] cL[1] cL[0] f cR[0] cR[1] cR[2]
        //     o     o     o   |   o     o     o
        // Note that the boundary faces already have ghost-cells attached
        // if the boundary-condition accommodates them.
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 1 .. niv) {
                    // Work from left to right, attaching cells to left of each face.
                    auto fL = get_ifi(i,j,k);
                    fL.left_cells ~= get_cell(i-1,j,k);
                    auto fL1 = get_ifi(i-1,j,k);
                    foreach (n; 0 .. n_ghost_cell_layers-1) {
                        if (fL1.left_cells.length > n) { fL.left_cells ~= fL1.left_cells[n]; }
                    }
                    // Work from right to left, attaching cells to right of each face.
                    size_t iR = niv-1-i;
                    auto fR = get_ifi(iR,j,k);
                    fR.right_cells ~= get_cell(iR,j,k);
                    auto fR1 = get_ifi(iR+1,j,k);
                    foreach (n; 0 .. n_ghost_cell_layers-1) {
                        if (fR1.right_cells.length > n) { fR.right_cells ~= fR1.right_cells[n]; }
                    }
                } // i loop
            } // j loop
        } // k loop
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                foreach (j; 1 .. njv) {
                    // Work from left to right, attaching cells to left of each face.
                    auto fL = get_ifj(i,j,k);
                    fL.left_cells ~= get_cell(i,j-1,k);
                    auto fL1 = get_ifj(i,j-1,k);
                    foreach (n; 0 .. n_ghost_cell_layers-1) {
                        if (fL1.left_cells.length > n) { fL.left_cells ~= fL1.left_cells[n]; }
                    }
                    // Work from right to left, attaching cells to right of each face.
                    size_t jR = njv-1-j;
                    auto fR = get_ifj(i,jR,k);
                    fR.right_cells ~= get_cell(i,jR,k);
                    auto fR1 = get_ifj(i,jR+1,k);
                    foreach (n; 0 .. n_ghost_cell_layers-1) {
                        if (fR1.right_cells.length > n) { fR.right_cells ~= fR1.right_cells[n]; }
                    }
                } // j loop
            } // i loop
        } // k loop
        if (myConfig.dimensions == 3) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    foreach (k; 1 .. nkv) {
                        // Work from left to right, attaching cells to left of each face.
                        auto fL = get_ifk(i,j,k);
                        fL.left_cells ~= get_cell(i,j,k-1);
                        auto fL1 = get_ifk(i,j,k-1);
                        foreach (n; 0 .. n_ghost_cell_layers-1) {
                            if (fL1.left_cells.length > n) { fL.left_cells ~= fL1.left_cells[n]; }
                        }
                        // Work from right to left, attaching cells to right of each face.
                        size_t kR = nkv-1-k;
                        auto fR = get_ifk(i,j,kR);
                        fR.right_cells ~= get_cell(i,j,kR);
                        auto fR1 = get_ifk(i,j,kR+1);
                        foreach (n; 0 .. n_ghost_cell_layers-1) {
                            if (fR1.right_cells.length > n) { fR.right_cells ~= fR1.right_cells[n]; }
                        }
                    } //k loop
                } // i loop
            } // j loop
        } // end if (myConfig.dimensions == 3)
        // Make the cell, vertex, and face id value consistent with the index in the array.
        // We will depend on this equality in other parts of the flow solver.
        // We also note that these cells are interior to the block (i.e. not ghost cells)
        foreach (i, c; cells) {
            c.contains_flow_data = true;
            c.is_interior_to_domain = true;
        }
        foreach (i, v; vertices) { v.id = to!int(i); }
        //
        // Set references to boundary faces in bc objects.
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                auto f = get_ifj(i, njc, k);
                bc[Face.north].faces ~= f;
                bc[Face.north].outsigns ~= 1;
                f.i_bndry = bc[Face.north].outsigns.length - 1;
            }
        }
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                auto f = get_ifi(nic, j, k);
                bc[Face.east].faces ~= f;
                bc[Face.east].outsigns ~= 1;
                f.i_bndry = bc[Face.east].outsigns.length - 1;
            }
        }
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                auto f = get_ifj(i, 0, k);
                bc[Face.south].faces ~= f;
                bc[Face.south].outsigns ~= -1;
                f.i_bndry = bc[Face.south].outsigns.length - 1;
            }
        }
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                auto f = get_ifi(0, j, k);
                bc[Face.west].faces ~= f;
                bc[Face.west].outsigns ~= -1;
                f.i_bndry = bc[Face.west].outsigns.length - 1;
            }
        }
        if (myConfig.dimensions == 3) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    auto f = get_ifk(i, j, nkc);
                    bc[Face.top].faces ~= f;
                    bc[Face.top].outsigns ~= 1;
                    f.i_bndry = bc[Face.top].outsigns.length - 1;
                }
            }
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    auto f = get_ifk(i, j, 0);
                    bc[Face.bottom].faces ~= f;
                    bc[Face.bottom].outsigns ~= -1;
                    f.i_bndry = bc[Face.bottom].outsigns.length - 1;
                }
            }
        } // end if dimensions == 3
        //
        // Bind interfaces vertices to cells.
        //
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    auto c = get_cell(i,j,k);
                    c.iface.length = (myConfig.dimensions == 3) ? 6 : 4;
                    c.outsign.length = c.iface.length;
                    c.iface[Face.west] = get_ifi(i,j,k); c.outsign[Face.west] = -1;
                    c.iface[Face.east] = get_ifi(i+1,j,k); c.outsign[Face.east] = 1;
                    c.iface[Face.south] = get_ifj(i,j,k); c.outsign[Face.south] = -1;
                    c.iface[Face.north] = get_ifj(i,j+1,k); c.outsign[Face.north] = 1;
                    // VTK order for vertices (on bottom face).
                    c.vtx.length = 0;
                    c.vtx ~= get_vtx(i,j,k);
                    c.vtx ~= get_vtx(i+1,j,k);
                    c.vtx ~= get_vtx(i+1,j+1,k);
                    c.vtx ~= get_vtx(i,j+1,k);
                    if (myConfig.dimensions == 3) {
                        c.iface[Face.bottom] = get_ifk(i,j,k); c.outsign[Face.bottom] = -1;
                        c.iface[Face.top] = get_ifk(i,j,k+1); c.outsign[Face.top] = 1;
                        // VTK order for vertices on top face.
                        c.vtx ~= get_vtx(i,j,k+1);
                        c.vtx ~= get_vtx(i+1,j,k+1);
                        c.vtx ~= get_vtx(i+1,j+1,k+1);
                        c.vtx ~= get_vtx(i,j+1,k+1);
                    }
                }
            }
        }
        // Now do the id based mapping to vertices used for the new datastructures
        size_t nvtxs = 4;
        if (myConfig.dimensions == 3) nvtxs = 8;
        celldata.c2v.length = nkc*njc*nic;
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t c = cell_index(i,j,k);
                    celldata.c2v[c].length = nvtxs;
                    celldata.c2v[c][0] = vertex_index(i,j,k);
                    celldata.c2v[c][1] = vertex_index(i+1,j,k);
                    celldata.c2v[c][2] = vertex_index(i+1,j+1,k);
                    celldata.c2v[c][3] = vertex_index(i,j+1,k);
                    if (myConfig.dimensions == 3) {
                        celldata.c2v[c][4] = vertex_index(i,j,k+1);
                        celldata.c2v[c][5] = vertex_index(i+1,j,k+1);
                        celldata.c2v[c][6] = vertex_index(i+1,j+1,k+1);
                        celldata.c2v[c][7] = vertex_index(i,j+1,k+1);
                    }
                }
            }
        }
        // Test out a index based cell to face mapping
        celldata.c2f.length = nkc*njc*nic;
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t c = cell_index(i,j,k);
                    celldata.c2f[c] ~= ifj_index(i,j+1,k);// north
                    celldata.c2f[c] ~= ifi_index(i+1,j,k);// east
                    celldata.c2f[c] ~= ifj_index(i,j,k);  // south
                    celldata.c2f[c] ~= ifi_index(i,j,k);  // west
                    if (myConfig.dimensions == 3) {
                        celldata.c2f[c] ~= ifk_index(i,j,k+1); // top
                        celldata.c2f[c] ~= ifk_index(i,j,k);   // bottom
                    }
                    celldata.nfaces[c] = celldata.c2f[c].length;
                }
            }
        }
        celldata.outsigns.length = nkc*njc*nic;
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    size_t c = cell_index(i,j,k);
                    celldata.outsigns[c] ~= 1;// north
                    celldata.outsigns[c] ~= 1;// east
                    celldata.outsigns[c] ~= -1;  // south
                    celldata.outsigns[c] ~= -1;  // west
                    if (myConfig.dimensions == 3) {
                        celldata.outsigns[c] ~= 1; // top
                        celldata.outsigns[c] ~= -1;   // bottom
                    }
                }
            }
        }

        //
        // for instances when a numerical Jacobian will be formed (precondition matrix or adjoint operator), it is helpful
        // to have the cell_cloud filled with references to the nearby cells that effect the convective fluxes for the given
        // cell. We may then later treat the structured and unstructured blocks in the same manner when constructing the Jacobian.
        // NB. currently only gathers nearest-neighbours for a (spatially) first-order Jacobian.
        foreach (i, c; cells) {
            int[] cloud_id_list = [];
            c.cell_cloud ~= c;
            cloud_id_list ~= c.id;
            bool cell_found;
            foreach (f; c.iface) {
                if (f.left_cells.length > 0) {
                    cell_found = cloud_id_list.canFind(f.left_cells[0].id);
                    if (!cell_found) {
                        c.cell_cloud ~= f.left_cells[0];
                        cloud_id_list ~= f.left_cells[0].id;
                    }
                }
                if (f.right_cells.length > 0) {
                    cell_found = cloud_id_list.canFind(f.right_cells[0].id);
                    if (!cell_found) {
                        c.cell_cloud ~= f.right_cells[0];
                        cloud_id_list ~= f.right_cells[0].id;
                    }
                }
            }
        }
        // Sometimes it is convenient for an interface to come complete
        // with information about the vertices that define it and also
        // the cells that adjoin it, as for the unstructured grid.
        //
        // ifi interfaces are west interfaces, with their unit normal pointing east.
        // In 2D, vtx0==p00, vtx1==p01.
        // In 3D, the cycle [vtx0,vtx1,vtx2,vtx3] progresses counter-clockwise around
        // the periphery of the face when the normal unit vector is pointing toward you.
        // t1 vector aligned with j-index direction
        // t2 vector aligned with k-index direction
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. niv) {
                    auto f = get_ifi(i,j,k);
                    f.vtx.length = 0;
                    if (myConfig.dimensions == 3) {
                        f.vtx ~= get_vtx(i,j,k);
                        f.vtx ~= get_vtx(i,j+1,k);
                        f.vtx ~= get_vtx(i,j+1,k+1);
                        f.vtx ~= get_vtx(i,j,k+1);
                    } else {
                        f.vtx ~= get_vtx(i,j);
                        f.vtx ~= get_vtx(i,j+1);
                    }
                    if (i == 0) {
                        f.is_on_boundary = true;
                        f.bc_id = Face.west;
                        if (bc[Face.west].ghost_cell_data_available) {
                            f.left_cell = f.left_cells[0];
                        }
                        f.right_cell = get_cell(i,j,k);
                    } else if (i == nic) {
                        f.is_on_boundary = true;
                        f.bc_id = Face.east;
                        f.left_cell = get_cell(i-1,j,k);
                        if (bc[Face.east].ghost_cell_data_available) {
                            f.right_cell = f.right_cells[0];
                        }
                    } else {
                        f.left_cell = get_cell(i-1,j,k);
                        f.right_cell = get_cell(i,j,k);
                    }
                } // i loop
            } // j loop
        } // for k
        // ifj interfaces are south interfaces, with their unit normal pointing north.
        // In 2D, vtx0==p10, vtx1==p00.
        // t1 vector aligned with k-index direction
        // t2 vector aligned with i-index direction
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njv) {
                    auto f = get_ifj(i,j,k);
                    f.vtx.length = 0;
                    if (myConfig.dimensions == 3) {
                        f.vtx ~= get_vtx(i,j,k);
                        f.vtx ~= get_vtx(i,j,k+1);
                        f.vtx ~= get_vtx(i+1,j,k+1);
                        f.vtx ~= get_vtx(i+1,j,k);
                    } else {
                        f.vtx ~= get_vtx(i+1,j);
                        f.vtx ~= get_vtx(i,j);
                    }
                    if (j == 0) {
                        f.is_on_boundary = true;
                        f.bc_id = Face.south;
                        if (bc[Face.south].ghost_cell_data_available) {
                            f.left_cell = f.left_cells[0];
                        }
                        f.right_cell = get_cell(i,j,k);
                    } else if (j == njc) {
                        f.is_on_boundary = true;
                        f.bc_id = Face.north;
                        f.left_cell = get_cell(i,j-1,k);
                        if (bc[Face.north].ghost_cell_data_available) {
                            f.right_cell = f.right_cells[0];
                        }
                    } else {
                        f.left_cell = get_cell(i,j-1,k);
                        f.right_cell = get_cell(i,j,k);
                    }
                } // j loop
            } // i loop
        } // for k
        if (myConfig.dimensions == 3) {
            // ifk interfaces are bottom interfaces, with unit normal pointing to top.
            // t1 vector aligned with i-index direction
            // t2 vector aligned with j-index direction
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njc) {
                    foreach (k; 0 .. nkv) {
                        auto f = get_ifk(i,j,k);
                        f.vtx.length = 0;
                        f.vtx ~= get_vtx(i,j,k);
                        f.vtx ~= get_vtx(i+1,j,k);
                        f.vtx ~= get_vtx(i+1,j+1,k);
                        f.vtx ~= get_vtx(i,j+1,k);
                        if (k == 0) {
                            f.is_on_boundary = true;
                            f.bc_id = Face.bottom;
                            if (bc[Face.bottom].ghost_cell_data_available) {
                                f.left_cell = f.left_cells[0];
                            }
                            f.right_cell = get_cell(i,j,k);
                        } else if (k == nkc) {
                            f.is_on_boundary = true;
                            f.bc_id = Face.top;
                            f.left_cell = get_cell(i,j,k-1);
                            if (bc[Face.top].ghost_cell_data_available) {
                                f.right_cell = f.right_cells[0];
                            }
                        } else {
                            f.left_cell = get_cell(i,j,k-1);
                            f.right_cell = get_cell(i,j,k);
                        }
                    } // for k
                } // j loop
            } // i loop
        }
        //
        store_references_for_derivative_calc(0);
        grid = new StructuredGrid(gridFileName, myConfig.grid_format);
        grid.sort_cells_into_bins();
        sync_vertices_from_underlying_grid(0);
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

    @nogc
    override void compute_primary_cell_geometric_data(size_t gtl)
    // Compute cell and interface geometric properties.
    {
        Vector3 dummy;
        Vector3 ds;
        if (myConfig.dimensions == 2) {
            foreach (c; cells) { c.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
            foreach (f; faces) { f.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
        } else { // 3D
            foreach (c; cells) { c.update_3D_geometric_data(gtl, myConfig.true_centroids); }
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
        //
        /* Extrapolate (with first-order) cell positions and volumes to ghost cells. */
        // TODO -- think about how to make these things consistent.
        @nogc
        void extrap(ref Vector3 pos, ref const(Vector3) p1, ref const(Vector3) p2)
        {
            pos.set(p1); pos *= 2.0; pos -= p2;
        }
        //
        auto option = CopyDataOption.cell_lengths_only;
        bool nghost3 = n_ghost_cell_layers == 3;
        //
        if (bc[Face.east].ghost_cell_data_available) {
            if (n_ghost_cell_layers > nic) {
                throw new FlowSolverException("Too few cells in i-direction for ghost cell copies.");
            }
            foreach (j; 0 .. njc) {
                foreach (k; 0 .. nkc) {
                    auto f = get_ifi(nic, j, k);
                    foreach (n; 0 .. n_ghost_cell_layers) {
                        f.right_cells[n].copy_values_from(get_cell(nic-1-n,j,k), option);
                    }
                    auto cell_1 = get_cell(nic-1,j,k);
                    auto cell_2 = get_cell(nic-2,j,k);
                    auto ghost_cell = f.right_cells[0];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = f.right_cells[1];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    if (nghost3) {
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = f.right_cells[2];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    }
                    foreach (n; 0 .. n_ghost_cell_layers) f.right_cells[n].update_celldata_geometry();
                }
            }
        }
        if (bc[Face.west].ghost_cell_data_available) {
            if (n_ghost_cell_layers > nic) {
                throw new FlowSolverException("Too few cells in i-direction for ghost cell copies.");
            }
            foreach (j; 0 .. njc) {
                foreach (k; 0 .. nkc) {
                    auto f = get_ifi(0, j, k);
                    foreach (n; 0 .. n_ghost_cell_layers) {
                        f.left_cells[n].copy_values_from(get_cell(n,j,k), option);
                    }
                    auto cell_1 = get_cell(0,j,k);
                    auto cell_2 = get_cell(1,j,k);
                    auto ghost_cell = f.left_cells[0];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = f.left_cells[1];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    if (nghost3) {
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = f.left_cells[2];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    }
                    foreach (n; 0 .. n_ghost_cell_layers) f.left_cells[n].update_celldata_geometry();
                }
            }
        }
        if (bc[Face.north].ghost_cell_data_available) {
            if (n_ghost_cell_layers > njc) {
                throw new FlowSolverException("Too few cells in j-direction for ghost cell copies.");
            }
            foreach (i; 0 .. nic) {
                foreach (k; 0 .. nkc) {
                    auto f = get_ifj(i, njc, k);
                    foreach (n; 0 .. n_ghost_cell_layers) {
                        f.right_cells[n].copy_values_from(get_cell(i,njc-1-n,k), option);
                    }
                    auto cell_1 = get_cell(i,njc-1,k);
                    auto cell_2 = get_cell(i,njc-2,k);
                    auto ghost_cell = f.right_cells[0];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = f.right_cells[1];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    if (nghost3) {
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = f.right_cells[2];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    }
                    foreach (n; 0 .. n_ghost_cell_layers) f.right_cells[n].update_celldata_geometry();
                }
            }
        }
        if (bc[Face.south].ghost_cell_data_available) {
            if (n_ghost_cell_layers > njc) {
                throw new FlowSolverException("Too few cells in j-direction for ghost cell copies.");
            }
            foreach (i; 0 .. nic) {
                foreach (k; 0 .. nkc) {
                    auto f = get_ifj(i, 0, k);
                    foreach (n; 0 .. n_ghost_cell_layers) {
                        f.left_cells[n].copy_values_from(get_cell(i,n,k), option);
                    }
                    auto cell_1 = get_cell(i,0,k);
                    auto cell_2 = get_cell(i,1,k);
                    auto ghost_cell = f.left_cells[0];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    cell_2 = cell_1;
                    cell_1 = ghost_cell;
                    ghost_cell = f.left_cells[1];
                    extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                    ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    if (nghost3) {
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = f.left_cells[2];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                    }
                    foreach (n; 0 .. n_ghost_cell_layers) f.left_cells[n].update_celldata_geometry();
                }
            }
        }
        if (myConfig.dimensions == 3) {
            if (bc[Face.top].ghost_cell_data_available) {
                if (n_ghost_cell_layers > nkc) {
                    throw new FlowSolverException("Too few cells in k-direction for ghost cell copies.");
                }
                foreach (i; 0 .. nic) {
                    foreach (j; 0 .. njc) {
                        auto f = get_ifk(i, j, nkc);
                        foreach (n; 0 .. n_ghost_cell_layers) {
                            f.right_cells[n].copy_values_from(get_cell(i,j,nkc-1-n), option);
                        }
                        auto cell_1 = get_cell(i,j,nkc-1);
                        auto cell_2 = get_cell(i,j,nkc-2);
                        auto ghost_cell = f.right_cells[0];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = f.right_cells[1];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                        if (nghost3) {
                            cell_2 = cell_1;
                            cell_1 = ghost_cell;
                            ghost_cell = f.right_cells[2];
                            extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                            ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                        }
                        foreach (n; 0 .. n_ghost_cell_layers) f.right_cells[n].update_celldata_geometry();
                    }
                }
            }
            if (bc[Face.bottom].ghost_cell_data_available) {
                if (n_ghost_cell_layers > nkc) {
                    throw new FlowSolverException("Too few cells in k-direction for ghost cell copies.");
                }
                foreach (i; 0 .. nic) {
                    foreach (j; 0 .. njc) {
                        auto f = get_ifk(i, j, 0);
                        foreach (n; 0 .. n_ghost_cell_layers) {
                            f.left_cells[n].copy_values_from(get_cell(i,j,n), option);
                        }
                        auto cell_1 = get_cell(i,j,0);
                        auto cell_2 = get_cell(i,j,1);
                        auto ghost_cell = f.left_cells[0];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                        cell_2 = cell_1;
                        cell_1 = ghost_cell;
                        ghost_cell = f.left_cells[1];
                        extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                        ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                        if (nghost3) {
                            cell_2 = cell_1;
                            cell_1 = ghost_cell;
                            ghost_cell = f.left_cells[2];
                            extrap(ghost_cell.pos[gtl], cell_1.pos[gtl], cell_2.pos[gtl]);
                            ghost_cell.volume[gtl] = 2.0*cell_1.volume[gtl] - cell_2.volume[gtl];
                        }
                        foreach (n; 0 .. n_ghost_cell_layers) f.left_cells[n].update_celldata_geometry();
                    }
                }
            }
        } // end if dimensions == 3

    } // end compute_primary_cell_geometric_data()

    @nogc
    override void precompute_stencil_data(size_t gtl)
    {
        bool second_order = myConfig.interpolation_order>=2;
        bool third_order  = myConfig.interpolation_order==3;
        // Update stencil information, or set it all for the first time WIP (NNG)
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. niv) {
                    size_t fid = ifi_index(i,j,k);
                    number lenL1 = celldata.lengths[facedata.stencil_idxs[fid].L1][0];
                    number lenL0 = celldata.lengths[facedata.stencil_idxs[fid].L0][0];
                    number lenR0 = celldata.lengths[facedata.stencil_idxs[fid].R0][0];
                    number lenR1 = celldata.lengths[facedata.stencil_idxs[fid].R1][0];
                    if (second_order) {
                        facedata.l2r2_interp_data[fid].set(lenL1, lenL0, lenR0, lenR1);
                    }
                    if (third_order) {
                        number lenL2 = celldata.lengths[facedata.stencil_idxs[fid].L2][0];
                        number lenR2 = celldata.lengths[facedata.stencil_idxs[fid].R2][0];
                        facedata.l3r3_interp_data[fid].set(lenL2, lenL1, lenL0, lenR0, lenR1, lenR2);
                    }
                } // i loop
            } // j loop
        } // k loop
        foreach (k; 0 .. nkc) {
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njv) {
                    size_t fid = ifj_index(i,j,k);
                    number lenL1 = celldata.lengths[facedata.stencil_idxs[fid].L1][1];
                    number lenL0 = celldata.lengths[facedata.stencil_idxs[fid].L0][1];
                    number lenR0 = celldata.lengths[facedata.stencil_idxs[fid].R0][1];
                    number lenR1 = celldata.lengths[facedata.stencil_idxs[fid].R1][1];
                    if (second_order) {
                        facedata.l2r2_interp_data[fid].set(lenL1, lenL0, lenR0, lenR1);
                    }
                    if (third_order) {
                        number lenL2 = celldata.lengths[facedata.stencil_idxs[fid].L2][1];
                        number lenR2 = celldata.lengths[facedata.stencil_idxs[fid].R2][1];
                        facedata.l3r3_interp_data[fid].set(lenL2, lenL1, lenL0, lenR0, lenR1, lenR2);
                    }
                } // j loop
            } // i loop
        } // k loop
        if (myConfig.dimensions == 3) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    foreach (k; 0 .. nkv) {
                        size_t fid = ifk_index(i,j,k);
                        number lenL1 = celldata.lengths[facedata.stencil_idxs[fid].L1][2];
                        number lenL0 = celldata.lengths[facedata.stencil_idxs[fid].L0][2];
                        number lenR0 = celldata.lengths[facedata.stencil_idxs[fid].R0][2];
                        number lenR1 = celldata.lengths[facedata.stencil_idxs[fid].R1][2];
                        if (second_order) {
                            facedata.l2r2_interp_data[fid].set(lenL1, lenL0, lenR0, lenR1);
                        }
                        if (third_order) {
                            number lenL2 = celldata.lengths[facedata.stencil_idxs[fid].L2][2];
                            number lenR2 = celldata.lengths[facedata.stencil_idxs[fid].R2][2];
                            facedata.l3r3_interp_data[fid].set(lenL2, lenL1, lenL0, lenR0, lenR1, lenR2);
                        }
                    } // k loop
                } // i loop
            } // j loop
        } // end 3d

    } // end precompute_stencil_data()

    @nogc
    override void compute_least_squares_setup(size_t gtl)
    {
        // Update the least-squares geometric weights and the workspaces, if appropriate.
        // The weights should be calculated when the grid is initialised or moved.
        // They are needed for flow gradient calculations that feed into the viscous fluxes.
        if (myConfig.viscous && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)) {
            final switch (myConfig.spatial_deriv_locn) {
            case SpatialDerivLocn.vertices:
                foreach(vtx; vertices) {
                    if (myConfig.viscous_least_squares_type == ViscousLeastSquaresType.weighted_normal ||
                        myConfig.viscous_least_squares_type == ViscousLeastSquaresType.unweighted_normal) {
                        vtx.grad.set_up_workspace_leastsq_via_normal(myConfig, vtx.cloud_pos, vtx.pos[gtl], true, *(vtx.ws_grad));
                    } else {
                        vtx.grad.set_up_workspace_leastsq_via_qr_factorization(myConfig, vtx.cloud_pos, vtx.pos[gtl], true, *(vtx.ws_grad));
                    }
                }
                break;
            case SpatialDerivLocn.faces:
                foreach(iface; faces) {
                    if (myConfig.viscous_least_squares_type == ViscousLeastSquaresType.weighted_normal ||
                        myConfig.viscous_least_squares_type == ViscousLeastSquaresType.unweighted_normal) {
                        iface.grad.set_up_workspace_leastsq_via_normal(myConfig, iface.cloud_pos, iface.pos, false, *(iface.ws_grad));
                    } else {
                        iface.grad.set_up_workspace_leastsq_via_qr_factorization(myConfig, iface.cloud_pos, iface.pos, false, *(iface.ws_grad));
                    }
                }
                break;
            case SpatialDerivLocn.cells:
                foreach(cell; cells) {
                    if (myConfig.viscous_least_squares_type == ViscousLeastSquaresType.weighted_normal ||
                        myConfig.viscous_least_squares_type == ViscousLeastSquaresType.unweighted_normal) {
                        cell.grad.set_up_workspace_leastsq_via_normal(myConfig, cell.cloud_pos, cell.pos[gtl], false, *(cell.ws_grad));
                    } else {
                        cell.grad.set_up_workspace_leastsq_via_qr_factorization(myConfig, cell.cloud_pos, cell.pos[gtl], false, *(cell.ws_grad));
                    }
                }
            } // end switch
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
            store_references_for_derivative_calc_at_cells(gtl);
            break;
        }
    } // end store_references_for_derivative_calc()

    void store_references_for_derivative_calc_at_cells(size_t gtl)
    {
        // This locations is only valid for the weighted least squares calculation.
        foreach (c; cells) {
            // First cell in the cloud is the cell itself.  Differences are taken about it.
            c.cloud_pos ~= &(c.pos[0]);
            c.cloud_fs ~= c.fs;
            // Subsequent cells are the surrounding interfaces.
            foreach (i, f; c.iface) {
                c.cloud_pos ~= &(f.pos);
                c.cloud_fs ~= f.fs;
            } // end foreach face
        }
        // Check that we have correctly assembled clouds.
        foreach (i; 0 .. nic) {
            foreach (j; 0 .. njc) {
                foreach (k; 0 .. nkc) {
                    auto c = get_cell(i,j,k);
                    auto ncloud = c.cloud_pos.length;
                    if (ncloud < 3) {
                        string msg = format("Too few points in cloud around cell centre[%d,%d,%d] "~
                                            "ncloud=%d nic=%d njc=%d nkc=%d",
                                            i, j, k, ncloud, nic, njc, nkc);
                        throw new FlowSolverException(msg);
                    }
                }
            }
        }
    } // end store_references_for_derivative_calc_at_cells

    void store_references_for_derivative_calc_at_faces(size_t gtl)
    {
        // The weighted least squares calculation is expecting the interface
        // at which the gradient is being calculated to be stored in position [0].
        // However the divergence calculation is expecting a specific ordering of
        // the cloud points, as such we must look up the spatial_deriv_calc type
        // to decide which cloud to use.
        if (myConfig.dimensions == 2) {
            // First, i-faces
            foreach (i; 0 .. niv) {
                foreach (j; 0 .. njc) {
                    FVInterface face = get_ifi(i,j);
                    // Points nearby.
                    if (i == 0) {
                        // west boundary
                        FVInterface D = get_ifj(i,j);
                        FluidFVCell E = get_cell(i,j);
                        FVInterface F = get_ifj(i,j+1);
                        // Retain locations and references to flow states for later.
                        face.cloud_pos = [&(face.pos), &(D.pos), &(E.pos[gtl]), &(F.pos)];
                        face.cloud_fs = [face.fs, D.fs, E.fs, F.fs];
                    } else if (i == nic) {
                        // east boundary
                        FVInterface A = get_ifj(i-1,j+1);
                        FluidFVCell B = get_cell(i-1,j);
                        FVInterface C = get_ifj(i-1,j);
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
                        FVInterface A = get_ifj(i-1,j+1);
                        FluidFVCell B = get_cell(i-1,j);
                        FVInterface C = get_ifj(i-1,j);
                        FVInterface D = get_ifj(i,j);
                        FluidFVCell E = get_cell(i,j);
                        FVInterface F = get_ifj(i,j+1);
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
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njv) {
                    FVInterface face = get_ifj(i,j);
                    // Points nearby.
                    if (j == 0) {
                        // south boundary
                        FVInterface D = get_ifi(i+1,j);
                        FluidFVCell E = get_cell(i,j);
                        FVInterface F = get_ifi(i,j);
                        // Retain locations and references to flow states for later.
                        face.cloud_pos = [&(face.pos), &(D.pos), &(E.pos[gtl]), &(F.pos)];
                        face.cloud_fs = [face.fs, D.fs, E.fs, F.fs];
                    } else if (j == njc) {
                        // north boundary
                        FVInterface A = get_ifi(i,j-1);
                        FluidFVCell B = get_cell(i,j-1);
                        FVInterface C = get_ifi(i+1,j-1);
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
                        FVInterface A = get_ifi(i,j-1);
                        FluidFVCell B = get_cell(i,j-1);
                        FVInterface C = get_ifi(i+1,j-1);
                        FVInterface D = get_ifi(i+1,j);
                        FluidFVCell E = get_cell(i,j);
                        FVInterface F = get_ifi(i,j);
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
            foreach (i; 0 .. niv) {
                foreach (j; 0 .. njc) {
                    foreach (k; 0 .. nkc) {
                        FVInterface face = get_ifi(i,j,k);
                        // Points nearby.
                        if (i == 0) {
                            // west boundary
                            FVInterface F = get_ifj(i,j+1,k);
                            FVInterface G = get_ifj(i,j,k);
                            FVInterface H = get_ifk(i,j,k+1);
                            FVInterface I = get_ifk(i,j,k);
                            FluidFVCell J = get_cell(i,j,k);
                            // Retain locations and references to flow states for later.
                            face.cloud_pos = [&(face.pos), &(F.pos), &(G.pos), &(H.pos),
                                              &(I.pos), &(J.pos[gtl])];
                            face.cloud_fs = [face.fs, F.fs, G.fs, H.fs, I.fs, J.fs];
                        } else if (i == nic) {
                            // east boundary
                            FVInterface A = get_ifj(i-1,j+1,k);
                            FVInterface B = get_ifj(i-1,j,k);
                            FVInterface C = get_ifk(i-1,j,k+1);
                            FVInterface D = get_ifk(i-1,j,k);
                            FluidFVCell E = get_cell(i-1,j,k);
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
                            FVInterface A = get_ifj(i-1,j+1,k);
                            FVInterface B = get_ifj(i-1,j,k);
                            FVInterface C = get_ifk(i-1,j,k+1);
                            FVInterface D = get_ifk(i-1,j,k);
                            FluidFVCell E = get_cell(i-1,j,k);
                            FVInterface F = get_ifj(i,j+1,k);
                            FVInterface G = get_ifj(i,j,k);
                            FVInterface H = get_ifk(i,j,k+1);
                            FVInterface I = get_ifk(i,j,k);
                            FluidFVCell J = get_cell(i,j,k);
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
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njv) {
                    foreach (k; 0 .. nkc) {
                        FVInterface face = get_ifj(i,j,k);
                        // Points nearby.
                        if (j == 0) {
                            // south boundary
                            FVInterface F = get_ifi(i+1,j,k);
                            FVInterface G = get_ifi(i,j,k);
                            FVInterface H = get_ifk(i,j,k+1);
                            FVInterface I = get_ifk(i,j,k);
                            FluidFVCell J = get_cell(i,j,k);
                            // Retain locations and references to flow states for later.
                            face.cloud_pos = [&(face.pos), &(F.pos), &(G.pos), &(H.pos),
                                              &(I.pos), &(J.pos[gtl])];
                            face.cloud_fs = [face.fs, F.fs, G.fs, H.fs, I.fs, J.fs];
                        } else if (j == njc) {
                            // north boundary
                            FVInterface A = get_ifi(i+1,j-1,k);
                            FVInterface B = get_ifi(i,j-1,k);
                            FVInterface C = get_ifk(i,j-1,k+1);
                            FVInterface D = get_ifk(i,j-1,k);
                            FluidFVCell E = get_cell(i,j-1,k);
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
                            FVInterface A = get_ifi(i+1,j-1,k);
                            FVInterface B = get_ifi(i,j-1,k);
                            FVInterface C = get_ifk(i,j-1,k+1);
                            FVInterface D = get_ifk(i,j-1,k);
                            FluidFVCell E = get_cell(i,j-1,k);
                            FVInterface F = get_ifi(i+1,j,k);
                            FVInterface G = get_ifi(i,j,k);
                            FVInterface H = get_ifk(i,j,k+1);
                            FVInterface I = get_ifk(i,j,k);
                            FluidFVCell J = get_cell(i,j,k);
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
            foreach (i; 0 .. nic) {
                foreach (j; 0 .. njc) {
                    foreach (k; 0 .. nkv) {
                        FVInterface face = get_ifk(i,j,k);
                        // Points nearby.
                        if (k == 0) {
                            // bottom boundary
                            FVInterface F = get_ifj(i,j+1,k);
                            FVInterface G = get_ifj(i,j,k);
                            FVInterface H = get_ifi(i+1,j,k);
                            FVInterface I = get_ifi(i,j,k);
                            FluidFVCell J = get_cell(i,j,k);
                            // Retain locations and references to flow states for later.
                            face.cloud_pos = [&(face.pos), &(F.pos), &(G.pos), &(H.pos),
                                              &(I.pos), &(J.pos[gtl])];
                            face.cloud_fs = [face.fs, F.fs, G.fs, H.fs, I.fs, J.fs];
                        } else if (k == nkc) {
                            // top boundary
                            FVInterface A = get_ifj(i,j+1,k-1);
                            FVInterface B = get_ifj(i,j,k-1);
                            FVInterface C = get_ifi(i+1,j,k-1);
                            FVInterface D = get_ifi(i,j,k-1);
                            FluidFVCell E = get_cell(i,j,k-1);
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
                            FVInterface A = get_ifj(i,j+1,k-1);
                            FVInterface B = get_ifj(i,j,k-1);
                            FVInterface C = get_ifi(i+1,j,k-1);
                            FVInterface D = get_ifi(i,j,k-1);
                            FluidFVCell E = get_cell(i,j,k-1);
                            FVInterface F = get_ifj(i,j+1,k);
                            FVInterface G = get_ifj(i,j,k);
                            FVInterface H = get_ifi(i+1,j,k);
                            FVInterface I = get_ifi(i,j,k);
                            FluidFVCell J = get_cell(i,j,k);
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
        //
        // Check that we have some points in all clouds.
        foreach (f; faces) {
            auto ncloud = f.cloud_pos.length;
            if (ncloud < 3) {
                string msg = format("Too few points in cloud around midpoint of face id=%d ncloud=%d", f.id, ncloud);
                throw new FlowSolverException(msg);
            }
        }
    } // end store_references_for_derivative_calc_at_faces()

    void store_references_for_derivative_calc_at_vertices(size_t gtl)
    {
        if (myConfig.dimensions == 2) {
            // First, do all of the internal secondary cells.
            // i.e. Those not on a boundary.
            foreach (i; 1 .. niv-1) {
                foreach (j; 1 .. njv-1) {
                    // Secondary-cell centre is a primary-cell vertex.
                    FVVertex vtx = get_vtx(i,j);
                    // These are the corners of the secondary cell.
                    FluidFVCell A = get_cell(i,j-1);
                    FluidFVCell B = get_cell(i,j);
                    FluidFVCell C = get_cell(i-1,j);
                    FluidFVCell D = get_cell(i-1,j-1);
                    // Retain locations and references to flow states for later.
                    vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos[gtl]), &(C.pos[gtl]), &(D.pos[gtl])];
                    vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
                } // j loop
            } // i loop
            // Half-cells along the edges of the block.
            // East boundary
            foreach (j; 1 .. njv-1) {
                size_t i = niv-1;
                FVVertex vtx = get_vtx(i,j);
                FVInterface A = get_ifi(i,j-1);
                FVInterface B = get_ifi(i,j);
                FluidFVCell C = get_cell(i-1,j);
                FluidFVCell D = get_cell(i-1,j-1);
                vtx.cloud_pos = [&(A.pos), &(B.pos), &(C.pos[gtl]), &(D.pos[gtl])];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // j loop
            // West boundary
            foreach (j; 1 .. njv-1) {
                size_t i = 0;
                FVVertex vtx = get_vtx(i,j);
                // These are the corners of the secondary cell.
                FluidFVCell A = get_cell(i,j-1);
                FluidFVCell B = get_cell(i,j);
                FVInterface C = get_ifi(i,j);
                FVInterface D = get_ifi(i,j-1);
                vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos[gtl]), &(C.pos), &(D.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // j loop
            // North boundary
            foreach (i; 1 .. niv-1) {
                size_t j = njv-1;
                FVVertex vtx = get_vtx(i,j);
                FluidFVCell A = get_cell(i,j-1);
                FVInterface B = get_ifj(i,j);
                FVInterface C = get_ifj(i-1,j);
                FluidFVCell D = get_cell(i-1,j-1);
                vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos), &(C.pos), &(D.pos[gtl])];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // i loop
            // South boundary
            foreach (i; 1 .. niv-1) {
                size_t j = 0;
                FVVertex vtx = get_vtx(i,j);
                FVInterface A = get_ifj(i,j);
                FluidFVCell B = get_cell(i,j);
                FluidFVCell C = get_cell(i-1,j);
                FVInterface D = get_ifj(i-1,j);
                vtx.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos[gtl]), &(D.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs, D.fs];
            } // i loop
            // For the corners, we are going to use the same divergence-theorem-based
            // gradient calculator and let one edge collapse to a point, thus giving
            // it a triangle to compute over.  This should be fine.
            // North-east corner
            {
                size_t i = niv-1; size_t j = njv-1;
                FVVertex vtx = get_vtx(i,j);
                FVInterface A = get_ifi(i,j-1);
                FVInterface B = get_ifj(i-1,j);
                FluidFVCell C = get_cell(i-1,j-1);
                vtx.cloud_pos = [&(A.pos), &(B.pos), &(C.pos[gtl])];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
            // South-east corner
            {
                size_t i = niv-1; size_t j = 0;
                FVVertex vtx = get_vtx(i,j);
                FVInterface A = get_ifi(i,j);
                FluidFVCell B = get_cell(i-1,j);
                FVInterface C = get_ifj(i-1,j);
                vtx.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
            // South-west corner
            {
                size_t i = 0; size_t j = 0;
                FVVertex vtx = get_vtx(i,j);
                FVInterface A = get_ifj(i,j);
                FluidFVCell B = get_cell(i,j);
                FVInterface C = get_ifi(i,j);
                vtx.cloud_pos = [&(A.pos), &(B.pos[gtl]), &(C.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
            // North-west corner
            {
                size_t i = 0; size_t j = njv-1;
                FVVertex vtx = get_vtx(i,j);
                FluidFVCell A = get_cell(i,j-1);
                FVInterface B = get_ifj(i,j);
                FVInterface C = get_ifi(i,j-1);
                vtx.cloud_pos = [&(A.pos[gtl]), &(B.pos), &(C.pos)];
                vtx.cloud_fs = [A.fs, B.fs, C.fs];
            }
        } else { // Flow quantity derivatives for 3D.
            // Internal secondary cell geometry information
            foreach (i; 0 .. niv-2) {
                foreach (j; 0 .. njv-2) {
                    foreach (k; 0 .. nkv-2) {
                        FVVertex vtx = get_vtx(i+1,j+1,k+1);
                        FluidFVCell c0 = get_cell(i,j,k);
                        FluidFVCell c1 = get_cell(i+1,j,k);
                        FluidFVCell c2 = get_cell(i+1,j+1,k);
                        FluidFVCell c3 = get_cell(i,j+1,k);
                        FluidFVCell c4 = get_cell(i,j,k+1);
                        FluidFVCell c5 = get_cell(i+1,j,k+1);
                        FluidFVCell c6 = get_cell(i+1,j+1,k+1);
                        FluidFVCell c7 = get_cell(i,j+1,k+1);
                        vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos[gtl]), &(c3.pos[gtl]),
                                         &(c4.pos[gtl]), &(c5.pos[gtl]), &(c6.pos[gtl]), &(c7.pos[gtl])];
                        vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                    }
                }
            }
            // East boundary secondary cell geometry information
            foreach (j; 0 .. njv-2) {
                foreach (k; 0 .. nkv-2) {
                    size_t i = niv-2;
                    FVVertex vtx = get_vtx(i+1,j+1,k+1);
                    FluidFVCell c0 = get_cell(i,j,k);
                    FVInterface c1 = get_ifi(i+1,j,k);
                    FVInterface c2 = get_ifi(i+1,j+1,k);
                    FluidFVCell c3 = get_cell(i,j+1,k);
                    FluidFVCell c4 = get_cell(i,j,k+1);
                    FVInterface c5 = get_ifi(i+1,j,k+1);
                    FVInterface c6 = get_ifi(i+1,j+1,k+1);
                    FluidFVCell c7 = get_cell(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos[gtl]),
                                     &(c4.pos[gtl]), &(c5.pos), &(c6.pos), &(c7.pos[gtl])];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // West boundary secondary cell geometry information
            foreach (j; 0 .. njv-2) {
                foreach (k; 0 .. nkv-2) {
                    int i = -1; // note negative value
                    FVVertex vtx = get_vtx(i+1,j+1,k+1);
                    FVInterface c0 = get_ifi(i+1,j,k);
                    FluidFVCell c1 = get_cell(i+1,j,k);
                    FluidFVCell c2 = get_cell(i+1,j+1,k);
                    FVInterface c3 = get_ifi(i+1,j+1,k);
                    FVInterface c4 = get_ifi(i+1,j,k+1);
                    FluidFVCell c5 = get_cell(i+1,j,k+1);
                    FluidFVCell c6 = get_cell(i+1,j+1,k+1);
                    FVInterface c7 = get_ifi(i+1,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos), &(c1.pos[gtl]), &(c2.pos[gtl]), &(c3.pos),
                                     &(c4.pos), &(c5.pos[gtl]), &(c6.pos[gtl]), &(c7.pos)];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // North boundary secondary cell geometry information
            foreach (i; 0 .. niv-2) {
                foreach (k; 0 .. nkv-2) {
                    size_t j = njv-2;
                    FVVertex vtx = get_vtx(i+1,j+1,k+1);
                    FluidFVCell c0 = get_cell(i,j,k);
                    FluidFVCell c1 = get_cell(i+1,j,k);
                    FVInterface c2 = get_ifj(i+1,j+1,k);
                    FVInterface c3 = get_ifj(i,j+1,k);
                    FluidFVCell c4 = get_cell(i,j,k+1);
                    FluidFVCell c5 = get_cell(i+1,j,k+1);
                    FVInterface c6 = get_ifj(i+1,j+1,k+1);
                    FVInterface c7 = get_ifj(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos), &(c3.pos),
                                     &(c4.pos[gtl]), &(c5.pos[gtl]), &(c6.pos), &(c7.pos)];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // South boundary secondary cell geometry information
            foreach (i; 0 .. niv-2) {
                foreach (k; 0 .. nkv-2) {
                    int j = -1; // note negative
                    FVVertex vtx = get_vtx(i+1,j+1,k+1);
                    FVInterface c0 = get_ifj(i,j+1,k);
                    FVInterface c1 = get_ifj(i+1,j+1,k);
                    FluidFVCell c2 = get_cell(i+1,j+1,k);
                    FluidFVCell c3 = get_cell(i,j+1,k);
                    FVInterface c4 = get_ifj(i,j+1,k+1);
                    FVInterface c5 = get_ifj(i+1,j+1,k+1);
                    FluidFVCell c6 = get_cell(i+1,j+1,k+1);
                    FluidFVCell c7 = get_cell(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos), &(c1.pos), &(c2.pos[gtl]), &(c3.pos[gtl]),
                                     &(c4.pos), &(c5.pos), &(c6.pos[gtl]), &(c7.pos[gtl])];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // Top boundary secondary cell geometry information
            foreach (i; 0 .. niv-2) {
                foreach (j; 0 .. njv-2) {
                    size_t k = nkv-2;
                    FVVertex vtx = get_vtx(i+1,j+1,k+1);
                    FluidFVCell c0 = get_cell(i,j,k);
                    FluidFVCell c1 = get_cell(i+1,j,k);
                    FluidFVCell c2 = get_cell(i+1,j+1,k);
                    FluidFVCell c3 = get_cell(i,j+1,k);
                    FVInterface c4 = get_ifk(i,j,k+1);
                    FVInterface c5 = get_ifk(i+1,j,k+1);
                    FVInterface c6 = get_ifk(i+1,j+1,k+1);
                    FVInterface c7 = get_ifk(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos[gtl]), &(c3.pos[gtl]),
                                     &(c4.pos), &(c5.pos), &(c6.pos), &(c7.pos)];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // Bottom boundary secondary cell geometry information
            foreach (i; 0 .. niv-2) {
                foreach (j; 0 .. njv-2) {
                    int k = -1; // note negative
                    FVVertex vtx = get_vtx(i+1,j+1,k+1);
                    FVInterface c0 = get_ifk(i,j,k+1);
                    FVInterface c1 = get_ifk(i+1,j,k+1);
                    FVInterface c2 = get_ifk(i+1,j+1,k+1);
                    FVInterface c3 = get_ifk(i,j+1,k+1);
                    FluidFVCell c4 = get_cell(i,j,k+1);
                    FluidFVCell c5 = get_cell(i+1,j,k+1);
                    FluidFVCell c6 = get_cell(i+1,j+1,k+1);
                    FluidFVCell c7 = get_cell(i,j+1,k+1);
                    vtx.cloud_pos = [&(c0.pos), &(c1.pos), &(c2.pos), &(c3.pos),
                                     &(c4.pos[gtl]), &(c5.pos[gtl]), &(c6.pos[gtl]), &(c7.pos[gtl])];
                    vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs, c6.fs, c7.fs];
                }
            }
            // Now, do the 4 edges around the bottom face.
            // Bottom-South edge [0]-->[1]
            foreach (i; 1 .. niv-1) {
                size_t j = 0; size_t k = 0;
                FVVertex vtx = get_vtx(i,j,k);
                FluidFVCell c0 = get_cell(i-1,j,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifj(i-1,j,k);
                FVInterface c3 = get_ifk(i-1,j,k);
                FVInterface c4 = get_ifj(i,j,k);
                FVInterface c5 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Bottom-North edge [3]-->[2]
            foreach (i; 1 .. niv-1) {
                size_t j = njv-2; size_t k = 0;
                FVVertex vtx = get_vtx(i,j+1,k);
                FluidFVCell c0 = get_cell(i-1,j,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifj(i-1,j+1,k);
                FVInterface c3 = get_ifk(i-1,j,k);
                FVInterface c4 = get_ifj(i,j+1,k);
                FVInterface c5 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Bottom-West edge [0]-->[3]
            foreach (j; 1 .. njv-1) {
                size_t i = 0; size_t k = 0;
                FVVertex vtx = get_vtx(i,j,k);
                FluidFVCell c0 = get_cell(i,j-1,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i,j-1,k);
                FVInterface c3 = get_ifk(i,j-1,k);
                FVInterface c4 = get_ifi(i,j,k);
                FVInterface c5 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Bottom-East edge [1]-->[2]
            foreach (j; 1 .. njv-1) {
                size_t i = niv-2; size_t k = 0;
                FVVertex vtx = get_vtx(i+1,j,k);
                FluidFVCell c0 = get_cell(i,j-1,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i+1,j-1,k);
                FVInterface c3 = get_ifk(i,j-1,k);
                FVInterface c4 = get_ifi(i+1,j,k);
                FVInterface c5 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // 4 edges around the top face.
            // Top-South edge [4]-->[5]
            foreach (i; 1 .. niv-1) {
                size_t j = 0; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i,j,k+1);
                FluidFVCell c0 = get_cell(i-1,j,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifj(i-1,j,k);
                FVInterface c3 = get_ifk(i-1,j,k+1);
                FVInterface c4 = get_ifj(i,j,k);
                FVInterface c5 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Top-North edge [7]-->[6]
            foreach (i; 1 .. niv-1) {
                size_t j = njv-2; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i,j+1,k+1);
                FluidFVCell c0 = get_cell(i-1,j,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifj(i-1,j+1,k);
                FVInterface c3 = get_ifk(i-1,j,k+1);
                FVInterface c4 = get_ifj(i,j+1,k);
                FVInterface c5 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Top-West edge [4]-->[7]
            foreach (j; 1 .. njv-1) {
                size_t i = 0; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i,j,k+1);
                FluidFVCell c0 = get_cell(i,j-1,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i,j-1,k);
                FVInterface c3 = get_ifk(i,j-1,k+1);
                FVInterface c4 = get_ifi(i,j,k);
                FVInterface c5 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Top-East edge [5]-->[6]
            foreach (j; 1 .. njv-1) {
                size_t i = niv-2; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i+1,j,k+1);
                FluidFVCell c0 = get_cell(i,j-1,k);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i+1,j-1,k);
                FVInterface c3 = get_ifk(i,j-1,k+1);
                FVInterface c4 = get_ifi(i+1,j,k);
                FVInterface c5 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // 4 edges running from bottom to top.
            // South-West edge [0]-->[4]
            foreach (k; 1 .. nkv-1) {
                size_t i = 0; size_t j = 0;
                FVVertex vtx = get_vtx(i,j,k);
                FluidFVCell c0 = get_cell(i,j,k-1);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i,j,k-1);
                FVInterface c3 = get_ifj(i,j,k-1);
                FVInterface c4 = get_ifi(i,j,k);
                FVInterface c5 = get_ifj(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // South-East edge [1]-->[5]
            foreach (k; 1 .. nkv-1) {
                size_t i = niv-2; size_t j = 0;
                FVVertex vtx = get_vtx(i+1,j,k);
                FluidFVCell c0 = get_cell(i,j,k-1);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i+1,j,k-1);
                FVInterface c3 = get_ifj(i,j,k-1);
                FVInterface c4 = get_ifi(i+1,j,k);
                FVInterface c5 = get_ifj(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // North-East edge [2]-->[6]
            foreach (k; 1 .. nkv-1) {
                size_t i = niv-2; size_t j = njv-2;
                FVVertex vtx = get_vtx(i+1,j+1,k);
                FluidFVCell c0 = get_cell(i,j,k-1);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i+1,j,k-1);
                FVInterface c3 = get_ifj(i,j+1,k-1);
                FVInterface c4 = get_ifi(i+1,j,k);
                FVInterface c5 = get_ifj(i,j+1,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // North-West edge [3]-->[7]
            foreach (k; 1 .. nkv-1) {
                size_t i = 0; size_t j = njv-2;
                FVVertex vtx = get_vtx(i,j+1,k);
                FluidFVCell c0 = get_cell(i,j,k-1);
                FluidFVCell c1 = get_cell(i,j,k);
                FVInterface c2 = get_ifi(i,j,k-1);
                FVInterface c3 = get_ifj(i,j+1,k-1);
                FVInterface c4 = get_ifi(i,j,k);
                FVInterface c5 = get_ifj(i,j+1,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos[gtl]), &(c2.pos),
                                 &(c3.pos), &(c4.pos), &(c5.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs, c4.fs, c5.fs];
            }
            // Finally, the 8 corners.
            // South-West-Bottom corner [0]
            {
                size_t i = 0; size_t j = 0; size_t k = 0;
                FVVertex vtx = get_vtx(i,j,k);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i,j,k);
                FVInterface c2 = get_ifj(i,j,k);
                FVInterface c3 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // South-East-Bottom corner [1]
            {
                size_t i = niv-2; size_t j = 0; size_t k = 0;
                FVVertex vtx = get_vtx(i+1,j,k);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i+1,j,k);
                FVInterface c2 = get_ifj(i,j,k);
                FVInterface c3 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-East-Bottom corner [2]
            {
                size_t i = niv-2; size_t j = njv-2; size_t k = 0;
                FVVertex vtx = get_vtx(i+1,j+1,k);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i+1,j,k);
                FVInterface c2 = get_ifj(i,j+1,k);
                FVInterface c3 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-West-Bottom corner [3]
            {
                size_t i = 0; size_t j = njv-2; size_t k = 0;
                FVVertex vtx = get_vtx(i,j+1,k);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i,j,k);
                FVInterface c2 = get_ifj(i,j+1,k);
                FVInterface c3 = get_ifk(i,j,k);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // South-West-Top corner [4]
            {
                size_t i = 0; size_t j = 0; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i,j,k+1);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i,j,k);
                FVInterface c2 = get_ifj(i,j,k);
                FVInterface c3 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // South-East-Top corner [5]
            {
                size_t i = niv-2; size_t j = 0; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i+1,j,k+1);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i+1,j,k);
                FVInterface c2 = get_ifj(i,j,k);
                FVInterface c3 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-East-Top corner [6]
            {
                size_t i = niv-2; size_t j = njv-2; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i+1,j+1,k+1);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i+1,j,k);
                FVInterface c2 = get_ifj(i,j+1,k);
                FVInterface c3 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
            // North-West-Top corner [7]
            {
                size_t i = 0; size_t j = njv-2; size_t k = nkv-2;
                FVVertex vtx = get_vtx(i,j+1,k+1);
                FluidFVCell c0 = get_cell(i,j,k);
                FVInterface c1 = get_ifi(i,j,k);
                FVInterface c2 = get_ifj(i,j+1,k);
                FVInterface c3 = get_ifk(i,j,k+1);
                vtx.cloud_pos = [&(c0.pos[gtl]), &(c1.pos), &(c2.pos), &(c3.pos)];
                vtx.cloud_fs = [c0.fs, c1.fs, c2.fs, c3.fs];
            }
        } // end if (myConfig.dimensions
        //
        // Check that we have correctly assembled clouds.
        foreach (i; 0 .. niv) {
            foreach (j; 0 .. njv) {
                foreach (k; 0 .. nkv) {
                    auto vtx = get_vtx(i,j,k);
                    auto ncloud = vtx.cloud_pos.length;
                    if (ncloud < 3) {
                        string msg = format("Too few points in cloud around vertex[%d,%d,%d] "~
                                            "ncloud=%d niv=%d njv=%d nkv=%d",
                                            i, j, k, ncloud, niv, njv, nkv);
                        throw new FlowSolverException(msg);
                    }
                }
            }
        }
    } // end store_references_for_derivative_calc_at_vertices()

    @nogc
    override void sync_vertices_from_underlying_grid(size_t gtl=0)
    {
        if (myConfig.dimensions == 3) {
            if (grid.niv != niv || grid.njv != njv || grid.nkv != nkv ) {
                string msg = "Mismatch in 3D grid size";
                debug {
                    msg ~= text("\nFor block[", id, "] we have read grid.niv=", grid.niv,
                                " njv=", grid.njv, " nkv=", grid.nkv);
                }
                throw new FlowSolverException(msg);
            }
            foreach (k; 0 .. nkv) {
                foreach (j; 0 .. njv) {
                    foreach (i; 0 .. niv) {
                        auto vtx = get_vtx(i,j,k);
                        auto src_vtx = grid[i,j,k];
                        vtx.pos[gtl].set(src_vtx);
                    }
                }
            }
        } else { // 2D case
            if (grid.niv != niv || grid.njv != njv || grid.nkv != 1) {
                string msg = "Mismatch in 2D grid size";
                debug {
                    msg ~= text("\nFor block[", id, "] we have read grid.niv=", grid.niv,
                                " njv=", grid.njv, " nkv=", grid.nkv);
                }
                throw new FlowSolverException(msg);
            }
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    auto vtx = get_vtx(i,j);
                    auto src_vtx = grid[i,j];
                    vtx.pos[gtl].set(src_vtx.x, src_vtx.y, to!number(0.0));
                }
            }
        }
        return;
    } // end sync_vertices_from_underlying_grid()

    @nogc
    override void sync_vertices_to_underlying_grid(size_t gtl=0)
    // Note that we reuse the StructuredGrid object that was created on the
    // use of init_grid_and_flow_arrays().
    {
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    auto vtx = get_vtx(i,j,k);
                    auto dest_vtx = grid[i,j,k];
                    dest_vtx.set(vtx.pos[gtl]);
                }
            }
        }
        return;
    } // end sync_vertices_to_underlying_grid()

    @nogc
    override void average_turbulent_transprops_to_faces()
    {
        if (!myConfig.turb_model.isTurbulent) return;

        // ifi interior faces
        size_t idx;
        foreach(k; 0 .. nkc) {
            foreach(j; 0 .. njc){
                foreach(i; 1 .. niv-1){
                    idx = ifi_index(i,j,k);
                    faces[idx].fs.mu_t = 0.5*(faces[idx].left_cell.fs.mu_t + faces[idx].right_cell.fs.mu_t);
                    faces[idx].fs.k_t = 0.5*(faces[idx].left_cell.fs.k_t + faces[idx].right_cell.fs.k_t);
                }
            }
        }
        // ifj interior faces
        foreach(k; 0 .. nkc) {
            foreach(j; 1 .. njv-1){
                foreach(i; 0 .. nic){
                    idx = ifj_index(i,j,k);
                    faces[idx].fs.mu_t = 0.5*(faces[idx].left_cell.fs.mu_t + faces[idx].right_cell.fs.mu_t);
                    faces[idx].fs.k_t = 0.5*(faces[idx].left_cell.fs.k_t + faces[idx].right_cell.fs.k_t);
                }
            }
        }

        // ifk interior faces
        foreach(k; 1 .. nkv-1) {
            foreach(j; 0 .. njc){
                foreach(i; 0 .. nic){
                    idx = ifk_index(i,j,k);
                    faces[idx].fs.mu_t = 0.5*(faces[idx].left_cell.fs.mu_t + faces[idx].right_cell.fs.mu_t);
                    faces[idx].fs.k_t = 0.5*(faces[idx].left_cell.fs.k_t + faces[idx].right_cell.fs.k_t);
                }
            }
        }

        // ifi edge faces
        foreach(k; 0 .. nkc){
            foreach(j; 0 .. njc){
                size_t i = 0;
                idx = ifi_index(i,j,k);
                faces[idx].average_turbulent_transprops();

                i = niv-1;
                idx = ifi_index(i,j,k);
                faces[idx].average_turbulent_transprops();
            }
        }

        // ifj edge faces
        foreach(k; 0 .. nkc){
            size_t j = 0;
            foreach(i; 0 .. nic){
                idx = ifj_index(i,j,k);
                faces[idx].average_turbulent_transprops();
            }
            j = njv-1;
            foreach(i; 0 .. nic){
                idx = ifj_index(i,j,k);
                faces[idx].average_turbulent_transprops();
            }
        }

        // Note that nkc is actually 1 in 2D, so we use nkv>1 to test for skipping
        if (nkv>1) {
            // ifk edge faces
            size_t k = 0;
            foreach(j; 0 .. njc){
                foreach(i; 0 .. nic){
                    idx = ifk_index(i,j,k);
                    faces[idx].average_turbulent_transprops();
                }
            }
            k = nkv-1;
            foreach(j; 0 .. njc){
                foreach(i; 0 .. nic){
                    idx = ifk_index(i,j,k);
                    faces[idx].average_turbulent_transprops();
                }
            }
        }
    }

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

    @nogc
    override void propagate_inflow_data_west_to_east()
    {
        // Assume that the west-face ghost cells have appropriate data.
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                auto src_cell = get_ifi(0,j,k).left_cells[0];
                foreach (i; 0 .. nic) {
                    auto dest_cell = get_cell(i,j,k);
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
    override void set_face_flowstates_to_averages_from_cells()
    {
        // It turns out that some shock-detectors need flow derivatives before the
        // convective-flux calculation is done.  That used to be the only place
        // that the face FlowState was filled in and it was done as a side-effect,
        // which has confused just about everyone at some time in their work on the code.
        foreach (f; faces) {
            FluidFVCell cL0 = (f.left_cells.length > 0) ? f.left_cells[0] : f.right_cells[0];
            FluidFVCell cR0 = (f.right_cells.length > 0) ? f.right_cells[0]: f.left_cells[0];
            f.fs.copy_average_values_from(*(cL0.fs), *(cR0.fs));
        }
    }

    @nogc
    void first_order_flux_calc(size_t gtl, size_t[] face_idxs)
    {
        immutable size_t neq = myConfig.cqi.n;
        Vector3 gvel;
        gvel.clear();
        foreach(idx; face_idxs){
            size_t l = facedata.stencil_idxs[idx].L0;
            size_t r = facedata.stencil_idxs[idx].R0;

            Lft.copy_values_from(celldata.flowstates[l]);
            Rght.copy_values_from(celldata.flowstates[r]);

            facedata.flowstates[idx].copy_average_values_from(*Lft, *Rght);

            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
    }

    @nogc
    void second_order_flux_calc(size_t gtl, size_t[] face_idxs)
    {
    /*
        Eilmer's classic second-order piece-wise parabolic reconstruction, using a
        4 cell symmetric stencil.
        Notes: This routine currently ignores the supress_reconstruction face parameter.
    */
        immutable bool hpl = myConfig.apply_heuristic_pressure_based_limiting;
        immutable size_t neq = myConfig.cqi.n;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool is3D = (myConfig.dimensions == 3);
        immutable bool MHD = myConfig.MHD;
        immutable bool apply_limiter = myConfig.apply_limiter;
        immutable bool extrema_clipping = myConfig.extrema_clipping;
        immutable InterpolateOption ti = myConfig.thermo_interpolator;

        number beta = 1.0;
        Vector3 gvel;
        gvel.clear();

        foreach(idx; face_idxs){
            size_t L1 = facedata.stencil_idxs[idx].L1;
            size_t L0 = facedata.stencil_idxs[idx].L0;
            size_t R0 = facedata.stencil_idxs[idx].R0;
            size_t R1 = facedata.stencil_idxs[idx].R1;
            if (hpl) beta =
                compute_heuristic_pressure_limiter(celldata.flowstates[L1].gas.p,
                                                   celldata.flowstates[L0].gas.p,
                                                   celldata.flowstates[R0].gas.p,
                                                   celldata.flowstates[R1].gas.p);

            Lft.copy_values_from(celldata.flowstates[L0]);
            Rght.copy_values_from(celldata.flowstates[R0]);

            interp_l2r2(celldata.flowstates[L1], celldata.flowstates[L0],
                        celldata.flowstates[R0], celldata.flowstates[R1],
                        facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                        facedata.l2r2_interp_data[idx], nsp, nmodes, nturb,
                        ti, MHD, apply_limiter, extrema_clipping,
                        myConfig, *Lft, *Rght, beta);
            facedata.flowstates[idx].copy_average_values_from(*Lft, *Rght);

            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
        return;
    }

    @nogc
    void third_order_flux_calc(size_t gtl, size_t[] face_idxs)
    {
    /*
        Experimental high-order parabolic reconstruction, using a 6 cell
        symmetric stencil, and possibly the unusual ASF flux calculator.
    */
        immutable bool hpl = myConfig.apply_heuristic_pressure_based_limiting;
        immutable size_t neq = myConfig.cqi.n;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool is3D = (myConfig.dimensions == 3);
        immutable bool MHD = myConfig.MHD;
        immutable bool apply_limiter = myConfig.apply_limiter;
        immutable bool extrema_clipping = myConfig.extrema_clipping;
        immutable InterpolateOption ti = myConfig.thermo_interpolator;

        number beta = 1.0;
        Vector3 gvel;
        gvel.clear();

        foreach(idx; face_idxs){
            size_t L2 = facedata.stencil_idxs[idx].L2;
            size_t L1 = facedata.stencil_idxs[idx].L1;
            size_t L0 = facedata.stencil_idxs[idx].L0;
            size_t R0 = facedata.stencil_idxs[idx].R0;
            size_t R1 = facedata.stencil_idxs[idx].R1;
            size_t R2 = facedata.stencil_idxs[idx].R2;
            if (hpl) beta =
                compute_heuristic_pressure_limiter(celldata.flowstates[L2].gas.p,
                                                   celldata.flowstates[L1].gas.p,
                                                   celldata.flowstates[L0].gas.p,
                                                   celldata.flowstates[R0].gas.p,
                                                   celldata.flowstates[R1].gas.p,
                                                   celldata.flowstates[R2].gas.p);

            Lft.copy_values_from(celldata.flowstates[L0]);
            Rght.copy_values_from(celldata.flowstates[R0]);

            interp_l3r3(celldata.flowstates[L2], celldata.flowstates[L1],
                        celldata.flowstates[L0], celldata.flowstates[R0],
                        celldata.flowstates[R1], celldata.flowstates[R2],
                        facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                        facedata.l3r3_interp_data[idx], nsp, nmodes, nturb,
                        ti, MHD, apply_limiter, extrema_clipping,
                        myConfig, *Lft, *Rght, beta);

            facedata.flowstates[idx].copy_average_values_from(*Lft, *Rght);

            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
        return;
    }

    @nogc
    void asf_flux_calc(size_t gtl, size_t[] face_idxs, bool adaptive=false)
    {
    /*
        Low dissipation flux calculator, using the unusual ASF method.
    */
        immutable size_t neq = myConfig.cqi.n;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable bool is3D = (myConfig.dimensions == 3);

        number factor = 1.0;
        foreach(idx; face_idxs){
            size_t L1 = facedata.stencil_idxs[idx].L1;
            size_t L0 = facedata.stencil_idxs[idx].L0;
            size_t R0 = facedata.stencil_idxs[idx].R0;
            size_t R1 = facedata.stencil_idxs[idx].R1;

            facedata.flowstates[idx].copy_average_values_from(celldata.flowstates[L0], celldata.flowstates[R0]);

            if (adaptive) factor = 1.0 - facedata.flowstates[idx].S;

            // This routine uses and reuses a "flux" vector that belongs to the SFluidBlock
            // This is because we need to transform it to the global face, without messing
            // with the components of facedata.fluxes
            flux.clear();
            ASF_242(celldata.flowstates[L1], celldata.flowstates[L0],
                    celldata.flowstates[R0], celldata.flowstates[R1],
                    myConfig, nsp, nmodes, nturb, 
                    facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                    flux, factor);
            foreach(i; 0 .. neq) {
                facedata.fluxes[idx*neq +i] += flux[i];
            }
        }
        return;
    }

    @nogc
    override void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0)
    // Compute the flux from flow-field data on either-side of the interface.
    {
        bool second_order = allow_high_order_interpolation && (myConfig.interpolation_order == 2);
        bool third_order = allow_high_order_interpolation && (myConfig.interpolation_order == 3);
        bool pure_asf = myConfig.flux_calculator == FluxCalculator.asf;
        bool adaptive_asf = myConfig.flux_calculator == FluxCalculator.adaptive_ausmdv_asf;

        //if (face_idxs.length==0) face_idxs = facedata.all_face_idxs;

        if (!pure_asf) {
            if (second_order) {
                second_order_flux_calc(0, facedata.all_face_idxs);
            } else if (third_order) {
                third_order_flux_calc(0, facedata.all_face_idxs);
            } else {
                first_order_flux_calc(0, facedata.all_face_idxs);
            }
        }

        if (pure_asf || adaptive_asf) asf_flux_calc(0, facedata.all_face_idxs, adaptive_asf);
        return;
    } // end convective_flux_phase0()

    @nogc
    override void convective_flux_phase0_legacy(bool allow_high_order_interpolation, size_t gtl=0,
                                         FluidFVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
    // Compute the flux from flow-field data on either-side of the interface.
    {
        // Barring exceptions at the block boundaries, the general process is:
        // (1) interpolate LEFT and RIGHT interface states from cell-center properties.
        // (2) save u, v, w, T for the viscous flux calculation by making a local average.
        // The values for u, v and T may be updated subsequently by the interface-flux function.
        // (3) Apply the flux calculator to the Lft,Rght flow states.
        //
        bool do_reconstruction = allow_high_order_interpolation && (myConfig.interpolation_order > 1);
        if (iface_list.length == 0) { iface_list = faces; }
        //
        // Low-order reconstruction just copies data from adjacent FV_Cell.
        // Note that ,even for high-order reconstruction, we depend upon this copy for
        // the viscous-transport and diffusion coefficients.
        //
        foreach (f; iface_list) {
            if (myConfig.high_order_flux_calculator && f.is_on_boundary && !bc[f.bc_id].ghost_cell_data_available) {
                throw new Error("ghost cell data missing");
            }
            // Typical code path, with interpolation for the flowstates to the left and right of the interface.
            if (do_reconstruction && !f.in_suppress_reconstruction_zone &&
                !(myConfig.suppress_reconstruction_at_shocks && (f.fs.S >= (1.0 - ESSENTIALLY_ZERO)))) {
                one_d.interp(f, *Lft, *Rght);
            } else {
                FluidFVCell cL0 = (f.left_cells.length > 0) ? f.left_cells[0] : f.right_cells[0];
                FluidFVCell cR0 = (f.right_cells.length > 0) ? f.right_cells[0]: f.left_cells[0];
                Lft.copy_values_from(cL0.fs);
                Rght.copy_values_from(cR0.fs);
            }
            f.fs.copy_average_values_from(*Lft, *Rght);
            if (f.is_on_boundary && bc[f.bc_id].convective_flux_computed_in_bc) continue;
            compute_interface_flux(*Lft, *Rght, f, myConfig, omegaz);
        }
        return;
    } // end convective_flux_phase0()

    @nogc
    override void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0,
                                         FluidFVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
    // Compute the flux from data on either-side of the interface.
    // For the structured-grid block, there is nothing to do.
    // The unstructured-grid block needs to work in two phases.
    {
        return;
    }

    @nogc
    override void convective_flux_phase2(bool allow_high_order_interpolation, size_t gtl=0,
                                         FluidFVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
    // Compute the flux from data on either-side of the interface.
    // For the structured-grid block, there is nothing to do.
    // The unstructured-grid block needs to work in two phases.
    {
        return;
    }

    @nogc void copy_current_corner_coords()
    {
        if (myConfig.dimensions == 2) {
            FVVertex vtx00 = get_vtx(0,0);
            corner_coords[0] = vtx00.pos[0].x.re;
            corner_coords[1] = vtx00.pos[0].y.re;
            corner_coords[2] = vtx00.pos[0].z.re;
            FVVertex vtx10 = get_vtx(nic,0);
            corner_coords[3] = vtx10.pos[0].x.re;
            corner_coords[4] = vtx10.pos[0].y.re;
            corner_coords[5] = vtx10.pos[0].z.re;
            FVVertex vtx11 = get_vtx(nic,njc);
            corner_coords[6] = vtx11.pos[0].x.re;
            corner_coords[7] = vtx11.pos[0].y.re;
            corner_coords[8] = vtx11.pos[0].z.re;
            FVVertex vtx01 = get_vtx(0,njc);
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
            FVVertex vtx000 = get_vtx(0,0,0);
            corner_coords[0] = vtx000.pos[0].x.re;
            corner_coords[1] = vtx000.pos[0].y.re;
            corner_coords[2] = vtx000.pos[0].z.re;
            FVVertex vtx100 = get_vtx(nic,0,0);
            corner_coords[3] = vtx100.pos[0].x.re;
            corner_coords[4] = vtx100.pos[0].y.re;
            corner_coords[5] = vtx100.pos[0].z.re;
            FVVertex vtx110 = get_vtx(nic,njc,0);
            corner_coords[6] = vtx110.pos[0].x.re;
            corner_coords[7] = vtx110.pos[0].y.re;
            corner_coords[8] = vtx110.pos[0].z.re;
            FVVertex vtx010 = get_vtx(0,njc,0);
            corner_coords[9] = vtx010.pos[0].x.re;
            corner_coords[10] = vtx010.pos[0].y.re;
            corner_coords[11] = vtx010.pos[0].z.re;
            FVVertex vtx001 = get_vtx(0,0,nkc);
            corner_coords[12] = vtx001.pos[0].x.re;
            corner_coords[13] = vtx001.pos[0].y.re;
            corner_coords[14] = vtx001.pos[0].z.re;
            FVVertex vtx101 = get_vtx(nic,0,nkc);
            corner_coords[15] = vtx101.pos[0].x.re;
            corner_coords[16] = vtx101.pos[0].y.re;
            corner_coords[17] = vtx101.pos[0].z.re;
            FVVertex vtx111 = get_vtx(nic,njc,nkc);
            corner_coords[18] = vtx111.pos[0].x.re;
            corner_coords[19] = vtx111.pos[0].y.re;
            corner_coords[20] = vtx111.pos[0].z.re;
            FVVertex vtx011 = get_vtx(0,njc,nkc);
            corner_coords[21] = vtx011.pos[0].x.re;
            corner_coords[22] = vtx011.pos[0].y.re;
            corner_coords[23] = vtx011.pos[0].z.re;
        }
    } // end copy_current_corner_coords()

    @nogc void set_current_corner_coords_to_infinity()
    {
        foreach(ref double cc; corner_coords) { cc = double.infinity; }
    }

    override size_t[] get_cell_write_indices() {
        size_t[] index;
        size_t[3] ijk;
        bool in_row, in_column, in_page;

        if ((myConfig.nic_write == 1) && (myConfig.njc_write == 1) && (myConfig.nkc_write == 1)) {
            index = iota(0, nic * njc * nkc).array();
        } else {
            foreach (indx; 0 .. (nic * njc * nkc)) {
                ijk = to_ijk_indices_for_cell(indx);
                // Lengthy conditional check here- spread it over a few lines
                // Either the cell falls on the specified indices, or is the last one in a specific row/column/page
                in_row = (ijk[0] % myConfig.nic_write) == 0;
                in_column = (ijk[1] % myConfig.njc_write) == 0;
                in_page = (ijk[2] % myConfig.nkc_write) == 0;
                if ((in_row && in_column && in_page) || // The specific indices condition, now for the end of row/column/page condition
                    (in_column && in_page && (ijk[0] == (nic - 1))) || (in_column && in_page && (ijk[1] == (njc - 1))) || (in_row && in_column && (ijk[0] == (nkc - 1)))) {
                    index ~= indx;
                }
            }
        }
        return index;
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
