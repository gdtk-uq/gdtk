// ufluidblock.d
// Class for unstructured blocks of cells, for use within Eilmer4.
// Peter J. 2014-11-07 first serious cut.

module ufluidblock;

import std.conv;
import std.file;
import std.json;
import std.stdio;
import std.format;
import std.string;
import std.array;
import std.math;
import std.algorithm;
import std.range;

import util.lua;
import util.lua_service;
import gas.luagas_model;
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
import fvvertex;
import fvinterface;
import fvcell;
import lsqinterp;
import fluidblock;
import fluidblockio_old;
import bc;
import grid_motion;
import grid_motion_udf;
import geom.luawrap.luausgrid;
import luaflowstate;
import fluidblockio_new;
import nm;
import user_defined_source_terms;


class UFluidBlock: FluidBlock {
public:
    size_t nvertices;
    size_t nboundaries;
    UnstructuredGrid grid;
    size_t ninteriorfaces;
    FVCell[] ghost_cells;
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.

public:
    this(in int id, JSONValue json_data)
    {
        label = getJSONstring(json_data, "label", "");
        ncells = getJSONint(json_data, "ncells", 0);
        size_t nghost = GlobalConfig.use_structured_reconstruction ? 2 : 1;
        super(id, Grid_t.unstructured_grid, ncells, nghost, label);
        nvertices = getJSONint(json_data, "nvertices", 0);
        nfaces = getJSONint(json_data, "nfaces", 0);
        nboundaries = getJSONint(json_data, "nboundaries", 0);
        active = getJSONbool(json_data, "active", true);
        omegaz = getJSONdouble(json_data, "omegaz", 0.0);
        may_be_turbulent = getJSONbool(json_data, "may_be_turbulent", true);
    } // end constructor from json

    this(lua_State* L)
    // Construct from a Lua interpreter state.
    // This particular constructor is only used in the prep stage.
    // Note that we assume the FlowState fs to be in a nonrotating frame.
    {
        auto grid = checkUnstructuredGrid(L, 2);
        double omegaz = luaL_checknumber(L, 5);
        ncells = grid.ncells;
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
        string msg = "get_unstructured_grid_fv(): ";
        Vector3 pos;
        number volume, xyplane_area, iLength, jLength, kLength, L_min;
        foreach(cell_idx, ref cell; cells) {
            // collect the vertices for this cell
            Vector3[] vertices;
            foreach (j; grid.cells[cell_idx].vtx_id_list) {
                vertices ~= grid.vertices[j];
            }
            if (GlobalConfig.dimensions == 2 ) {
                switch (vertices.length) {
                case 3:
                    xyplane_triangle_cell_properties(vertices[0], vertices[1], vertices[2],
                                                     pos, xyplane_area, iLength, jLength, L_min);
                    break;
                case 4:
                    xyplane_quad_cell_properties(vertices[0], vertices[1],
                                                 vertices[2], vertices[3],
                                                 pos, xyplane_area, iLength, jLength, L_min);
                    break;
                default:
                    debug { msg ~= format("Unhandled number of vertices: %d", vertices.length); }
                    throw new FlowSolverException(msg);
                } // end switch
                volume = xyplane_area * ((GlobalConfig.axisymmetric) ? pos.y : to!number(1.0) );
            } else if (GlobalConfig.dimensions == 3 ) {
                switch (vertices.length) {
                case 4:
                    tetrahedron_properties(vertices[0], vertices[1],
                                           vertices[2], vertices[3],
                                           pos, volume, L_min);
                    break;
                case 8:
                    hex_cell_properties(vertices[0], vertices[1], vertices[2], vertices[3],
                                        vertices[4], vertices[5], vertices[6], vertices[7],
                                        GlobalConfig.true_centroids, pos, volume, iLength, jLength, kLength);
                    break;
                case 5:
                    pyramid_properties(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4],
                                       GlobalConfig.true_centroids, pos, volume, L_min);
                    break;
                case 6:
                    wedge_properties(vertices[0], vertices[1], vertices[2],
                                     vertices[3], vertices[4], vertices[5],
                                     GlobalConfig.true_centroids, pos, volume, L_min);
                    break;
                default:
                    debug { msg ~= format("Unhandled number of vertices: %d", vertices.length); }
                    throw new FlowSolverException(msg);
                }
            }
            if (omegaz != 0.0) {
                throw new Error("Oops, we have not yet implemented rotating-frame code here.");
            }
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
                    errMsg ~= "The returned object is not a proper _FlowState handle or suitable table.";
                    luaL_error(L, errMsg.toStringz);
                }
            }
            celldata.flowstates ~= *myfs; // Copy the myfs flowstate into the celldata structure
            cells[cell_idx] = new FVCell(myConfig, pos, &celldata, volume, to!int(cell_idx));
        }
        block_io = get_fluid_block_io(this);
        if (lua_fs) { lua_settop(L, 0); }
    } // end constructor from Lua state

    override JSONValue get_header()
    // return information in JSON format that describes this block
    {
        JSONValue header = ["structured": false];
        header["nic"] = to!int(ncells);
        header["njc"] = 1;
        header["nkc"] = 1;

        return header;
    }

    override void init_workspace()
    {
        super.init_workspace();
        // Workspace for flux_calc method.
    }

    @nogc override int get_interpolation_order()
    {
        return myConfig.interpolation_order;
    }

    @nogc override void set_interpolation_order(int order)
    {
        myConfig.interpolation_order = order;
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
        // Although we make the helper functions available within
        // the block-specific Lua interpreter, we should use
        // those functions only in the context of the master thread.
        setSampleHelperFunctions(myL);
        setGridMotionHelperFunctions(myL);
        lua_pushinteger(myL, ncells); lua_setglobal(myL, "ncells");
        lua_pushinteger(myL, nvertices); lua_setglobal(myL, "nvertices");
        lua_pushinteger(myL, nfaces); lua_setglobal(myL, "nfaces");
        lua_pushinteger(myL, nboundaries); lua_setglobal(myL, "nboundaries");
        lua_pushinteger(myL, Face.north); lua_setglobal(myL, "north");
        lua_pushinteger(myL, Face.east); lua_setglobal(myL, "east");
        lua_pushinteger(myL, Face.south); lua_setglobal(myL, "south");
        lua_pushinteger(myL, Face.west); lua_setglobal(myL, "west");
        lua_pushinteger(myL, Face.top); lua_setglobal(myL, "top");
        lua_pushinteger(myL, Face.bottom); lua_setglobal(myL, "bottom");
        lua_pushinteger(myL, n_ghost_cell_layers); lua_setglobal(myL, "n_ghost_cell_layers");
    } // end init_lua_globals()

    override void init_boundary_conditions(JSONValue json_data)
    // Initialize boundary conditions after the blocks are constructed,
    // because we want access to the full collection of valid block references.
    {
        foreach (boundary; 0 .. nboundaries) {
            string json_key = format("boundary_%d", boundary);
            auto bc_json_data = json_data[json_key];
            bc ~= make_BC_from_json(bc_json_data, id, to!int(boundary));
        }
    } // end init_boundary_conditions()

    override string toString() const
    {
        char[] repr;
        repr ~= "UFluidBlock(unstructured_grid, ";
        repr ~= "id=" ~ to!string(id);
        repr ~= " label=\"" ~ label ~ "\"";
        repr ~= ", active=" ~ to!string(active);
        repr ~= ", grid_type=\"" ~ gridTypeName(grid_type) ~ "\"";
        repr ~= ", omegaz=" ~ to!string(omegaz);
        repr ~= ", may_be_turbulent=" ~ to!string(may_be_turbulent);
        repr ~= ", ncells=" ~ to!string(ncells);
        repr ~= ", nvertices=" ~ to!string(nvertices);
        repr ~= ", nfaces=" ~ to!string(nfaces);
        repr ~= ", \n    bc=[b_" ~ to!string(0) ~ "=" ~ to!string(bc[0]);
        foreach (i; 1 .. bc.length) {
            repr ~= ",\n        b_" ~ to!string(i) ~ "=" ~ to!string(bc[i]);
        }
        repr ~= "\n       ]"; // end bc list
        repr ~= ")";
        return to!string(repr);
    }

    @nogc
    override void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
        grid.find_enclosing_cell(p, indx, found); // delegate to the grid object
    }

    override void init_grid_and_flow_arrays(string gridFileName)
    {
        grid = new UnstructuredGrid(gridFileName, myConfig.grid_format, myConfig.true_centroids);
        if (grid.nvertices != nvertices) {
            string msg = format("UnstructuredGrid: incoming grid has %d vertices " ~
                                "but expected %d vertices.", grid.nvertices, nvertices);
            throw new FlowSolverException(msg);
        }
        if (grid.nfaces != nfaces) {
            string msg = format("UnstructuredGrid: incoming grid has %d faces " ~
                                "but expected %d faces.", grid.nfaces, nfaces);
            throw new FlowSolverException(msg);
        }
        if (grid.ncells != ncells) {
            string msg = format("UnstructuredGrid: incoming grid has %d cells " ~
                                "but expected %d cells.", grid.ncells, ncells);
            throw new FlowSolverException(msg);
        }
        grid.sort_cells_into_bins();
        ninteriorfaces = grid.ninteriorfaces;

        // Assemble array storage for finite-volume cells, etc.
        bool lsq_workspace_at_vertices = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.vertices);
        foreach (i, v; grid.vertices) {
            auto new_vtx = new FVVertex(myConfig, lsq_workspace_at_vertices, to!int(i));
            new_vtx.pos[0] = v;
            vertices ~= new_vtx;
        }

        // Allocate densified face and cell data
        auto gmodel = myConfig.gmodel;
        size_t neq    = myConfig.cqi.n;
        size_t nsp    = myConfig.n_species;
        size_t nmodes = myConfig.n_modes;
        size_t nturb  = myConfig.turb_model.nturb;
        size_t nftl   = myConfig.n_flow_time_levels;

        bool lsq_workspace_at_cells = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.cells);
        bool lsq_workspace_at_faces = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces);

        size_t nghost = 0;
        foreach (bndry; grid.boundaries) nghost += bndry.face_id_list.length;

        nghost *= n_ghost_cell_layers;

        allocate_dense_celldata(ncells, nghost, neq, nftl);
        allocate_dense_facedata(nfaces, nghost, neq, nftl);

        // We have some unstructured specific stuff that also needs to be allocated,
        // though only when using the fully unstructured convective gradients
        if (n_ghost_cell_layers<2){
            celldata.face_distances.length = ncells;
            celldata.lsqws.length = ncells + nghost;
            celldata.lsqgradients.reserve(ncells + nghost);
            foreach (i; 0 .. ncells + nghost) celldata.lsqgradients ~= LSQInterpGradients(nsp, nmodes, nturb); // TODO: skip if first order
            facedata.dL.length = nfaces;
            facedata.dR.length = nfaces;
            celldata.c2v.length = ncells;
        }

        if (myConfig.unstructured_limiter == UnstructuredLimiter.venkat_mlp) {
            vertexdata.n_indices.length = grid.vertices.length;
            vertexdata.cell_cloud_indices.length = grid.vertices.length;
        }

        foreach (i, f; grid.faces) {
            auto new_face = new FVInterface(myConfig, IndexDirection.none, &facedata, to!int(i));
            faces ~= new_face;
        }

        cells.reserve(ncells);
        ghost_cells.reserve(nghost);
        foreach (i, c; grid.cells) {
            // Note that the cell id and the index in the cells array are the same.
            // We will reply upon this connection in other parts of the flow code.
            auto new_cell = new FVCell(myConfig, &celldata, to!int(i));
            new_cell.contains_flow_data = true;
            new_cell.is_interior_to_domain = true;
            cells ~= new_cell;
        }
        // Bind the interfaces, vertices and cells together,
        // using the indices stored in the unstructured grid.
        foreach (i, f; faces) {
            f.left_cells.length = n_ghost_cell_layers;
            f.right_cells.length = n_ghost_cell_layers;
            foreach (j; grid.faces[i].vtx_id_list) {
                f.vtx ~= vertices[j];
            }
        }
        foreach (i, c; cells) {
            foreach (j; grid.cells[i].vtx_id_list) {
                c.vtx ~= vertices[j];
            }
            auto nf = grid.cells[i].face_id_list.length;
            if (nf != grid.cells[i].outsign_list.length) {
                string msg = format("Mismatch in face_id_list, outsign_list lengths: %d %d",
                                    grid.cells[i].face_id_list.length,
                                    grid.cells[i].outsign_list.length);
                throw new FlowSolverException(msg);
            }
            foreach (j; 0 .. nf) {
                auto my_face = faces[grid.cells[i].face_id_list[j]];
                int my_outsign = grid.cells[i].outsign_list[j];
                c.iface ~= my_face;
                celldata.c2f[i] ~= grid.cells[i].face_id_list[j];
                c.outsign ~= my_outsign;
                celldata.outsigns[i] ~= my_outsign;
                if (my_outsign == 1) {
                    if (my_face.left_cell) {
                        string msg = format("Already have cell %d attached to left-of-face %d. Attempt to add cell %d.",
                                            my_face.left_cell.id, my_face.id, c.id);
                        throw new FlowSolverException(msg);
                    } else {
                        my_face.left_cell = c;
                        my_face.left_cells[0] = c;
                        facedata.f2c[grid.cells[i].face_id_list[j]].left = i;
                    }
                } else {
                    if (my_face.right_cell) {
                        string msg = format("Already have cell %d attached to right-of-face %d. Attempt to add cell %d.",
                                            my_face.right_cell.id, my_face.id, c.id);
                        throw new FlowSolverException(msg);
                    } else {
                        my_face.right_cell = c;
                        my_face.right_cells[0] = c;
                        facedata.f2c[grid.cells[i].face_id_list[j]].right = i;
                    }
                }
            }
            celldata.nfaces[i] = celldata.c2f[i].length;
        } // end foreach cells
        // If we are doing fully unstructured reconstruction, we need to set up some dense data
        if (n_ghost_cell_layers<2) {
            foreach (i, c; cells) {
                foreach (j; grid.cells[i].vtx_id_list) {
                    celldata.c2v[i] ~= j;
                }
                celldata.face_distances[i].length = celldata.nfaces[i];
            }
        }
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
        //
        // Work through the faces on the boundaries and add ghost cells.
        if (nboundaries != grid.nboundaries) {
            string msg = format("Mismatch in number of boundaries: %d %d",
                                nboundaries, grid.nboundaries);
            throw new FlowSolverException(msg);
        }
        int ghost_cell_id = to!int(cells.length);
        foreach (i, bndry; grid.boundaries) {
            auto nf = bndry.face_id_list.length;
            if (nf != bndry.outsign_list.length) {
                string msg = format("Mismatch in face_id_list, outsign_list lengths: %d %d",
                                    bndry.face_id_list.length, bndry.outsign_list.length);
                throw new FlowSolverException(msg);
            }
            foreach (j; 0 .. nf) {
                FVInterface my_face = faces[bndry.face_id_list[j]];
                my_face.is_on_boundary = true;
                my_face.bc_id = i; // note which boundary this face is on
                int my_outsign = bndry.outsign_list[j];
                bc[i].faces ~= my_face;
                bc[i].outsigns ~= my_outsign;
                my_face.i_bndry = bc[i].outsigns.length - 1;
                if (bc[i].ghost_cell_data_available) {
                    FVCell ghost0 = new FVCell(myConfig, &celldata, ghost_cell_id);
                    ghost_cell_id++;
                    ghost0.contains_flow_data = bc[i].ghost_cell_data_available;
                    ghost0.is_ghost_cell = true;
                    bc[i].ghostcells ~= ghost0;
                    ghost_cells ~= ghost0;
                    if (my_outsign == 1) {
                        if (my_face.right_cell) {
                            string msg = format("Already have cell %d attached to right-of-face %d."
                                                ~" Attempt to add ghost cell %d.",
                                                my_face.right_cell.id, my_face.id, ghost0.id);
                            throw new FlowSolverException(msg);
                        } else {
                            my_face.right_cell = ghost0;
                            my_face.right_cells[0] = ghost0;
                            facedata.f2c[bndry.face_id_list[j]].right = ghost0.id;
                        }
                    } else {
                        if (my_face.left_cell) {
                            string msg = format("Already have cell %d attached to left-of-face %d."
                                                ~" Attempt to add ghost cell %d.",
                                                my_face.left_cell.id, my_face.id, ghost0.id);
                            throw new FlowSolverException(msg);
                        } else {
                            my_face.left_cell = ghost0;
                            my_face.left_cells[0] = ghost0;
                            facedata.f2c[bndry.face_id_list[j]].left = ghost0.id;
                        }
                    }
                    // For nick's experimental pseudo-structured stencils, we need a second layer of
                    // ghost cells. Put them in the stencils here, for now.
                    if (n_ghost_cell_layers>1) {
                        FVCell ghost1 = new FVCell(myConfig, &celldata, ghost_cell_id);
                        ghost_cell_id++;
                        ghost1.contains_flow_data = bc[i].ghost_cell_data_available;
                        ghost1.is_ghost_cell = true;
                        bc[i].ghostcells ~= ghost1;
                        ghost_cells ~= ghost1;
                        size_t fid = bndry.face_id_list[j];
                        if (my_outsign == 1) {
                            facedata.stencil_idxs[fid].R1 = ghost1.id;
                            my_face.right_cells[1] = ghost1;
                        } else {
                            facedata.stencil_idxs[fid].L1 = ghost1.id;
                            my_face.left_cells[1] = ghost1;
                        }
                    }
                } // end if (bc[i].ghost_cell_data_available
            } // end foreach j
        } // end foreach i
        // At this point, all faces should have either one finite-volume cell
        // or one ghost cell attached to each side -- check that this is true.

        // Setup dense arrays for fast checks of left or right data available (NNG)
        foreach (i, bndry; grid.boundaries) {

            // For any other kind of boundary we need to mark which side the interior cell is
            auto nf = bndry.face_id_list.length;
            foreach (j; 0 .. nf) {
                size_t my_id = faces[bndry.face_id_list[j]].id;
                size_t bid = faces[bndry.face_id_list[j]].id;
                int my_outsign = bndry.outsign_list[j];

                // For a shared boundary both sides are okay, so we skip
                if (startsWith(bc[i].type, "exchange_")) { continue; }

                if (my_outsign == 1) { // The ghost cell is a right cell
                    facedata.left_interior_only[bid] = true;
                } else {               // the ghost cell is a left cell
                    facedata.right_interior_only[bid] = true;
                }
            }
        }
        foreach (f; faces) {
            bool ok = true;
            string msg = " ";
            if (f.is_on_boundary) {
                if (bc[f.bc_id].ghost_cell_data_available) {
                    if (f.left_cell && f.right_cell) {
                        ok = true;
                    } else {
                        ok = false;
                        msg ~= "Boundary face does not have a cell attached to each side.";
                    }
                } else {
                    if ((f.left_cell is null) != (f.right_cell is null)) {
                        ok = true;
                    } else {
                        ok = false;
                        msg ~= "Boundary face does not have exactly one cell attached.";
                    }
                }
            } else {
                // not on a boundary, should have one interior cell per side.
                if (f.left_cell && f.right_cell) {
                    ok = true;
                } else {
                    ok = false;
                    msg ~= "Non-boundary face does not have a cell attached to each side.";
                }
            }
            if (!ok) {
                msg = format("After adding ghost cells to face %d: ", f.id) ~ msg;
                writeln("Oops... ", msg);
                writeln("grid.faces[f.id].vtx_id_list=", grid.faces[f.id].vtx_id_list);
                if (f.left_cell) writeln("f.left_cell=", f.left_cell); else writeln("no left cell");
                if (f.right_cell) writeln("f.right_cell=", f.right_cell); else writeln("no right cell");
                throw new FlowSolverException(msg);
            }
        } // end foreach f

        //
        // Set up the cell clouds for least-squares derivative estimation for use in
        // interpolation/reconstruction of flow quantities at left- and right- sides
        // of cell faces.
        // (Will be used for the convective fluxes).
        if (myConfig.use_extended_stencil) {
            string err = "Extended cell clouds are not enabled by default for memory reasons.";
            err       ~= "Please change cloud_nmax in lsqinterp to 112 and recompile.";
            throw new Error(err);
            //allocate_extended_cell_cloud_references();
        } else {
            allocate_regular_cell_cloud_references();
        } // end else
        // We will now store the cloud of points in cloud_pos for viscous derivative calcualtions.
        // This is equivalent to store_references_for_derivative_calc(size_t gtl) in sblock.d
        if (myConfig.viscous) {
            if (myConfig.spatial_deriv_calc ==  SpatialDerivCalc.divergence) {
                throw new Error("Divergence theorem not implemented for unstructured grid");
            }
            // else continue on to fill least-squares/finite-difference cloud
            final switch (myConfig.spatial_deriv_locn) {
            case SpatialDerivLocn.vertices:
                throw new Error("spatial_deriv_locn at vertices not implemented for unstructured grid");
                // no need for break;
            case SpatialDerivLocn.faces:
                throw new Error("spatial_deriv_locn at interfaces not implemented for unstructured grid");
                // no need for break;
            case SpatialDerivLocn.cells:
                // We will now store the cloud of points in cloud_pos for viscous derivative calcualtions.
                foreach (c; cells) {
                    // First cell in the cloud is the cell itself.  Differences are taken about it.
                    c.cloud_pos ~= &(c.pos[0]);
                    c.cloud_fs ~= c.fs;
                    // Subsequent cells are the surrounding cells.
                    foreach (i, f; c.iface) {
                        c.cloud_pos ~= &(f.pos);
                        c.cloud_fs ~= f.fs;
                    } // end foreach face
                }
            } // end switch (myConfig.spatial_deriv_locn)
        } // if myConfig.viscous

        // Experimental extra stencil size for using the structured style reconstruction
        if (n_ghost_cell_layers>1) {
            foreach (fid; 0 .. nfaces) {
                auto face = faces[fid];
                auto left = face.left_cells[0];
                auto right= face.right_cells[0];
                // Even in an unstructured grid, we know that the vertices in a cell
                // have a required order. This means we can find the opposite face
                // by checking its vertices.

                if (left.id<ncells){
                    FVInterface oface = left.get_opposite_face(face);

                    if (oface.right_cells[0].id == left.id) {
                        face.left_cells[1]            = oface.left_cells[0];
                    } else if (oface.left_cells[0].id == left.id ){
                        face.left_cells[1]            = oface.right_cells[0];
                    } else {
                        throw new Error("Error fixing second left cell for face.");
                    }
                }

                if (right.id<ncells){
                    FVInterface oface = right.get_opposite_face(face);
                    if (oface.left_cells[0].id == right.id) {
                        face.right_cells[1]           = oface.right_cells[0];
                    } else if (oface.right_cells[0].id == right.id ){
                        face.right_cells[1]           = oface.left_cells[0];
                    } else {
                        throw new Error("Error fixing second right cell for face.");
                    }
                }
            }
            // The structured-style reconstruction requires this array to define its stencil.
            foreach (fid; 0 .. nfaces) {
                facedata.stencil_idxs[fid].L1 = faces[fid].left_cells[1].id;
                facedata.stencil_idxs[fid].L0 = faces[fid].left_cells[0].id;
                facedata.stencil_idxs[fid].R0 = faces[fid].right_cells[0].id;
                facedata.stencil_idxs[fid].R1 = faces[fid].right_cells[1].id;
            }
        } // End n_ghost_cell_layers>1 special init code.
    } // end init_grid_and_flow_arrays()

    void allocate_extended_cell_cloud_references()
    {
        foreach (c; cells) {
            // First cell in the cloud is the cell itself.  Differences are taken about it.
            c.cell_cloud ~= c;
            size_t[] cell_ids;
            cell_ids ~= c.id;
            // next we add the nearest neighbour cells
            // Note that we assume these cells are next in the array when creating the low order Jacobian
            foreach (i, f; c.iface) {
                if (c.outsign[i] == 1) {
                    if (f.right_cell && f.right_cell.contains_flow_data) {
                        c.cell_cloud ~= f.right_cell;
                        cell_ids ~= f.right_cell.id;
                        celldata.cell_cloud_indices[c.id] ~= f.right_cell.id;
                    }
                } else {
                    if (f.left_cell && f.left_cell.contains_flow_data) {
                        c.cell_cloud ~= f.left_cell;
                        cell_ids ~= f.left_cell.id;
                        celldata.cell_cloud_indices[c.id] ~= f.left_cell.id;
                    }
                }
            } // end foreach face
            // finally we add the remaining cells that share a node with the central cell
            bool is_on_boundary = false;
            foreach(face; c.iface) if (face.is_on_boundary) is_on_boundary = true;
            if (is_on_boundary) {
                // apply nearest-face neighbour
                foreach (i, f; c.iface) {
                    if (c.outsign[i] == 1) {
                        if (f.right_cell && cell_ids.canFind(f.right_cell.id) == false && f.right_cell.contains_flow_data) {
                            c.cell_cloud ~= f.right_cell;
                            cell_ids ~= f.right_cell.id;
                            celldata.cell_cloud_indices[c.id] ~= f.right_cell.id;
                        }
                    } else {
                        if (f.left_cell && cell_ids.canFind(f.left_cell.id) == false && f.left_cell.contains_flow_data) {
                            c.cell_cloud ~= f.left_cell;
                            cell_ids ~= f.left_cell.id;
                            celldata.cell_cloud_indices[c.id] ~= f.left_cell.id;
                        }
                    }
                } // end foreach face
            } else {
                // apply nearest-node neighbour
                foreach(vtx; c.vtx) {
                    foreach(cid; cellIndexListPerVertex[vtx.id])
                        if (cell_ids.canFind(cid) == false && cells[cid].contains_flow_data)
                            { c.cell_cloud ~= cells[cid]; cell_ids ~= cid; celldata.cell_cloud_indices[c.id] ~= cid; }
                }
            }
        } // end foreach cell
    }

    void allocate_regular_cell_cloud_references()
    {
        foreach (c; cells) {
            // First cell in the cloud is the cell itself.  Differences are taken about it.
            c.cell_cloud ~= c;
            // Subsequent cells are the surrounding cells.
            foreach (i, f; c.iface) {
                if (c.outsign[i] == 1) {
                    if (f.right_cell && f.right_cell.contains_flow_data) {
                        c.cell_cloud ~= f.right_cell;
                        celldata.cell_cloud_indices[c.id] ~= f.right_cell.id;
                    }
                } else {
                    if (f.left_cell && f.left_cell.contains_flow_data) {
                        c.cell_cloud ~= f.left_cell;
                        celldata.cell_cloud_indices[c.id] ~= f.left_cell.id;
                    }
                }
            } // end foreach face
        } // end foreach cell
    }

    void build_cloud_of_cell_references_at_each_vertex()
    {
        // For the MLP limiter, we need access to the gradients stored in the cells.
        foreach (vtx; vertices) {
            foreach(cid; cellIndexListPerVertex[vtx.id]) { vtx.cell_cloud ~= cells[cid]; }
        }

        if (myConfig.unstructured_limiter == UnstructuredLimiter.venkat_mlp) {
            foreach (i, vtx; vertices) {
                foreach(cid; cellIndexListPerVertex[vtx.id]) { vertexdata.cell_cloud_indices[i] ~= cells[cid].id; }
                vertexdata.n_indices[i] = vertexdata.cell_cloud_indices[i].length;
            }
        }
    }

    @nogc
    override void precompute_stencil_data(size_t gtl)
    {
    /*
        Used in the experimental quasi-structured reconstruction mode.
    */
        if (myConfig.interpolation_order!=2) return;

        // Having completed the LLRR stencil indices, we can now precompute the stencil
        // coefficients, using the cell sizes, kind of like how the structured code does it.
        // I should note that the lengths generated by this routine are slightly different
        // to the ones used by the real structured code. If this becomes a problem, we can
        // change it later.
        foreach (fid; 0 .. nfaces) {
            size_t L1 = facedata.stencil_idxs[fid].L1;
            size_t L0 = facedata.stencil_idxs[fid].L0;
            size_t R0 = facedata.stencil_idxs[fid].R0;
            size_t R1 = facedata.stencil_idxs[fid].R1;
            //                    dL0
            //                   <---
            //         <------dL1----
            //   .---------.---------.---------.---------.
            //   |         |         |         |         |
            //   |    o    |    o    |    o    |    o    |
            //   |         |         |         |         |
            //   .---------.---------.---------.---------.
            //    <-lenL1-> <-lenL0->

            double dL0 = distance_between(facedata.positions[fid], celldata.positions[L0]);
            double dL1 = distance_between(facedata.positions[fid], celldata.positions[L1]);
            number lenL0 = 2.0*dL0;
            number lenL1 = 2.0*(dL1-lenL0);

            double dR0 = distance_between(facedata.positions[fid], celldata.positions[R0]);
            double dR1 = distance_between(facedata.positions[fid], celldata.positions[R1]);

            number lenR0 = 2.0*dR0;
            number lenR1 = 2.0*(dR1-lenR0);
            facedata.l2r2_interp_data[fid].set(lenL1, lenL0, lenR0, lenR1);

        }
    }

    @nogc
    override void compute_primary_cell_geometric_data(size_t gtl)
    {
        if (myConfig.dimensions == 2) {
            foreach (c; cells) { c.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
            foreach (f; faces) { f.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
        } else { // 3D
            foreach (c; cells) { c.update_3D_geometric_data(gtl, myConfig.true_centroids); }
            foreach (f; faces) { f.update_3D_geometric_data(gtl); }
        }

        if (celldata.face_distances) {
            foreach(i, ifaces; celldata.c2f){
                foreach(j, jface; ifaces){
                    celldata.face_distances[i][j] = facedata.positions[jface] - celldata.positions[i];
                }
            }
        }
        //
        // Guess the position ghost-cell centres and copy cross-cell lengths.
        // Copy without linear extrapolation for the moment.
        //
        // 25-Feb-2014 Note copied from Eilmer3
        // Jason Qin and Paul Petrie-Repar have identified the lack of exact symmetry in
        // the reconstruction process at the wall as being a cause of the leaky wall
        // boundary conditions.  Note that the symmetry is not consistent with the
        // linear extrapolation used for the positions and volumes in Eilmer3.
        // 2018-06-10: Maybe we can fix some of this by eliminating ghost-cells for
        // boundary conditions that do not really have gas regions backing the
        // ghost-cells.
        foreach (i, bndry; grid.boundaries) {
            if (bc[i].ghost_cell_data_available == false) { continue; }
            auto nf = bndry.face_id_list.length;
            foreach (j; 0 .. nf) {
                auto my_face = faces[bndry.face_id_list[j]];
                auto my_outsign = bndry.outsign_list[j];
                if (my_outsign == 1) {
                    auto inside0 = my_face.left_cell;
                    Vector3 delta; delta = my_face.pos; delta -= inside0.pos[gtl];
                    auto ghost0 = my_face.right_cell;
                    ghost0.pos[gtl] = my_face.pos; ghost0.pos[gtl] += delta;
                    ghost0.iLength = inside0.iLength;
                    ghost0.jLength = inside0.jLength;
                    ghost0.kLength = inside0.kLength;
                    ghost0.L_min = inside0.L_min;
                    ghost0.L_max = inside0.L_max;
                    ghost0.update_celldata_geometry();

                    if (n_ghost_cell_layers>1){
                        // Sometimes we get situations where inside1 is actually a
                        // ghost cell, which might not have been init'd yet.
                        // These should get fixed below in update_nonshared_ghost_cell_positions
                        auto ghost1  = my_face.right_cells[1];
                        auto inside1 = my_face.left_cells[1];
                        delta = my_face.pos; delta -= inside1.pos[gtl];
                        ghost1.pos[gtl] = my_face.pos; ghost1.pos[gtl] += delta;
                        ghost1.iLength = ghost0.iLength;
                        ghost1.jLength = ghost0.jLength;
                        ghost1.kLength = ghost0.kLength;
                        ghost1.L_min = ghost0.L_min;
                        ghost1.L_max = ghost0.L_max;
                        ghost1.update_celldata_geometry();
                    }
                } else {
                    auto inside0 = my_face.right_cell;
                    Vector3 delta; delta = my_face.pos; delta -= inside0.pos[gtl];
                    auto ghost0 = my_face.left_cell;
                    ghost0.pos[gtl] = my_face.pos; ghost0.pos[gtl] += delta;
                    ghost0.iLength = inside0.iLength;
                    ghost0.jLength = inside0.jLength;
                    ghost0.kLength = inside0.kLength;
                    ghost0.L_min = inside0.L_min;
                    ghost0.L_max = inside0.L_max;
                    ghost0.update_celldata_geometry();

                    if (n_ghost_cell_layers>1){
                        auto ghost1  = my_face.left_cells[1];
                        auto inside1 = my_face.right_cells[1];
                        delta = my_face.pos; delta -= inside1.pos[gtl];
                        ghost1.pos[gtl] = my_face.pos; ghost1.pos[gtl] += delta;
                        ghost1.iLength = ghost0.iLength;
                        ghost1.jLength = ghost0.jLength;
                        ghost1.kLength = ghost0.kLength;
                        ghost1.L_min = ghost0.L_min;
                        ghost1.L_max = ghost0.L_max;
                        ghost1.update_celldata_geometry();
                    }
                } // end if my_outsign
            } // end foreach j
        } // end foreach bndry
    } // end compute_primary_cell_geometric_data()

    @nogc
    override void update_nonshared_ghost_cell_positions(size_t gtl)
    {
        // Jason Qin and Paul Petrie-Repar have identified the lack of exact symmetry in
        // the reconstruction process at the wall as being a cause of the leaky wall
        // boundary conditions.
        if (n_ghost_cell_layers<2) return;

        foreach (i, bndry; grid.boundaries) {
            if (bc[i].ghost_cell_data_available == false) { continue; }
            if (startsWith(bc[i].type, "exchange_")) { continue; }

            auto nf = bndry.face_id_list.length;
            foreach (j; 0 .. nf) {
                auto my_face = faces[bndry.face_id_list[j]];
                auto my_outsign = bndry.outsign_list[j];
                if (my_outsign == 1) {
                    auto ghost1  = my_face.right_cells[1];
                    auto inside1 = my_face.left_cells[1];

                    Vector3 delta; delta = my_face.pos; delta -= inside1.pos[gtl];
                    ghost1.pos[gtl] = my_face.pos; ghost1.pos[gtl] += delta;
                    ghost1.iLength = inside1.iLength;
                    ghost1.jLength = inside1.jLength;
                    ghost1.kLength = inside1.kLength;
                    ghost1.L_min = inside1.L_min;
                    ghost1.L_max = inside1.L_max;
                    ghost1.update_celldata_geometry();
                } else {
                    auto ghost1  = my_face.left_cells[1];
                    auto inside1 = my_face.right_cells[1];

                    Vector3 delta; delta = my_face.pos; delta -= inside1.pos[gtl];
                    ghost1.pos[gtl] = my_face.pos; ghost1.pos[gtl] += delta;
                    ghost1.iLength = inside1.iLength;
                    ghost1.jLength = inside1.jLength;
                    ghost1.kLength = inside1.kLength;
                    ghost1.L_min = inside1.L_min;
                    ghost1.L_max = inside1.L_max;
                    ghost1.update_celldata_geometry();
                } // end if my_outsign
            } // end foreach j
        } // end foreach bndry
    } // end compute_primary_cell_geometric_data()

    @nogc
    override void compute_least_squares_setup(size_t gtl)
    {
        // Update the least-squares geometric weights and the workspaces, if appropriate.
        // The weights should be calculated when the grid is initialised or moved.
        //
        if (myConfig.viscous && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)) {
            // LSQ weights are used in the calculation of flow gradients for the viscous terms.
            if (myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                foreach(cell; cells) {
                    cell.grad.set_up_workspace_leastsq(myConfig, cell.cloud_pos, cell.pos[gtl], false, *(cell.ws_grad));
                }
            } else if (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces) {
                foreach(iface; faces) {
                    iface.grad.set_up_workspace_leastsq(myConfig, iface.cloud_pos, iface.pos, false, *(iface.ws_grad));
                }
            } else { // myConfig.spatial_deriv_locn == vertices
                foreach(vtx; vertices) {
                    vtx.grad.set_up_workspace_leastsq(myConfig, vtx.cloud_pos, vtx.pos[gtl], true, *(vtx.ws_grad));
                }
            }
        }
        // The LSQ linear model for the flow field reconstruction is fitted using information
        // on the locations of the points.
        // This model is used later, as part of the convective flux calculation.
        if ((myConfig.interpolation_order > 1) && (n_ghost_cell_layers<2)) {
            foreach (c; cells) {
                try {
                    c.ws.assemble_and_invert_normal_matrix(c.cell_cloud, myConfig.dimensions, gtl);
                } catch (Exception e) {
                    debug {
                        writefln("In compute_least_squares_setup()," ~
                                 " we have failed to assemble and invert normal matrix for cell id=%d",
                                 c.id);
                    }
                    throw e;
                }
            }
            // Densified lsqgradients needs distances from the face to the left and right cells
            foreach(idx; 0 .. nfaces){              // TODO: What about no ghost boundaries?
                size_t l = facedata.f2c[idx].left;
                size_t r = facedata.f2c[idx].right;
                facedata.dL[idx] = facedata.positions[idx] - celldata.positions[l];
                facedata.dR[idx] = facedata.positions[idx] - celldata.positions[r];
            } // Done interior faces, next we do the boundaries of this block
        }
    } // end compute_least_squares_setup()

    @nogc
    override void sync_vertices_from_underlying_grid(size_t gtl=0)
    {
        foreach (i; 0 .. vertices.length) { vertices[i].pos[gtl].set(grid[i]); }
    }

    @nogc
    override void sync_vertices_to_underlying_grid(size_t gtl=0)
    {
        foreach (i; 0 .. vertices.length) { grid[i].set(vertices[i].pos[gtl]); }
    }

    override void read_new_underlying_grid(string fileName)
    {
        if (myConfig.verbosity_level > 1) { writeln("read_new_underlying_grid() for block ", id); }
        grid = new UnstructuredGrid(fileName, myConfig.grid_format, myConfig.true_centroids);
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
        throw new FlowSolverException("function not implemented for unstructured grid.");
    }

    @nogc
    void prepare_unstructured_gradients(bool allow_high_order_interpolation, size_t[] cell_idxs, size_t[] face_idxs)
    {
        if (!allow_high_order_interpolation) return;
        if (myConfig.interpolation_order == 1) return;

        if (cell_idxs.length==0) cell_idxs = celldata.all_cell_idxs;

        immutable bool needs_pressure_gradient =
                myConfig.unstructured_limiter == UnstructuredLimiter.park ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvan_albada ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvenkat ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp ||
                myConfig.unstructured_limiter == UnstructuredLimiter.hnishikawa;
        immutable size_t nsp = myConfig.n_species;
        immutable size_t nmodes = myConfig.n_modes;
        immutable size_t nturb = myConfig.turb_model.nturb;
        immutable size_t is3d = myConfig.dimensions==3;

        if (myConfig.unstructured_limiter == UnstructuredLimiter.venkat_mlp) {
            foreach(cid; cell_idxs){
                celldata.lsqgradients[cid].reset_max_min_values(celldata.flowstates[cid], nsp, nmodes, nturb, myConfig);

                foreach(vid; celldata.c2v[cid]){
                    foreach(ciid; vertexdata.cell_cloud_indices[vid]){
                        celldata.lsqgradients[cid].accumulate_max_min_values(celldata.flowstates[ciid], nsp, nmodes, nturb, myConfig);
                    }
                }
            }
        } else {
            foreach(cid; cell_idxs){
                celldata.lsqgradients[cid].reset_max_min_values(celldata.flowstates[cid], nsp, nmodes, nturb, myConfig);
                foreach(ciid; celldata.cell_cloud_indices[cid]){
                    celldata.lsqgradients[cid].accumulate_max_min_values(celldata.flowstates[ciid], nsp, nmodes, nturb, myConfig);
                }
            }
        }

        foreach(cid; cell_idxs){
            celldata.lsqgradients[cid].compute_lsq_values(celldata.flowstates[cid], celldata.lsqws[cid],
                                                          celldata.flowstates, celldata.cell_cloud_indices[cid],
                                                          myConfig, nsp, nmodes, nturb, is3d, needs_pressure_gradient);
        }
        return;
    }

    @nogc
    override void convective_flux_phase0new(bool allow_high_order_interpolation, size_t[] cell_idxs=[], size_t[] face_idxs=[])
    {
        if (n_ghost_cell_layers > 1) {
            // Use the structured stencils to do the inviscid fluxes
            if (face_idxs.length==0) face_idxs = facedata.all_face_idxs;

            bool second_order = allow_high_order_interpolation && (myConfig.interpolation_order == 2);
            if (second_order) {
                second_order_flux_calc(0, face_idxs);
            } else {
                first_order_flux_calc(0, face_idxs);
            }
        } else {
            // Compute gradients of flow quantities for higher-order reconstruction, if required.
            // To be used, later, in the convective flux calculation.
            prepare_unstructured_gradients(allow_high_order_interpolation, cell_idxs, face_idxs);
        }
        return;

    } // end convective_flux-phase0()

    @nogc
    override void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl,
                                         FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
    // Compute gradients of flow quantities for higher-order reconstruction, if required.
    // To be used, later, in the convective flux calculation.
    {

        if (cell_list.length == 0) { cell_list = cells; }
        if (vertex_list.length == 0) { vertex_list = vertices; }

        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            foreach (c; cell_list) {
                if (myConfig.unstructured_limiter == UnstructuredLimiter.venkat_mlp) {
                    c.gradients.store_max_min_values_for_extended_stencil(c.cell_cloud, myConfig);
                } else {
                    c.gradients.store_max_min_values_for_compact_stencil(c.cell_cloud, myConfig);
                }
                c.gradients.compute_lsq_values(c.cell_cloud, *(c.ws), myConfig);
            }
        } // end if interpolation_order > 1
    } // end convective_flux-phase0()

    @nogc
    override void convective_flux_phase1new(bool allow_high_order_interpolation, size_t[] cell_idxs=[], size_t[] face_idxs=[])
        // Compute limiter values of flow quantities for higher-order reconstruction, if required.
        // To be used, later, in the convective flux calculation.
    {
        if (!allow_high_order_interpolation) return;
        if (myConfig.interpolation_order == 1) return;
        if (GlobalConfig.frozen_limiter) return;
        if (n_ghost_cell_layers > 1) return;

        if (cell_idxs.length==0) cell_idxs = celldata.all_cell_idxs;

        switch (myConfig.unstructured_limiter) {
            case UnstructuredLimiter.svan_albada:
                // do nothing now
                break;
            case UnstructuredLimiter.min_mod:
                // do nothing now
                break;
            case UnstructuredLimiter.barth:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.barth_limit(c.cell_cloud, *(c.ws), myConfig); }
                break;
            case UnstructuredLimiter.park:
                immutable bool is3d = myConfig.dimensions == 3;
                foreach (i; cell_idxs) {
                    celldata.lsqgradients[i].park_limit2(celldata.lsqgradients[i], celldata.flowstates, i, facedata.f2c, is3d,
                                                         celldata.face_distances[i], celldata.c2f[i], celldata.nfaces[i], myConfig);
                }
                break;
            case UnstructuredLimiter.hvan_albada:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.van_albada_limit(c.cell_cloud, *(c.ws), true, myConfig); }
                break;
            case UnstructuredLimiter.van_albada:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.van_albada_limit(c.cell_cloud, *(c.ws), false, myConfig); }
                break;
            case UnstructuredLimiter.hnishikawa:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.nishikawa_limit(c.cell_cloud, *(c.ws), true, myConfig, 0); }
                break;
            case UnstructuredLimiter.nishikawa:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.nishikawa_limit(c.cell_cloud, *(c.ws), false, myConfig, 0); }
                break;
            case UnstructuredLimiter.hvenkat_mlp:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.venkat_mlp_limit(c.cell_cloud, *(c.ws), true, myConfig, 0); }
                break;
            case UnstructuredLimiter.venkat_mlp:
                foreach (i; cell_idxs) { auto c = cells[i]; c.gradients.venkat_mlp_limit(c.cell_cloud, *(c.ws), false, myConfig, 0); }
                break;
            case UnstructuredLimiter.hvenkat:
                immutable bool is3d = myConfig.dimensions == 3;
                foreach (i; cell_idxs) {
                    number phi_hp = park_equation(celldata.lsqgradients[i], celldata.flowstates, i, facedata.f2c, is3d,
                                                  celldata.face_distances[i], celldata.c2f[i], celldata.nfaces[i]);
                    celldata.lsqgradients[i].venkat_limit2(celldata.flowstates[i],
                        celldata.volumes[i], celldata.face_distances[i],
                        celldata.nfaces[i], phi_hp, myConfig);
                }
                break;
            case UnstructuredLimiter.venkat:
                number phi_hp = to!number(1.0);
                foreach (i; cell_idxs) {
                    celldata.lsqgradients[i].venkat_limit2(celldata.flowstates[i],
                        celldata.volumes[i], celldata.face_distances[i],
                        celldata.nfaces[i], phi_hp, myConfig);
                }
                break;
            default:
                throw new Error("Bad limiter selected");
        }
    } // end convective_flux-phase1()

    @nogc
    override void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0,
                                         FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
        // Compute limiter values of flow quantities for higher-order reconstruction, if required.
        // To be used, later, in the convective flux calculation.
    {
        if (!allow_high_order_interpolation) return;
        if (myConfig.interpolation_order == 1) return;
        if (GlobalConfig.frozen_limiter) return;

        if (cell_list.length == 0) { cell_list = cells; }
        if (vertex_list.length == 0) { vertex_list = vertices; }

        // It is more efficient to determine limiting factor here for some usg limiters.
        final switch (myConfig.unstructured_limiter) {
            case UnstructuredLimiter.svan_albada:
                // do nothing now
                break;
            case UnstructuredLimiter.min_mod:
                // do nothing now
                break;
            case UnstructuredLimiter.barth:
                foreach (c; cell_list) c.gradients.barth_limit(c.cell_cloud, *(c.ws), myConfig);
                break;
            case UnstructuredLimiter.park:
                foreach (c; cell_list) c.gradients.park_limit(c.cell_cloud, *(c.ws), myConfig);
                break;
            case UnstructuredLimiter.hvan_albada:
                foreach (c; cell_list) c.gradients.van_albada_limit(c.cell_cloud, *(c.ws), true, myConfig);
                break;
            case UnstructuredLimiter.van_albada:
                foreach (c; cell_list) c.gradients.van_albada_limit(c.cell_cloud, *(c.ws), false, myConfig);
                break;
            case UnstructuredLimiter.hnishikawa:
                foreach (c; cell_list) c.gradients.nishikawa_limit(c.cell_cloud, *(c.ws), true, myConfig, gtl);
                break;
            case UnstructuredLimiter.nishikawa:
                foreach (c; cell_list) c.gradients.nishikawa_limit(c.cell_cloud, *(c.ws), false, myConfig, gtl);
                break;
            case UnstructuredLimiter.hvenkat_mlp:
                foreach (c; cell_list) c.gradients.venkat_mlp_limit(c.cell_cloud, *(c.ws), true, myConfig, gtl);
                break;
            case UnstructuredLimiter.venkat_mlp:
                foreach (c; cell_list) c.gradients.venkat_mlp_limit(c.cell_cloud, *(c.ws), false, myConfig, gtl);
                break;
            case UnstructuredLimiter.hvenkat:
                foreach (c; cell_list) c.gradients.venkat_limit(c.cell_cloud, *(c.ws), true, myConfig, gtl);
                break;
            case UnstructuredLimiter.venkat:
                number phi_hp = to!number(1.0);
                foreach (c; cell_list) {
                    size_t i = c.id;
                    celldata.lsqgradients[i].venkat_limit2(celldata.flowstates[i],
                        celldata.volumes[i], celldata.face_distances[i],
                        celldata.nfaces[i], phi_hp, myConfig);
                }
                break;
        }
    } // end convective_flux-phase1()

    @nogc
    void copy_cell_data_to_ghost_cells(size_t[] face_idxs)
    {
        // Fill in gradients for ghost cells so that left- and right- cells at all faces,
        // including those along block boundaries, have the latest gradient values.
        foreach(idx; face_idxs){
            if (idx<ninteriorfaces) continue; // skip if given an interior face

            if (facedata.left_interior_only[idx]) {
                size_t l = facedata.f2c[idx].left;
                size_t r = facedata.f2c[idx].right;
                celldata.lsqgradients[r].copy_values_from(celldata.lsqgradients[l]);
            } else if (facedata.right_interior_only[idx]) {
                size_t l = facedata.f2c[idx].left;
                size_t r = facedata.f2c[idx].right;
                celldata.lsqgradients[l].copy_values_from(celldata.lsqgradients[r]);
            } else {
                // Nothing needs to be done for a shared interface,
                continue;
            }
        }
    }

    @nogc
    override void convective_flux_phase2new(bool allow_high_order_interpolation, size_t[] cell_idxs=[], size_t[] face_idxs=[])
    // Make use of the flow gradients to actually do the high-order reconstruction
    // and then compute fluxes of conserved quantities at all faces.
    {
        if (n_ghost_cell_layers > 1) return;

        immutable size_t neq = myConfig.cqi.n;
        immutable bool allow_reconstruction_anywhere = allow_high_order_interpolation && (myConfig.interpolation_order > 1);
        immutable bool suppress_reconstruction_at_boundaries = myConfig.suppress_reconstruction_at_boundaries;
        immutable bool do_all_faces = face_idxs.length==0;

        if (allow_reconstruction_anywhere && !suppress_reconstruction_at_boundaries) {
            if (do_all_faces) {
                copy_cell_data_to_ghost_cells(facedata.all_face_idxs[ninteriorfaces .. nfaces]);
            } else {
                copy_cell_data_to_ghost_cells(face_idxs);
            }
        }

        if (do_all_faces) face_idxs = facedata.all_face_idxs;
        Vector3 gvel;
        gvel.clear();

        // At this point, we should have all gradient values up to date and we are now ready
        // to reconstruct field values and compute the convective fluxes.
        foreach(idx; face_idxs){ // TODO: No ghost cell version
            size_t l = facedata.f2c[idx].left;
            size_t r = facedata.f2c[idx].right;
            Lft.copy_values_from(celldata.flowstates[l]);
            Rght.copy_values_from(celldata.flowstates[r]);
            facedata.flowstates[idx].copy_average_values_from(Lft, Rght);

            // Sort out logic for whether to do recontruction at this face
            bool do_reconstruction = allow_reconstruction_anywhere;
            if (allow_reconstruction_anywhere) {
                bool on_boundary = (idx>=ninteriorfaces);
                if (on_boundary && suppress_reconstruction_at_boundaries){
                    // Check for a real boundary, i.e. one that isn't an exchange boundary
                    if (facedata.left_interior_only[idx] || facedata.right_interior_only[idx]) {
                        do_reconstruction = false;
                    }
                }
            }

            if (do_reconstruction) {
                interp_both(facedata.dL[idx], facedata.dR[idx],
                            celldata.lsqgradients[l], celldata.lsqgradients[r],
                            celldata.flowstates[l], celldata.flowstates[r],
                            myConfig, *Lft, *Rght);
            }
            compute_interface_flux_interior(*Lft, *Rght, facedata.flowstates[idx], myConfig, gvel,
                                            facedata.positions[idx], facedata.normals[idx], facedata.tangents1[idx], facedata.tangents2[idx],
                                            facedata.fluxes[idx*neq .. (idx+1)*neq]);
        }
    } // end convective_flux-phase2()

    @nogc
    override void convective_flux_phase2(bool allow_high_order_interpolation, size_t gtl=0,
                                         FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
    // Make use of the flow gradients to actually do the high-order reconstruction
    // and then compute fluxes of conserved quantities at all faces.
    {
        if (cell_list.length == 0) { cell_list = cells; }
        if (iface_list.length == 0) { iface_list = faces; }

        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            // Fill in gradients for ghost cells so that left- and right- cells at all faces,
            // including those along block boundaries, have the latest gradient values.
            foreach (bcond; bc) {
                if (bcond.ghost_cell_data_available == false) { continue; }
                // Proceed to do some work only if we have ghost cells in which to insert gradients.
                bool found_mapped_cell_bc = false;
                foreach (gce; bcond.preReconAction) {
                    auto mygce = cast(GhostCellMappedCellCopy)gce;
                    if (mygce) {
                        found_mapped_cell_bc = true;
                        // we have already transferred the cell gradients for mapped cells
                    } // end if (mygce)
                } // end foreach gce
                if (!found_mapped_cell_bc) {
                    // There are no other cells backing the ghost cells on this boundary.
                    // Fill in ghost-cell gradients from the other side of the face.
                    foreach (i, f; bcond.faces) {
                        if (bcond.outsigns[i] == 1) {
                            f.right_cell.gradients.copy_values_from(*(f.left_cell.gradients));
                        } else {
                            f.left_cell.gradients.copy_values_from(*(f.right_cell.gradients));
                        }
                    } // end foreach f
                } // end if !found_mapped_cell_bc
            } // end foreach bcond
        } // end if interpolation_order > 1
        //
        // At this point, we should have all gradient values up to date and we are now ready
        // to reconstruct field values and compute the convective fluxes.
        foreach (f; iface_list) {
            bool do_reconstruction = allow_high_order_interpolation && !f.in_suppress_reconstruction_zone;
            if (f.left_cell && f.right_cell) {
                interp_both(myConfig, f, gtl, *Lft, *Rght, do_reconstruction);
                compute_interface_flux_interior(*Lft, *Rght, f, myConfig, omegaz);
            } else if (f.right_cell) {
                interp_right(myConfig, f, gtl, *Rght, do_reconstruction);
                compute_flux_at_left_wall(*Rght, f, myConfig, omegaz);
            } else if (f.left_cell) {
                interp_left(myConfig, f, gtl, *Lft, do_reconstruction);
                compute_flux_at_right_wall(*Lft, f, myConfig, omegaz);
            } else {
                assert(0, "oops, a face without attached cells");
            }
        } // end foreach face
    } // end convective_flux-phase2()

    override size_t[] get_cell_write_indices() {
        size_t[] index;
        index = iota(0, ncells, myConfig.nic_write).array();
        return index;
    }

    override void eval_udf_source_vectors(double simTime, size_t[] cell_idxs=[])
    {
        if (myConfig.udf_source_terms) {
            if (cell_idxs.length==0) cell_idxs = celldata.all_cell_idxs;
            foreach (i; cell_idxs) {
                auto cell = cells[i];
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                getUDFSourceTermsForCell(myL, cell, 0, simTime, myConfig, id, i_cell, j_cell, k_cell);
                cell.add_udf_source_vector();
            }
        }
    }
} // end class UFluidBlock
