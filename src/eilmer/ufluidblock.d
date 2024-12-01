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
import util.json_helper;
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
import block;

class UFluidBlock: FluidBlock {
public:
    size_t ncells;
    size_t nvertices;
    size_t nfaces;
    size_t nboundaries;
    UnstructuredGrid grid;
    size_t ninteriorfaces;
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.

public:
    this(in int id, JSONValue json_data)
    {
        label = getJSONstring(json_data, "label", "");
        ncells = getJSONint(json_data, "ncells", 0);
        // For an unstructured-grid, n_ghost_cell_layers=1
        super(id, Grid_t.unstructured_grid, ncells, 1, label);
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
            cells[cell_idx] = new FVCell(myConfig, pos, *myfs, volume, to!int(cell_idx));
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
        uint nsp = myConfig.n_species;
        foreach (i, v; grid.vertices) {
            auto new_vtx = new FVVertex(myConfig, lsq_workspace_at_vertices, to!int(i));
            new_vtx.pos[0] = v;
            vertices ~= new_vtx;
        }
        // sync_vertices_from_underlying_grid(0); // redundant, if done just above
        bool lsq_workspace_at_faces = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces);
        foreach (i, f; grid.faces) {
            auto new_face = new FVInterface(myConfig, IndexDirection.none, lsq_workspace_at_faces, to!int(i));
            faces ~= new_face;
        }
        foreach (i, c; grid.cells) {
            // Note that the cell id and the index in the cells array are the same.
            // We will reply upon this connection in other parts of the flow code.
            bool lsq_workspace_at_cells = (myConfig.viscous) && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
                && (myConfig.spatial_deriv_locn == SpatialDerivLocn.cells);
            auto new_cell = new FVCell(myConfig, lsq_workspace_at_cells, to!int(i));
            new_cell.contains_flow_data = true;
            new_cell.is_interior_to_domain = true;
            cells ~= new_cell;
        }
        // Bind the interfaces, vertices and cells together,
        // using the indices stored in the unstructured grid.
        foreach (i, f; faces) {
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
                c.outsign ~= my_outsign;
                if (my_outsign == 1) {
                    if (my_face.left_cell) {
                        string msg = format("Already have cell %d attached to left-of-face %d. Attempt to add cell %d.",
                                            my_face.left_cell.id, my_face.id, c.id);
                        throw new FlowSolverException(msg);
                    } else {
                        my_face.left_cell = c;
                    }
                } else {
                    if (my_face.right_cell) {
                        string msg = format("Already have cell %d attached to right-of-face %d. Attempt to add cell %d.",
                                            my_face.right_cell.id, my_face.id, c.id);
                        throw new FlowSolverException(msg);
                    } else {
                        my_face.right_cell = c;
                    }
                }
            }
        } // end foreach cells
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
        int ghost_cell_count = 0;
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
                    // Make ghost-cell id values distinct from FVCell ids so that
                    // the warning/error messages are somewhat informative.
                    FVCell ghost0 = new FVCell(myConfig, false, ghost_cell_start_id+ghost_cell_count);
                    ghost_cell_count++;
                    ghost0.contains_flow_data = bc[i].ghost_cell_data_available;
                    bc[i].ghostcells ~= ghost0;
                    if (my_outsign == 1) {
                        if (my_face.right_cell) {
                            string msg = format("Already have cell %d attached to right-of-face %d."
                                                ~" Attempt to add ghost cell %d.",
                                                my_face.right_cell.id, my_face.id, ghost0.id);
                            throw new FlowSolverException(msg);
                        } else {
                            my_face.right_cell = ghost0;
                        }
                    } else {
                        if (my_face.left_cell) {
                            string msg = format("Already have cell %d attached to left-of-face %d."
                                                ~" Attempt to add ghost cell %d.",
                                                my_face.left_cell.id, my_face.id, ghost0.id);
                            throw new FlowSolverException(msg);
                        } else {
                            my_face.left_cell = ghost0;
                        }
                    }
                } // end if (bc[i].ghost_cell_data_available
            } // end foreach j
        } // end foreach i
        // At this point, all faces should have either one finite-volume cell
        // or one ghost cell attached to each side -- check that this is true.
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
        auto nmodes = myConfig.n_modes;
        auto nturb = myConfig.turb_model.nturb;
        if (myConfig.use_extended_stencil) {
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
                        }
                    } else {
                        if (f.left_cell && f.left_cell.contains_flow_data) {
                            c.cell_cloud ~= f.left_cell;
                            cell_ids ~= f.left_cell.id;
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
                            }
                        } else {
                            if (f.left_cell && cell_ids.canFind(f.left_cell.id) == false && f.left_cell.contains_flow_data) {
                                c.cell_cloud ~= f.left_cell;
                                cell_ids ~= f.left_cell.id;
                            }
                        }
                    } // end foreach face
                } else {
                    // apply nearest-node neighbour
                    foreach(vtx; c.vtx) {
                        foreach(cid; cellIndexListPerVertex[vtx.id])
                            if (cell_ids.canFind(cid) == false && cells[cid].contains_flow_data)
                                { c.cell_cloud ~= cells[cid]; cell_ids ~= cid; }
                    }
                }
                c.ws = new LSQInterpWorkspace();
                c.gradients = new LSQInterpGradients(nsp, nmodes, nturb);
            } // end foreach cell
        } else {
            foreach (c; cells) {
                // First cell in the cloud is the cell itself.  Differences are taken about it.
                c.cell_cloud ~= c;
                // Subsequent cells are the surrounding cells.
                foreach (i, f; c.iface) {
                    if (c.outsign[i] == 1) {
                        if (f.right_cell && f.right_cell.contains_flow_data) {
                            c.cell_cloud ~= f.right_cell;
                        }
                    } else {
                        if (f.left_cell && f.left_cell.contains_flow_data) {
                            c.cell_cloud ~= f.left_cell;
                        }
                    }
                } // end foreach face
                c.ws = new LSQInterpWorkspace();
                c.gradients = new LSQInterpGradients(nsp, nmodes, nturb);
            } // end foreach cell
        } // end else
        // We will also need derivative storage in ghostcells because the
        // reconstruction functions will expect to be able to access the gradients
        // either side of each interface.
        // We will be able to get these gradients from the mapped-cells
        // in an adjoining block.
        // [TODO] think about this for the junction of usgrid and sgrid blocks.
        // The sgrid blocks will not have the gradients stored within the cells.
        foreach (bci; bc) {
            if (bci.ghost_cell_data_available) {
                foreach (c; bci.ghostcells) {
                    c.gradients = new LSQInterpGradients(nsp, nmodes, nturb);
                }
            }
        }
        //
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
                    c.cloud_fs ~= &(c.fs);
                    // Subsequent cells are the surrounding cells.
                    foreach (i, f; c.iface) {
                        c.cloud_pos ~= &(f.pos);
                        c.cloud_fs ~= &(f.fs);
                    } // end foreach face
                }
                // We will also need derivative storage in ghostcells because the
                // reconstruction functions will expect to be able to access the gradients
                // either side of each interface.
                // We will be able to get these gradients from the mapped-cells
                // in an adjoining block.
                foreach (bci; bc) {
                    if (bci.ghost_cell_data_available) {
                        foreach (c; bci.ghostcells) {
                            c.grad = new FlowGradients(myConfig);
                        }
                    }
                }
            } // end switch (myConfig.spatial_deriv_locn)
        } // if myConfig.viscous
    } // end init_grid_and_flow_arrays()

    void build_cloud_of_cell_references_at_each_vertex()
    {
        // For the MLP limiter, we need access to the gradients stored in the cells.
        foreach (vtx; vertices) {
            foreach(cid; cellIndexListPerVertex[vtx.id]) { vtx.cell_cloud ~= cells[cid]; }
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
                    cell.grad.set_up_workspace_leastsq(cell.cloud_pos, cell.pos[gtl], false, *(cell.ws_grad));
                }
            } else if (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces) {
                foreach(iface; faces) {
                    iface.grad.set_up_workspace_leastsq(iface.cloud_pos, iface.pos, false, *(iface.ws_grad));
                }
            } else { // myConfig.spatial_deriv_locn == vertices
                foreach(vtx; vertices) {
                    vtx.grad.set_up_workspace_leastsq(vtx.cloud_pos, vtx.pos[gtl], true, *(vtx.ws_grad));
                }
            }
        }
        // The LSQ linear model for the flow field reconstruction is fitted using information
        // on the locations of the points.
        // This model is used later, as part of the convective flux calculation.
        if (myConfig.interpolation_order > 1) {
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

    @nogc
    override void average_turbulent_transprops_to_faces()
    {
        if (!myConfig.turb_model.isTurbulent) return;

        foreach(idx; 0 .. ninteriorfaces){
            faces[idx].fs.mu_t = 0.5*(faces[idx].left_cell.fs.mu_t + faces[idx].right_cell.fs.mu_t);
            faces[idx].fs.k_t = 0.5*(faces[idx].left_cell.fs.k_t + faces[idx].right_cell.fs.k_t);
        }

        foreach(idx; ninteriorfaces .. nfaces){
            faces[idx].average_turbulent_transprops();
        }
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
    override void set_face_flowstates_to_averages_from_cells()
    {
        // It turns out that some shock-detectors need flow derivatives before the
        // convective-flux calculation is done.  That used to be the only place
        // that the face FlowState was filled in and it was done as a side-effect,
        // which has confused just about everyone at some time in their work on the code.
        foreach (f; faces) {
            if (f.left_cell && f.right_cell) {
                f.fs.copy_average_values_from(f.left_cell.fs, f.right_cell.fs);
            } else if (f.right_cell) {
                f.fs.copy_values_from(f.right_cell.fs);
            } else if (f.left_cell) {
                f.fs.copy_values_from(f.left_cell.fs);
            } else {
                assert(0, "oops, a face without attached cells");
            }
        }
    }

    @nogc
    override void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0,
                                         FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
    // Compute gradients of flow quantities for higher-order reconstruction, if required.
    // To be used, later, in the convective flux calculation.
    {

        if (cell_list.length == 0) { cell_list = cells; }
        if (vertex_list.length == 0) { vertex_list = vertices; }

        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            foreach (c; cell_list) {
                if (GlobalConfig.frozen_limiter == false) {
                    if (myConfig.unstructured_limiter == UnstructuredLimiter.venkat_mlp) {
                        c.gradients.store_max_min_values_for_extended_stencil(c.cell_cloud, myConfig);
                    } else {
                        c.gradients.store_max_min_values_for_compact_stencil(c.cell_cloud, myConfig);
                    }
                }
                c.gradients.compute_lsq_values(c.cell_cloud, *(c.ws), myConfig);
            }
        } // end if interpolation_order > 1
    } // end convective_flux-phase0()

        @nogc
    override void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0,
                                         FVCell[] cell_list = [], FVInterface[] iface_list = [], FVVertex[] vertex_list = [])
        // Compute limiter values of flow quantities for higher-order reconstruction, if required.
        // To be used, later, in the convective flux calculation.
    {

        if (cell_list.length == 0) { cell_list = cells; }
        if (vertex_list.length == 0) { vertex_list = vertices; }

        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            if (GlobalConfig.frozen_limiter == false) {
                foreach (c; cell_list) {
                    // It is more efficient to determine limiting factor here for some usg limiters.
                    final switch (myConfig.unstructured_limiter) {
                        case UnstructuredLimiter.svan_albada:
                            // do nothing now
                            break;
                        case UnstructuredLimiter.min_mod:
                            // do nothing now
                            break;
                        case UnstructuredLimiter.barth:
                            c.gradients.barth_limit(c.cell_cloud, *(c.ws), myConfig);
                            break;
                        case UnstructuredLimiter.park:
                            c.gradients.park_limit(c.cell_cloud, *(c.ws), myConfig);
                            break;
                        case UnstructuredLimiter.hvan_albada:
                            c.gradients.van_albada_limit(c.cell_cloud, *(c.ws), true, myConfig);
                            break;
                        case UnstructuredLimiter.van_albada:
                            c.gradients.van_albada_limit(c.cell_cloud, *(c.ws), false, myConfig);
                            break;
                        case UnstructuredLimiter.hnishikawa:
                            c.gradients.nishikawa_limit(c.cell_cloud, *(c.ws), true, myConfig, gtl);
                            break;
                        case UnstructuredLimiter.nishikawa:
                            c.gradients.nishikawa_limit(c.cell_cloud, *(c.ws), false, myConfig, gtl);
                            break;
                        case UnstructuredLimiter.hvenkat_mlp:
                            c.gradients.venkat_mlp_limit(c.cell_cloud, *(c.ws), true, myConfig, gtl);
                            break;
                        case UnstructuredLimiter.venkat_mlp:
                            c.gradients.venkat_mlp_limit(c.cell_cloud, *(c.ws), false, myConfig, gtl);
                            break;
                        case UnstructuredLimiter.hvenkat:
                            c.gradients.venkat_limit(c.cell_cloud, *(c.ws), true, myConfig, gtl);
                            break;
                        case UnstructuredLimiter.venkat:
                            c.gradients.venkat_limit(c.cell_cloud, *(c.ws), false, myConfig, gtl);
                            break;
                    } // end switch
                } // end foreach c
            } // if (frozen_limter == false)
        } // end if interpolation_order > 1
    } // end convective_flux-phase0()

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
    } // end convective_flux-phase1()

    override size_t[] get_cell_write_indices() {
        size_t[] index;
        index = iota(0, ncells, myConfig.nic_write).array();
        return index;
    }
} // end class UFluidBlock
