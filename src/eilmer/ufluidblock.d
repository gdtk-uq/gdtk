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
import lsqinterp;
import fluidblock;
import bc;


class UFluidBlock: FluidBlock {
public:
    size_t ncells;
    size_t nvertices;
    size_t nfaces;
    size_t nboundaries;
    UnstructuredGrid grid;
    // Work-space that gets reused.
    // The following objects are used in the convective_flux method.
    LsqInterpolator lsq;

public:
    this(in int id, JSONValue json_data)
    {
        label = getJSONstring(json_data, "label", "");
        ncells = getJSONint(json_data, "ncells", 0);
        super(id, Grid_t.unstructured_grid, ncells, label);
        nvertices = getJSONint(json_data, "nvertices", 0);
        nfaces = getJSONint(json_data, "nfaces", 0);
        nboundaries = getJSONint(json_data, "nboundaries", 0);
        active = getJSONbool(json_data, "active", true);
        omegaz = getJSONdouble(json_data, "omegaz", 0.0);
        // Workspace for flux_calc method.
        lsq = new LsqInterpolator(dedicatedConfig[id]);
    } // end constructor from json

    @nogc override int get_interpolation_order()
    {
        return lsq.get_interpolation_order();
    }

    @nogc override void set_interpolation_order(int order)
    {
        lsq.set_interpolation_order(order);
    }

    override void init_lua_globals()
    {
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
    } // end init_lua_globals()

    override void init_boundary_conditions(JSONValue json_data)
    // Initialize boundary conditions after the blocks are fully constructed,
    // because we want access to the full collection of valid block references.
    {
        foreach (boundary; 0 .. nboundaries) {
            string json_key = format("boundary_%d", boundary);
            auto bc_json_data = json_data[json_key];
            bc ~= make_BC_from_json(bc_json_data, id, to!int(boundary));
        }
        foreach (bci; bc) bci.post_bc_construction();
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

    // The following 5 access methods are here to match the structured-grid API
    // but they're really not intended for serious use on the unstructured-grid.
    @nogc 
    override ref FVCell get_cell(size_t i, size_t j, size_t k=0) 
    {
        return cells[i]; // j, k ignored
    }
    @nogc 
    override ref FVInterface get_ifi(size_t i, size_t j, size_t k=0) 
    {
        return faces[i];
    }
    @nogc
    override ref FVInterface get_ifj(size_t i, size_t j, size_t k=0)
    {
        return faces[i];
    }
    @nogc
    override ref FVInterface get_ifk(size_t i, size_t j, size_t k=0)
    {
        return faces[i];
    }
    @nogc
    override ref FVVertex get_vtx(size_t i, size_t j, size_t k=0)
    {
        return vertices[i];
    }

    override void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
        grid.find_enclosing_cell(p, indx, found); // delegate to the grid object
    }

    override void init_grid_and_flow_arrays(string gridFileName)
    {
        grid = new UnstructuredGrid(gridFileName, myConfig.grid_format);
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
        // Assemble array storage for finite-volume cells, etc.
        bool lsq_workspace_at_vertices = (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.vertices);
        foreach (i, v; grid.vertices) {
            auto new_vtx = new FVVertex(myConfig, lsq_workspace_at_vertices, to!int(i));
            if (myConfig.unstructured_limiter == UnstructuredLimiter.mlp)
                new_vtx.gradients = new LSQInterpGradients(myConfig.gmodel.n_species, myConfig.gmodel.n_modes);
            new_vtx.pos[0] = v;
            vertices ~= new_vtx;
        }
        // sync_vertices_from_underlying_grid(0); // redundant, if done just above
        bool lsq_workspace_at_faces = (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
            && (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces);
        foreach (i, f; grid.faces) {
            auto new_face = new FVInterface(myConfig, lsq_workspace_at_faces, to!int(i));
            faces ~= new_face;
        }
        foreach (i, c; grid.cells) {
            // Note that the cell id and the index in the cells array are the same.
            // We will reply upon this connection in other parts of the flow code.
            auto new_cell = new FVCell(myConfig, to!int(i));
            new_cell.will_have_valid_flow = true;
            new_cell.is_interior = true;
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
                auto my_outsign = grid.cells[i].outsign_list[j];
                c.iface ~= my_face;
                c.outsign ~= to!double(my_outsign);
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
                if (bc[i].ghost_cell_data_available) {
                    // Make ghost-cell id values distinct from FVCell ids so that
                    // the warning/error messages are somewhat informative. 
                    FVCell ghost0 = new FVCell(myConfig, ghost_cell_start_id+ghost_cell_count);
                    ghost_cell_count++;
                    ghost0.will_have_valid_flow = bc[i].ghost_cell_data_available;
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
        auto nsp = myConfig.gmodel.n_species;
        auto nmodes = myConfig.gmodel.n_modes;
        if (myConfig.use_extended_stencil) {
            foreach (c; cells) {
                // First cell in the cloud is the cell itself.  Differences are taken about it.
                c.cell_cloud ~= c;
                size_t[] cell_ids;
                cell_ids ~= c.id;
                // Subsequent cells are the surrounding cells.
                bool is_on_boundary = false;
                foreach(face; c.iface) if (face.is_on_boundary) is_on_boundary = true;
                if (is_on_boundary) {
                    // apply nearest-face neighbour
                    foreach (i, f; c.iface) {
                        if (c.outsign[i] > 0.0) {
                            if (f.right_cell && f.right_cell.will_have_valid_flow) {
                                c.cell_cloud ~= f.right_cell;
                            }
                        } else {
                            if (f.left_cell && f.left_cell.will_have_valid_flow) {
                                c.cell_cloud ~= f.left_cell;
                            }
                        }
                    } // end foreach face
                } else {
                    // apply nearest-node neighbour
                    foreach(vtx; c.vtx) {
                        foreach(cid; cellIndexListPerVertex[vtx.id])
                            if (cell_ids.canFind(cid) == false && cells[cid].will_have_valid_flow)
                                { c.cell_cloud ~= cells[cid]; cell_ids ~= cid; } 
                    }
                } 
                c.ws = new LSQInterpWorkspace();
                c.gradients = new LSQInterpGradients(nsp, nmodes);
            } // end foreach cell
        } else {
            foreach (c; cells) {
                // First cell in the cloud is the cell itself.  Differences are taken about it.
                c.cell_cloud ~= c;
                // Subsequent cells are the surrounding cells.
                foreach (i, f; c.iface) {
                    if (c.outsign[i] > 0.0) {
                        if (f.right_cell && f.right_cell.will_have_valid_flow) {
                            c.cell_cloud ~= f.right_cell;
                        }
                    } else {
                        if (f.left_cell && f.left_cell.will_have_valid_flow) {
                            c.cell_cloud ~= f.left_cell;
                        }
                    }
                } // end foreach face
                c.ws = new LSQInterpWorkspace();
                c.gradients = new LSQInterpGradients(nsp, nmodes);
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
                    c.gradients = new LSQInterpGradients(nsp, nmodes);
                }
            }
        }
        //
        // We will now store the cloud of points in cloud_pos for viscous derivative calcualtions.
        // This is equivalent to store_references_for_derivative_calc(size_t gtl) in sblock.d
        // first we define some functions to make the code cleaner (assume gtl=0)
        void boundary_face_cloud_wall(FVInterface f, size_t i, BoundaryCondition bc, size_t bndary_idx, BoundaryFaceSet boundary) {
            FVCell cell;
            // store the interface
            f.cloud_pos ~= &(f.pos);
            f.cloud_fs ~= f.fs;
            if (bc.outsigns[i] == 1) {
                // store internal neighbour cell
                f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
                f.cloud_fs ~= f.left_cell.fs;
                cell = f.left_cell;
            } else {
                // store internal neighbour cell
                f.cloud_pos ~= &(f.right_cell.pos[0]); // assume gtl = 0
                f.cloud_fs ~= f.right_cell.fs;
                cell = f.right_cell;
            }
            // now grab the remaining cloud points i.e. the faces of the neighbour cell
            foreach (other_face; cell.iface) { // loop around cell interfaces
                if(other_face.id != f.id) { // skip boundary interface
                    // store interface
                    f.cloud_pos ~= &(other_face.pos);
                    f.cloud_fs ~= other_face.fs;
                }
            }
        } // end void boundary_face_cloud_wall()
        void boundary_face_cloud_generic(FVInterface f, size_t i, BoundaryCondition bc, size_t bndary_idx, BoundaryFaceSet boundary){
            double[] cell_cloud_face_ids;
            double[] cell_cloud_cell_ids;
            FVCell[] cell_list;
            // store interface
            f.cloud_pos ~= &(f.pos);
            f.cloud_fs ~= f.fs;
            cell_cloud_face_ids ~= f.id;
            if (bc.outsigns[i] == 1) {
                // left cell should be interior
                // store as a neighbour cell
                f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
                f.cloud_fs ~= f.left_cell.fs;
                cell_list~= f.left_cell.cell_cloud;
                cell_cloud_cell_ids ~= f.left_cell.id;
                // store ghost0
                if (f.right_cell && f.right_cell.will_have_valid_flow) {
                    f.cloud_pos ~= &(f.right_cell.pos[0]);
                    f.cloud_fs ~= f.right_cell.fs;
                    cell_cloud_cell_ids ~= f.right_cell.id;
                }
            } else {
                // right_cell should be interior
                // store as a neighbour cell
                f.cloud_pos ~= &(f.right_cell.pos[0]); // assume gtl = 0
                f.cloud_fs ~= f.right_cell.fs;
                cell_list ~= f.right_cell.cell_cloud;
                cell_cloud_cell_ids ~= f.right_cell.id;
                // store ghost0
                if (f.left_cell && f.left_cell.will_have_valid_flow) { 
                    f.cloud_pos ~= &(f.left_cell.pos[0]);
                    f.cloud_fs ~= f.left_cell.fs;
                    cell_cloud_cell_ids ~= f.left_cell.id;
                }
            }
            foreach (fvtx; f.vtx) {
                foreach (other_face_id; faceIndexListPerVertex[fvtx.id]) {
                    FVInterface other_face = faces[other_face_id];
                    if (other_face.is_on_boundary &&
                        other_face.bc_id == bndary_idx &&
                        cell_cloud_face_ids.canFind(other_face.id) == false) {
                        // store left and right cell
                        if (other_face.left_cell && cell_cloud_cell_ids.canFind(other_face.left_cell.id) == false) {
                            f.cloud_pos ~= &(other_face.left_cell.pos[0]); // assume gtl = 0
                            f.cloud_fs ~= other_face.left_cell.fs;
                            cell_cloud_cell_ids ~= other_face.left_cell.id;
                        }
                        if (other_face.right_cell && cell_cloud_cell_ids.canFind(other_face.right_cell.id) == false) {
                            f.cloud_pos ~= &(other_face.right_cell.pos[0]); // assume gtl = 0
                            f.cloud_fs ~= other_face.right_cell.fs;
                            cell_cloud_cell_ids ~= other_face.right_cell.id;
                        }
                        cell_cloud_face_ids ~= other_face.id;
                    }
                }
            } // end foreach (fvtx; f.vtx)
            // in some instances we might have a degenerate cloud, let's make sure we have enough points
            if (f.cloud_pos.length <= 3) {
                foreach (fvtx; f.vtx) {
                    foreach (other_face_id; faceIndexListPerVertex[fvtx.id]) {
                        FVInterface other_face = faces[other_face_id];
                        if (other_face.is_on_boundary && cell_cloud_face_ids.canFind(other_face.id) == false) {
                            // store left and right cell
                            if (other_face.left_cell && cell_cloud_cell_ids.canFind(other_face.left_cell.id) == false) {
                                f.cloud_pos ~= &(other_face.left_cell.pos[0]); // assume gtl = 0
                                f.cloud_fs ~= other_face.left_cell.fs;
                                cell_cloud_cell_ids ~= other_face.left_cell.id;
                            }
                            if (other_face.right_cell && cell_cloud_cell_ids.canFind(other_face.right_cell.id) == false) {
                                f.cloud_pos ~= &(other_face.right_cell.pos[0]); // assume gtl = 0
                                f.cloud_fs ~= other_face.right_cell.fs;
                                cell_cloud_cell_ids ~= other_face.right_cell.id;
                            }
                            cell_cloud_face_ids ~= other_face.id;
                        }
                    }
                } // end foreach (fvtx; f.vtx)
            } // end if (f.cloud_pos.length <= 3) {
        } // end void boundary_face_cloud_generic()
        void internal_face_cloud_interfaces(FVInterface f, size_t i) {
            if ( f.is_on_boundary) { return;} // boundaries already set
            // store interface
            f.cloud_pos ~= &(f.pos);
            f.cloud_fs ~= f.fs;
            FVCell[] cell_list;
            cell_list ~= f.left_cell;
            cell_list ~= f.right_cell;
            foreach (cell; cell_list) {
                // store cell
                f.cloud_pos ~= &(cell.pos[0]); // assume gtl = 0
                f.cloud_fs ~= cell.fs;
                // now grab the cell interfaces
                foreach (other_face; cell.iface) {
                    if (other_face.id == f.id) continue;
                    // store interface
                    f.cloud_pos ~= &(other_face.pos);
                    f.cloud_fs ~= other_face.fs;
                } // end foreach (other_face; cell.iface)
            } // end foreach (cell; cell_list)
        } // end void internal_face_cloud_interfaces()
        void internal_face_cloud_cells(FVInterface f, size_t i) {
            if ( f.is_on_boundary) { return;} // boundaries already set
            double[] cell_cloud_face_ids;
            double[] cell_cloud_cell_ids;
            cell_cloud_face_ids ~= f.id;
            // store interface
            f.cloud_pos ~= &(f.pos);
            f.cloud_fs ~= f.fs;
            cell_cloud_face_ids ~= f.id;
            // store neighbour cell
            f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
            f.cloud_fs ~= f.left_cell.fs;
            cell_cloud_cell_ids ~= f.left_cell.id;
            // store neighbour cell
            f.cloud_pos ~= &(f.right_cell.pos[0]);
            f.cloud_fs ~= f.right_cell.fs;
            cell_cloud_cell_ids ~= f.right_cell.id;
            foreach (fvtx; f.vtx) {
                foreach (other_face_id; faceIndexListPerVertex[fvtx.id]) {
                    FVInterface other_face = faces[other_face_id];
                    if (cell_cloud_face_ids.canFind(other_face.id) == false) {
                        // store left and right cell
                        if (cell_cloud_cell_ids.canFind(other_face.left_cell.id) == false &&
                            other_face.left_cell.id < ghost_cell_start_id) {
                            f.cloud_pos ~= &(other_face.left_cell.pos[0]); // assume gtl = 0
                            f.cloud_fs ~= other_face.left_cell.fs;
                            cell_cloud_cell_ids ~= other_face.left_cell.id;
                        }
                        if (cell_cloud_cell_ids.canFind(other_face.right_cell.id) == false &&
                            other_face.right_cell.id < ghost_cell_start_id) {
                            f.cloud_pos ~= &(other_face.right_cell.pos[0]); // assume gtl = 0
                            f.cloud_fs ~= other_face.right_cell.fs;
                            cell_cloud_cell_ids ~= other_face.right_cell.id;
                        }
                        cell_cloud_face_ids ~= other_face.id;
                    }
                }
            } // end foreach (fvtx; f.vtx)
        } // end void internal_face_cloud_cells()
        if (myConfig.spatial_deriv_calc ==  SpatialDerivCalc.divergence) {
            throw new Error("Divergence theorem not implemented for unstructured grid");
        }
        // else continue on to fill least-squares/finite-difference cloud
        final switch (myConfig.spatial_deriv_locn) {
        case SpatialDerivLocn.vertices:
            throw new Error("spatial_deriv_locn at vertices not implemented for unstructured grid");
            // no need for break;
        case SpatialDerivLocn.faces:
            if (myConfig.spatial_deriv_calc ==  SpatialDerivCalc.least_squares) {  
                // this cloud allows for consistent gradients at block
                // boundaries, necessary to have 2nd order convergence for multiblock simulations.
                // 1. BOUNDARY INTERFACE CLOUDS
                foreach(bndary_idx, boundary; grid.boundaries) { // set boundary clouds first
                    BoundaryCondition bc = this.bc[bndary_idx];
                    // TO_DO: currently we have to check whether the boundary interface is on a wall,
                    //        using wall ghost cell data is causing issues for the viscous gradient calculations
                    //        at these interfaces. We are working towards removing ghost cell data from
                    //        wall type boundaries. In the mean time we will apply a temporary fix by using the
                    //        point cloud stencil which does not use any ghost cell data at walls. KD 17/06/2016
                    if (bc.is_wall_with_viscous_effects ||
                        myConfig.include_ghost_cells_in_spatial_deriv_clouds == false) {
                        // If the bc type is WALL OR we're not using ghost cell data then
                        // only use internal data 
                        foreach (i, f; bc.faces) { 
                            boundary_face_cloud_wall(f, i, bc, bndary_idx, boundary);
                        } // end foreach (i, f; bc.faces)
                    } else {
                        foreach (i, f; bc.faces) {
                            boundary_face_cloud_generic(f, i, bc, bndary_idx, boundary);
                        }
                    }
                } // end (bndary_idx, boundary; grid.boundaries)
                // 2. INTERNAL INTERFACE CLOUDS
                foreach (i, f; faces) {
                    internal_face_cloud_interfaces(f, i);
                }
            } // end if (myConfig.include_ghost_cells_in_spatial_deriv_clouds)
        } // end switch (myConfig.spatial_deriv_locn)
    } // end init_grid_and_flow_arrays()
    
    override void compute_primary_cell_geometric_data(size_t gtl)
    {
        if (myConfig.dimensions == 2) {
            foreach (c; cells) { c.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
            foreach (f; faces) { f.update_2D_geometric_data(gtl, myConfig.axisymmetric); }
        } else { // 3D
            foreach (c; cells) { c.update_3D_geometric_data(gtl); }
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
                    Vector3 delta = my_face.pos - inside0.pos[gtl];
                    auto ghost0 = my_face.right_cell;
                    ghost0.pos[gtl] = my_face.pos + delta;
                    ghost0.iLength = inside0.iLength;
                    ghost0.jLength = inside0.jLength;
                    ghost0.kLength = inside0.kLength;
                    ghost0.L_min = inside0.L_min;
                } else {
                    auto inside0 = my_face.right_cell;
                    Vector3 delta = my_face.pos - inside0.pos[gtl];
                    auto ghost0 = my_face.left_cell;
                    ghost0.pos[gtl] = my_face.pos + delta;
                    ghost0.iLength = inside0.iLength;
                    ghost0.jLength = inside0.jLength;
                    ghost0.kLength = inside0.kLength;
                    ghost0.L_min = inside0.L_min;
                } // end if my_outsign
            } // end foreach j
        } // end foreach bndry
    } // end compute_primary_cell_geometric_data()

    override void compute_least_squares_setup(size_t gtl)
    {
        // Update the least-squares geometric weights and the workspaces, if appropriate.
        // The weights should be calculated when the grid is initialised or moved.
        //
        if (myConfig.viscous && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)) {
            // LSQ weights are used in the calculation of flow gradients for the viscous terms.
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
        // The LSQ linear model for the flow field reconstruction is fitted using information
        // on the locations of the points.
        // This model is used later, as part of the convective flux calculation.
        if (myConfig.interpolation_order > 1) {
            foreach (c; cells) {
                try {
                    c.ws.assemble_and_invert_normal_matrix(c.cell_cloud, myConfig.dimensions, gtl);
                } catch (Exception e) {
                    writefln("In compute_least_squares_setup()," ~
                             " we have failed to assemble and invert normal matrix for cell id=%d",
                             c.id);
                    throw e;
                }
            }
        }
    } // end compute_least_squares_setup()

    override void sync_vertices_from_underlying_grid(size_t gtl=0)
    {
        foreach (i; 0 .. vertices.length) { vertices[i].pos[gtl].set(grid[i]); }
    }
    
    override void sync_vertices_to_underlying_grid(size_t gtl=0)
    {
        foreach (i; 0 .. vertices.length) { grid[i].set(vertices[i].pos[gtl]); }
    }

    override void read_new_underlying_grid(string fileName)
    {
        if (myConfig.verbosity_level > 1) { writeln("read_new_underlying_grid() for block ", id); }
        grid = new UnstructuredGrid(fileName, myConfig.grid_format);
        grid.sort_cells_into_bins();
    }
    
    override void write_underlying_grid(string fileName)
    {
        if (myConfig.verbosity_level > 1) { writeln("write_underlying_grid() for block ", id); }
        grid.write(fileName, myConfig.grid_format);
    }

    override double read_solution(string filename, bool overwrite_geometry_data)
    // Note that this function needs to be kept in sync with the BlockFlow class
    // over in flowsolution.d and with write_solution() below and with
    // write_initial_usg_flow_file_from_lua() in luaflowstate.d. 
    // Returns sim_time from file.
    {
        if (myConfig.verbosity_level > 1) { writeln("read_solution(): Start block ", id); }
        double sim_time; // to be read from file
        string[] expected_variable_list = myConfig.flow_variable_list.dup();
        string myLabel;
        size_t nvariables;
        int my_dimensions;
        size_t nc;
        switch (myConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File fin = File(filename, "rb");
            string expected_header = "unstructured_grid_flow 1.0";
            char[] found_header = new char[expected_header.length];
            fin.rawRead(found_header);
            if (found_header != expected_header) {
                throw new FlowSolverException("UFluidBlock.read_solution from raw_binary_file: " ~
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
            int[2] int2; fin.rawRead(int2); my_dimensions = int2[0]; nc = int2[1];
            if (my_dimensions != myConfig.dimensions) {
                string msg = text("dimensions found: " ~ to!string(my_dimensions));
                throw new FlowSolverException(msg);
            }
            if (nc != ncells) {
                string msg = text("For block[", id, "] we have a mismatch in solution size.",
                                  " Have read nc=", nc, " ncells=", ncells);
                throw new FlowSolverException(msg);
            }   
            foreach (i; 0 .. ncells) {
                cells[i].read_values_from_raw_binary(fin, overwrite_geometry_data);
            }
            break;
        default:
            auto byLine = new GzipByLine(filename);
            auto line = byLine.front; byLine.popFront();
            string format_version;
            formattedRead(line, "unstructured_grid_flow %s", &format_version);
            if (format_version != "1.0") {
                string msg = text("UFluidBlock.read_solution(): " ~ "format version found: "
                                  ~ format_version);
                throw new FlowSolverException(msg); 
            }
            line = byLine.front; byLine.popFront();
            formattedRead(line, "label: %s", &myLabel);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "sim_time: %g", &sim_time);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "variables: %d", &nvariables);
            line = byLine.front; byLine.popFront();
            // ingore variableNames = line.strip().split();
            line = byLine.front; byLine.popFront();
            formattedRead(line, "dimensions: %d", &my_dimensions);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "ncells: %d", &nc);
            if (nc != ncells) {
                string msg = text("For block[", id, "] we have a mismatch in solution size.",
                                  " Have read nc=", nc, " ncells=", ncells);
                throw new FlowSolverException(msg);
            }   
            foreach (i; 0 .. ncells) {
                line = byLine.front; byLine.popFront();
                cells[i].scan_values_from_string(line, overwrite_geometry_data);
            }
        } // end switch flow_format
        return sim_time;
    } // end read_solution()

    override void write_solution(string filename, double sim_time)
    // Write the flow solution (i.e. the primary variables at the cell centers)
    // for a single block.
    // Keep this function in sync with
    // write_initial_usg_flow_file_from_lua() from luaflowstate.d and
    // write_initial_flow_file() from flowstate.d.
    {
        if (myConfig.verbosity_level > 1) { writeln("write_solution(): Start block ", id); }
        switch (myConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File outfile = File(filename, "wb");
            int[1] int1; int[2] int2; double[1] dbl1; // buffer arrays
            string header = "unstructured_grid_flow 1.0";
            outfile.rawWrite(to!(char[])(header));
            int1[0] = to!int(label.length); outfile.rawWrite(int1);
            if (label.length > 0) { outfile.rawWrite(to!(char[])(label)); }
            dbl1[0] = sim_time; outfile.rawWrite(dbl1);
            int1[0] = to!int(myConfig.flow_variable_list.length); outfile.rawWrite(int1);
            foreach(varname; myConfig.flow_variable_list) {
                int1[0] = to!int(varname.length); outfile.rawWrite(int1);
                outfile.rawWrite(to!(char[])(varname));
            }
            int2[0] = myConfig.dimensions; int2[1] = to!int(ncells); outfile.rawWrite(int2);
            foreach(cell; cells) { cell.write_values_to_raw_binary(outfile); }
            outfile.close();
            break;
        default:
            auto outfile = new GzipOut(filename);
            auto writer = appender!string();
            formattedWrite(writer, "unstructured_grid_flow 1.0\n");
            formattedWrite(writer, "label: %s\n", label);
            formattedWrite(writer, "sim_time: %.18e\n", sim_time);
            formattedWrite(writer, "variables: %d\n", myConfig.flow_variable_list.length);
            // Variable list for cell on one line.
            foreach(varname; myConfig.flow_variable_list) {
                formattedWrite(writer, " \"%s\"", varname);
            }
            formattedWrite(writer, "\n");
            // Numbers of cells
            formattedWrite(writer, "dimensions: %d\n", myConfig.dimensions);
            formattedWrite(writer, "ncells: %d\n", ncells);
            outfile.compress(writer.data);
            // The actual cell data.
            foreach(cell; cells) {
                outfile.compress(" " ~ cell.write_values_to_string() ~ "\n");
            }
            outfile.finish();
        } // end switch flow_format
    } // end write_solution()

    override void propagate_inflow_data_west_to_east()
    {
        string msg = "propagate_inflow_data_west_to_east() " ~ 
            "function not implemented for unstructured grid.";
        throw new FlowSolverException(msg);
    }

    override void convective_flux_phase0(bool allow_high_order_interpolation, size_t gtl=0)
    // Compute gradients of flow quantities for higher-order reconstruction, if required.
    // To be used, later, in the convective flux calculation.
    {
        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {

            // for the MLP limiter we need to first loop over the vertices
            if (myConfig.unstructured_limiter == UnstructuredLimiter.mlp) {
                foreach (vtx; vertices) {
                    FVCell[] cell_cloud;
                    foreach(cid; cellIndexListPerVertex[vtx.id]) cell_cloud ~= cells[cid];
                    vtx.gradients.store_max_min_values_for_mlp_limiter(cell_cloud, myConfig);
                }
            }
            
            foreach (c; cells) {
                c.gradients.compute_lsq_values(c.cell_cloud, c.ws, myConfig);
            }
            if (GlobalConfig.frozen_limiter == false) {
                foreach (c; cells) {
                    // It is more efficient to determine limiting factor here for some usg limiters.
                    final switch (myConfig.unstructured_limiter) {
                        case UnstructuredLimiter.van_albada:
                            // do nothing now
                            break;
                        case UnstructuredLimiter.min_mod:
                            // do nothing now
                            break;
                        case UnstructuredLimiter.mlp:
                            c.gradients.mlp_limit(c.cell_cloud, c.ws, myConfig);
                            break;
                        case UnstructuredLimiter.barth:
                            c.gradients.barth_limit(c.cell_cloud, c.ws, myConfig);
                            break;
                        case UnstructuredLimiter.heuristic_van_albada:
                            c.gradients.heuristic_van_albada_limit(c.cell_cloud, c.ws, myConfig);
                            break;
                        case UnstructuredLimiter.venkat:
                            c.gradients.venkat_limit(c.cell_cloud, c.ws, myConfig, gtl);
                            break;
                    } // end switch
                } // end foreach c
            } // if (frozen_limter == false)
        } // end if interpolation_order > 1
    } // end convective_flux-phase0()

    override void convective_flux_phase1(bool allow_high_order_interpolation, size_t gtl=0)
    // Make use of the flow gradients to actually do the high-order reconstruction
    // and then compute fluxes of conserved quantities at all faces.
    {
        if (allow_high_order_interpolation && (myConfig.interpolation_order > 1)) {
            // Fill in gradients for ghost cells so that left- and right- cells at all faces,
            // including those along block boundaries, have the latest gradient values.
            foreach (bcond; bc) {
                if (bcond.ghost_cell_data_available == false) { continue; }
                // Proceed to do some work only if we have ghost cells in which to insert gradients.
                bool found_mapped_cell_bc = false;
                foreach (gce; bcond.preReconAction) {
                    auto mygce = cast(GhostCellMappedCellCopy)gce;
                    if (mygce && !myConfig.in_mpi_context) {
                        found_mapped_cell_bc = true;
                        // There is a mapped-cell backing the ghost cell, so we can copy its gradients.
                        foreach (i, f; bcond.faces) {
                            // Only FVCell objects in an unstructured-grid are expected to have
                            // precomputed gradients.  There will be an initialized reference
                            // in the FVCell object of a structured-grid block, so we need to
                            // test and avoid copying from such a reference.
                            auto mapped_cell_grad = mygce.get_mapped_cell(i).gradients;
                            if (bcond.outsigns[i] == 1) {
                                if (mapped_cell_grad) {
                                    f.right_cell.gradients.copy_values_from(mapped_cell_grad);
                                } else {
                                    // Fall back to looking over the face for suitable gradient data.
                                    f.right_cell.gradients.copy_values_from(f.left_cell.gradients);
                                }
                            } else {
                                if (mapped_cell_grad) {
                                    f.left_cell.gradients.copy_values_from(mapped_cell_grad);
                                } else {
                                    f.left_cell.gradients.copy_values_from(f.right_cell.gradients);
                                }
                            }
                        } // end foreach f
                    } // end if (mygce)
                } // end foreach gce
                if (!found_mapped_cell_bc) {
                    // There are no other cells backing the ghost cells on this boundary.
                    // Fill in ghost-cell gradients from the other side of the face.
                    foreach (i, f; bcond.faces) {
                        if (bcond.outsigns[i] == 1) {
                            f.right_cell.gradients.copy_values_from(f.left_cell.gradients);
                        } else {
                            f.left_cell.gradients.copy_values_from(f.right_cell.gradients);
                        }
                    } // end foreach f
                } // end if !found_mapped_cell_bc
            } // end foreach bcond
        } // end if interpolation_order > 1
        //
        // At this point, we should have all gradient values up to date and we are now ready
        // to reconstruct field values and compute the convective fluxes.
        foreach (f; faces) {
            if (f.left_cell && f.right_cell) {
                lsq.interp_both(f, gtl, Lft, Rght, allow_high_order_interpolation);
                f.fs.copy_average_values_from(Lft, Rght);
                compute_interface_flux(Lft, Rght, f, myConfig, omegaz);
            } else if (f.right_cell) {
                lsq.interp_right(f, gtl, Rght, allow_high_order_interpolation);
                f.fs.copy_values_from(Rght);
                compute_flux_at_left_wall(Rght, f, myConfig, omegaz);
            } else if (f.left_cell) {
                lsq.interp_left(f, gtl, Lft, allow_high_order_interpolation);
                f.fs.copy_values_from(Lft);
                compute_flux_at_right_wall(Lft, f, myConfig, omegaz);
            } else {
                assert(0, "oops, a face without attached cells");
            }
        } // end foreach face
    } // end convective_flux-phase1()

} // end class UFluidBlock
