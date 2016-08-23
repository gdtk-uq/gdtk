// ublock.d
// Class for unstructured blocks of cells, for use within Eilmer4.
// Peter J. 2014-11-07 first serious cut.

module ublock;

import std.conv;
import std.file;
import std.json;
import std.stdio;
import std.format;
import std.string;
import std.array;
import std.math;

import util.lua;
import json_helper;
import lua_helper;
import gzip;
import geom;
import grid;
import sgrid;
import usgrid;
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
import block;
import bc;

class UBlock: Block {
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
	super(id, Grid_t.unstructured_grid, label);
	ncells = getJSONint(json_data, "ncells", 0);
	nvertices = getJSONint(json_data, "nvertices", 0);
	nfaces = getJSONint(json_data, "nfaces", 0);
	nboundaries = getJSONint(json_data, "nboundaries", 0);
	active = getJSONbool(json_data, "active", true);
	omegaz = getJSONdouble(json_data, "omegaz", 0.0);
	// Workspace for flux_calc method.
	lsq = new LsqInterpolator(dedicatedConfig[id]);
    } // end constructor from json

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
	repr ~= "UBlock(unstructured_grid, ";
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

    override void find_enclosing_cell(double x, double y, double z, ref size_t indx, ref bool found)
    {
	grid.find_enclosing_cell(x, y, z, indx, found); // delegate to the grid object
    }

    override void init_grid_and_flow_arrays(string gridFileName)
    {
	grid = new UnstructuredGrid(gridFileName, "gziptext");
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
	// Assemble array storage for finite-volume cells, etc.
	bool lsq_workspace_at_vertices = (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
	    && (myConfig.spatial_deriv_locn == SpatialDerivLocn.vertices);
	foreach (i, v; grid.vertices) {
	    auto new_vtx = new FVVertex(myConfig, lsq_workspace_at_vertices, i);
	    new_vtx.pos[0] = v;
	    vertices ~= new_vtx;
	}
	bool lsq_workspace_at_faces = (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)
	    && (myConfig.spatial_deriv_locn == SpatialDerivLocn.faces);
	foreach (i, f; grid.faces) {
	    auto new_face = new FVInterface(myConfig, lsq_workspace_at_faces, i);
	    faces ~= new_face;
	}
	foreach (i, c; grid.cells) {
	    auto new_cell = new FVCell(myConfig, i);
	    new_cell.will_have_valid_flow = true;
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
	// Work through the faces on the boundaries and add ghost cells.
	if (nboundaries != grid.nboundaries) {
	    string msg = format("Mismatch in number of boundaries: %d %d",
				nboundaries, grid.nboundaries);
	    throw new FlowSolverException(msg);
	}
	size_t ghost_cell_count = 0;
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
		// Make ghost-cell id values distinct from FVCell ids so that
		// the warning/error messages are somewhat informative. 
		FVCell ghost0 = new FVCell(myConfig, 1000000+ghost_cell_count);
		ghost_cell_count++;
		ghost0.will_have_valid_flow = bc[i].ghost_cell_data_available;
		bc[i].faces ~= my_face;
		bc[i].outsigns ~= my_outsign;
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
	    }
	}
	// At this point, all faces should have either one finite-volume cell
	// or one ghost cells attached to a side -- check that this is true.
	foreach (f; faces) {
	    bool ok = true;
	    string msg = " ";
	    if (f.is_on_boundary) {
		if (f.left_cell && f.right_cell) {
		    ok = true;
		} else {
		    ok = false;
		    msg ~= "Boundary face does not have correct number of cells per side.";
		}
	    } else {
		// not on a boundary, should have one cell only per side.
		if (f.left_cell && f.right_cell) {
		    ok = true;
		} else {
		    ok = false;
		    msg ~= "Non-boundary face does not have correct number of cells per side.";
		}
	    }
	    if (!ok) {
		msg = format("After adding ghost cells to face %d: ", f.id) ~ msg;
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
	//equivalent to store_references_for_derivative_calc(size_t gtl) in sblock.d
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
		if (myConfig.include_ghost_cells_in_spatial_deriv_clouds) {
		    // this cloud allows for consistent gradients at block
		    // boundaries, necessary to have 2nd order convergence for multiblock simulations.
		    foreach(bndary_idx, boundary; grid.boundaries) { // set boundary clouds first
			BoundaryCondition bc = this.bc[bndary_idx];
			// TO_DO: currently we have to check whether the boundary interface is on a wall,
			//        using wall ghost cell data is causing issues for the viscous gradient calculations
			//        at these interfaces. We are working towards removing ghost cell data from
			//        wall type boundaries. In the mean time we will apply a temporary fix by using the
			//        point cloud stencil which does not use any ghost cell data at walls. KD 17/06/2016
			if (bc.is_wall) {
			    foreach (i, f; bc.faces) { 
				FVCell cell;
				// store the interface information
				f.cloud_pos ~= &(f.pos);
				f.cloud_fs ~= f.fs;
				if (bc.outsigns[i] == 1) {
				    // store "neighbour cell" information
				    f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.left_cell.fs;
				    cell = f.left_cell;
				} else {
				    // store "neighbour cell" information
				    f.cloud_pos ~= &(f.right_cell.pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.right_cell.fs;
				    cell = f.right_cell;
				}
				// now grab the remaining cloud points
				foreach (other_face; cell.iface) { // loop around cell interfaces
				    if(other_face.id == f.id) continue; // skip current interface
				    bool prevent_dup = false; // prevents double dipping in 3D simulations
				    foreach (vtx_other_face; other_face.vtx) { // loop around other_face vertices
					foreach (vtx_f; f.vtx) { // loop around f vertices
					    if (vtx_other_face.id == vtx_f.id && prevent_dup == false) {
						f.cloud_pos ~= &(other_face.pos); // store interface information
						f.cloud_fs ~= other_face.fs;
						prevent_dup = true;
					    } // end if (vtx_other_face.id == vtx_f.id && prevent_dup == false)
					} // end foreach (vtx_f; f.vtx)
				    } // end foreach (vtx_other_face; other_face.vtx)
				} // end foreach (other_face; cell.iface)
			    } // end foreach (i, f; bc.faces)
			} else {
			    foreach (i, f; bc.faces) {
				FVCell[] cell_list;
				// store interface
				f.cloud_pos ~= &(f.pos);
				f.cloud_fs ~= f.fs;
				if (bc.outsigns[i] == 1) {
				    // store "neighbour cell" information
				    f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.left_cell.fs;
				    cell_list~= f.left_cell.cell_cloud;
				    // store ghost0 information
				    if (f.right_cell.will_have_valid_flow) {
					f.cloud_pos ~= &(f.right_cell.pos[0]);
					f.cloud_fs ~= f.right_cell.fs;
				    }
				} else {
				    // store "neighbour cell" information
				    f.cloud_pos ~= &(f.right_cell.pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.right_cell.fs;
				    cell_list ~= f.right_cell.cell_cloud;
				    // store ghost0 information
				    if (f.left_cell.will_have_valid_flow) { 
					f.cloud_pos ~= &(f.left_cell.pos[0]);
					f.cloud_fs ~= f.left_cell.fs;
				    }
				}
				// now grab the remaining cells
				foreach (cell; cell_list) { // for each cell
				    foreach (other_face; cell.iface) { // loop around cell interfaces
					if (other_face.is_on_boundary && other_face.bc_id == bndary_idx) {
					    if(other_face.id == f.id) continue; // skip current interface
					    // store cell information
					    f.cloud_pos ~= &(cell.pos[0]); // assume gtl = 0
					    f.cloud_fs ~= cell.fs;
					    // and store it's neighbouring ghost0 information
					    if (bc.outsigns[i] == 1) {
						if (f.right_cell.will_have_valid_flow) {
						    f.cloud_pos ~= &(other_face.right_cell.pos[0]); // assume gtl = 0
						    f.cloud_fs ~= other_face.right_cell.fs;
						}
					    } else {
						if (f.left_cell.will_have_valid_flow) {
						    f.cloud_pos ~= &(other_face.left_cell.pos[0]); // assume gtl = 0
						    f.cloud_fs ~= other_face.left_cell.fs;
						} 
					    } // end else 
					}// end if (other_face.is_on_boundary && other_face.bc_id == bndary_idx)
				    } // end foreach (other_face; cell.iface)
				} // foreach (cell; cell_list)
			    } // end foreach (i, f; bc.faces)
			} // end else
		    } // end (bndary_idx, boundary; grid.boundaries)
		    foreach (i, f; faces) { // set internal face clouds
			if ( f.is_on_boundary == true) continue; // boundaries already set
			// store interface
			f.cloud_pos ~= &(f.pos);
			f.cloud_fs ~= f.fs;
			FVCell[] cell_list;
			cell_list ~= f.left_cell;
			cell_list ~= f.right_cell;
			foreach (cell; cell_list) {
			    // grab cell
			    f.cloud_pos ~= &(cell.pos[0]); // assume gtl = 0
			    f.cloud_fs ~= cell.fs;
			    // now grab the cell interfaces
			    foreach (other_face; cell.iface) {
				if (other_face.id == f.id) continue;
				bool prevent_dup = false; // prevents storing a face twice in 3D simulations
				foreach (vtx_other_face; other_face.vtx) { // loop around other_face vertices
				    foreach (vtx_f; f.vtx) { // loop around f vertices
					if (vtx_other_face.id == vtx_f.id && prevent_dup == false) {
					    // store interface information
					    f.cloud_pos ~= &(other_face.pos);
					    f.cloud_fs ~= other_face.fs;
					    prevent_dup = true;
					} // end if( vtx_other_face.id == vtx_f.id && prevent_dup == false)
				    } // end foreach (vtx_f; f.vtx)
				} // end foreach (vtx_other_face; other_face.vtx)
			    } // end foreach (other_face; cell.iface)
			} // end foreach (cell; cell_list)
		    } // end foreach (i, f; faces)
		} // end if (myConfig.include_ghost_cells_in_spatial_deriv_clouds)
		else { // else set similar cloud as structured solver
		    foreach(bndry_idx, boundary; grid.boundaries) { // set boundary clouds
			BoundaryCondition bc = this.bc[bndry_idx];
			foreach (i, f; bc.faces) { 
			    FVCell cell;
			    // store the interface information
			    f.cloud_pos ~= &(f.pos);
			    f.cloud_fs ~= f.fs;
			    if (bc.outsigns[i] == 1) {
				// store "neighbour cell" information
				f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
				f.cloud_fs ~= f.left_cell.fs;
				cell = f.left_cell;
			    } else {
				// store "neighbour cell" information
				f.cloud_pos ~= &(f.right_cell.pos[0]); // assume gtl = 0
				f.cloud_fs ~= f.right_cell.fs;
				cell = f.right_cell;
			    }
			    // now grab the remainding points
			    foreach (other_face; cell.iface) { // loop around cell interfaces
				if(other_face.id == f.id) continue; // skip current interface
				bool prevent_dup = false; // prevents double dipping in 3D simulations
				foreach (vtx_other_face; other_face.vtx) { // loop around other_face vertices
				    foreach (vtx_f; f.vtx) { // loop around f vertices
					if (vtx_other_face.id == vtx_f.id && prevent_dup == false) {
					    f.cloud_pos ~= &(other_face.pos); // store interface information
					    f.cloud_fs ~= other_face.fs;
					    prevent_dup = true;
					} // end if (vtx_other_face.id == vtx_f.id && prevent_dup == false)
				    } // end foreach (vtx_f; f.vtx)
				} // end foreach (vtx_other_face; other_face.vtx)
			    } // end foreach (other_face; cell.iface)
			} // end foreach (i, f; bc.faces)
		    } // end foreach(l, boundary; grid.boundaries) {
		    foreach (i, f; faces) { // set internal face clouds
			if ( f.is_on_boundary == true) continue; // boundaries already set
			// store interface
			f.cloud_pos ~= &(f.pos);
			f.cloud_fs ~= f.fs;
			FVCell[] cell_list;
			cell_list ~= f.left_cell;
			cell_list ~= f.right_cell;
			foreach (cell; cell_list) {
			    // grab cell
			    f.cloud_pos ~= &(cell.pos[0]); // assume gtl = 0
			    f.cloud_fs ~= cell.fs;
			    // now grab the interfaces
			    foreach (other_face; cell.iface) {
				if (other_face.id == f.id) continue;
				bool prevent_dup = false; // prevents double dipping in 3D simulations
				foreach (vtx_other_face; other_face.vtx) { // loop around other_face vertices
				    foreach (vtx_f; f.vtx) { // loop around f vertices
					if (vtx_other_face.id == vtx_f.id && prevent_dup == false) {
					    f.cloud_pos ~= &(other_face.pos); // store interface information
					    f.cloud_fs ~= other_face.fs;
					    prevent_dup = true;
					} // end if (vtx_other_face.id == vtx_f.id && prevent_dup == false)
				    } // end foreach (vtx_f; f.vtx)
				} // end foreach (vtx_other_face; other_face.vtx)
			    } // end foreach (other_face; cell.iface)
			} // end foreach (cell; cell_list)
		    } // end foreach (i, f; faces)
		} // end else
	    }  // if (myConfig.spatial_deriv_calc ==  SpatialDerivCalc.least_squares)
	    else { // set finite-difference cloud
		foreach (i, f; faces) {
		    if ( f.right_cell.will_have_valid_flow) {
			f.cloud_pos ~= &(f.right_cell.pos[0]); // assume gtl = 0
			f.cloud_fs ~= f.right_cell.fs;
		    }
		    if ( f.left_cell.will_have_valid_flow) {
			f.cloud_pos ~= &(f.left_cell.pos[0]); // assume gtl = 0
			f.cloud_fs ~= f.left_cell.fs;
		    }
		} // end foreach (i, f; faces)		
	    } // end else
	} // end switch (myConfig.spatial_deriv_locn)
    } // end init_grid_and_flow_arrays()
    
    override void compute_primary_cell_geometric_data(int gtl)
    {
	if (myConfig.dimensions == 2) {
	    // Primary cell geometry in the xy-plane.
	    foreach (i, cell; cells) {
		double vol, xyplane_area;
		switch (cell.vtx.length) {
		case 3:
		    xyplane_triangle_cell_properties(cell.vtx[0].pos[gtl], cell.vtx[1].pos[gtl],
						     cell.vtx[2].pos[gtl],
						     cell.pos[gtl], xyplane_area,
						     cell.iLength, cell.jLength, cell.L_min);
		    break;
		case 4:
		    xyplane_quad_cell_properties(cell.vtx[0].pos[gtl], cell.vtx[1].pos[gtl],
						 cell.vtx[2].pos[gtl], cell.vtx[3].pos[gtl],
						 cell.pos[gtl], xyplane_area,
						 cell.iLength, cell.jLength, cell.L_min);
		    break;
		default:
		    string msg = "compute_primary_cell_geometric_data(): ";
		    msg ~= format("Unhandled number of vertices: %d", cell.vtx.length);
		    throw new FlowSolverException(msg);
		} // end switch
		// Cell Volume.
		if ( myConfig.axisymmetric ) {
		    // Volume per radian = centroid y-ordinate * cell area
		    vol = xyplane_area * cell.pos[gtl].y;
		} else {
		    // Assume unit depth in the z-direction.
		    vol = xyplane_area;
		}
		if (vol < 0.0) {
		    string msg = text("Negative cell volume: Block ", id,
				      " vol for cell[", i, "]= ", vol);
		    throw new FlowSolverException(msg);
		}
		cell.volume[gtl] = vol;
		cell.areaxy[gtl] = xyplane_area;
		cell.kLength = 0.0;
	    }
	    // Face geometry in the xy-plane.
	    foreach (f; faces) {
		double xA = f.vtx[0].pos[gtl].x;
		double yA = f.vtx[0].pos[gtl].y;
		double xB = f.vtx[1].pos[gtl].x;
		double yB = f.vtx[1].pos[gtl].y;
		double LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
		// Direction cosines for the unit normal and two tangential directions.
		if (LAB > 1.0e-12) {
		    f.n.refx = (yB - yA) / LAB;
		    f.n.refy = -(xB - xA) / LAB;
		    f.n.refz = 0.0; // normal purely in xy-plane
		    f.t2 = Vector3(0.0, 0.0, 1.0);
		    f.t1 = cross(f.n, f.t2);
		    f.length = LAB; // Length in the XY-plane.
		} else {
		    f.n = Vector3(1.0, 0.0, 0.0); // Arbitrary direction
		    f.t2 = Vector3(0.0, 0.0, 1.0);
		    f.t1 = Vector3(0.0, 1.0, 0.0);
		    f.length = 0.0; // Zero length in the xy-plane
		}
		// Mid-point and area.
		// [TODO] think about using a better estimate for Ybar.
		f.Ybar = 0.5 * (yA + yB);
		if ( myConfig.axisymmetric ) {
		    f.area[gtl] = LAB * f.Ybar; // Face area per radian.
		} else {
		    f.area[gtl] = LAB; // Assume unit depth in the Z-direction.
		}
		f.pos = (f.vtx[0].pos[gtl] + f.vtx[1].pos[gtl])/2.0;
	    } // end foreach f
	} else {
	    // Primary cell geometry in 3D.
	    foreach (i, cell; cells) {
		switch (cell.vtx.length) {
		case 4:
		    tetrahedron_properties(cell.vtx[0].pos[gtl], cell.vtx[1].pos[gtl],
					   cell.vtx[2].pos[gtl], cell.vtx[3].pos[gtl],
					   cell.pos[gtl], cell.volume[gtl]);
		    cell.L_min = pow(cell.volume[gtl], 1.0/3.0);
		    cell.iLength = cell.L_min;
		    cell.jLength = cell.L_min;
		    cell.kLength = cell.L_min;
		    break;
		case 8:
		    hex_cell_properties(cell.vtx[0].pos[gtl], cell.vtx[1].pos[gtl],
					cell.vtx[2].pos[gtl], cell.vtx[3].pos[gtl],
					cell.vtx[4].pos[gtl], cell.vtx[5].pos[gtl],
					cell.vtx[6].pos[gtl], cell.vtx[7].pos[gtl],
					cell.pos[gtl], cell.volume[gtl],
					cell.iLength, cell.jLength, cell.kLength);
		    cell.L_min = cell.iLength;
		    if (cell.jLength < cell.L_min) cell.L_min = cell.jLength;
		    if (cell.kLength < cell.L_min) cell.L_min = cell.kLength;
		    break;
		case 5:
		    throw new Exception("need to implement pyramid cell properties");
		case 6:
		    throw new Exception("need to implement wedge cell properties");
		default:
		    string msg = "compute_primary_cell_geometric_data() cells: ";
		    msg ~= format("Unhandled number of vertices: %d", cell.vtx.length);
		    throw new FlowSolverException(msg);
		} // end switch
	    } // end foreach cell
	    // Face geometry in 3D.
	    foreach (f; faces) {
		switch (f.vtx.length) {
		case 3:
		    triangle_properties(f.vtx[0].pos[gtl], f.vtx[1].pos[gtl],
					f.vtx[2].pos[gtl],
					f.pos, f.n, f.t1, f.t2, f.area[gtl]);
		    break;
		case 4:
		    quad_properties(f.vtx[0].pos[gtl], f.vtx[1].pos[gtl],
				    f.vtx[2].pos[gtl], f.vtx[3].pos[gtl],
				    f.pos, f.n, f.t1, f.t2, f.area[gtl]);
		    break;
		default:
		    string msg = "compute_primary_cell_geometric_data(), faces: ";
		    msg ~= format("Unhandled number of vertices: %d", f.vtx.length);
		    throw new FlowSolverException(msg);
		} // end switch	    
	    } // end foreach f
	} // end if myConfig.dimensions
	//
	// Position ghost-cell centres and copy cross-cell lengths.
	// Copy without linear extrapolation for the moment.
	//
	// 25-Feb-2014 Note copied from Eilmer3
	// Jason Qin and Paul Petrie-Repar have identified the lack of exact symmetry in
	// the reconstruction process at the wall as being a cause of the leaky wall
	// boundary conditions.  Note that the symmetry is not consistent with the 
	// linear extrapolation used for the positions and volumes in Eilmer3.
	foreach (bndry; grid.boundaries) {
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
	if (myConfig.viscous && (myConfig.spatial_deriv_calc == SpatialDerivCalc.least_squares)) {
	    // LSQ weights are used in the calculation of flow gradients
	    // for the viscous terms.
	    compute_leastsq_weights(gtl);
	}
	if (myConfig.interpolation_order > 1) {
	    // The LSQ linear model for the flow field is fitted using 
	    // information on the locations of the points. 
	    foreach (c; cells) {
		c.ws.assemble_and_invert_normal_matrix(c.cell_cloud, myConfig.dimensions, gtl);
	    }
	}
    } // end compute_primary_cell_geometric_data()

    override void read_grid(string filename, size_t gtl=0)
    {
	throw new FlowSolverException("read_grid function NOT implemented for unstructured grid.");
    }

    override void write_grid(string filename, double sim_time, size_t gtl=0)
    // Note that we reuse the StructuredGrid object that was created on the
    // use of read_grid().
    {
	throw new FlowSolverException("write_grid function not yet implemented for unstructured grid.");
	// [TODO]
    } // end write_grid()

    override double read_solution(string filename, bool overwrite_geometry_data)
    // Note that this function needs to be kept in sync with the BlockFlow class
    // over in flowsolution.d and with write_solution() below and with
    // write_initial_usg_flow_file_from_lua() in luaflowstate.d. 
    // Returns sim_time from file.
    {
	if (myConfig.verbosity_level >= 1) {
	    writeln("read_solution(): Start block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	string format_version;
	formattedRead(line, "unstructured_grid_flow %s", &format_version);
	if (format_version != "1.0") {
	    string msg = text("UBlock.read_solution(): " ~ "format version found: "
			      ~ format_version);
	    throw new FlowSolverException(msg); 
	}
	string myLabel;
	line = byLine.front; byLine.popFront();
	formattedRead(line, "label: %s", &myLabel);
	double sim_time;
	line = byLine.front; byLine.popFront();
	formattedRead(line, "sim_time: %g", &sim_time);
	size_t nvariables;
	line = byLine.front; byLine.popFront();
	formattedRead(line, "variables: %d", &nvariables);
	line = byLine.front; byLine.popFront();
	// ingore variableNames = line.strip().split();
	line = byLine.front; byLine.popFront();
	int my_dimensions;
	formattedRead(line, "dimensions: %d", &my_dimensions);
	line = byLine.front; byLine.popFront();
	size_t nc;
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
	return sim_time;
    } // end read_solution()

    override void write_solution(string filename, double sim_time)
    // Write the flow solution (i.e. the primary variables at the cell centers)
    // for a single block.
    // Keep this function in sync with
    // write_initial_usg_flow_file_from_lua() from luaflowstate.d and
    // write_initial_flow_file() from flowstate.d.
    {
	if (myConfig.verbosity_level >= 1) {
	    writeln("write_solution(): Start block ", id);
	}
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	formattedWrite(writer, "unstructured_grid_flow 1.0\n");
	formattedWrite(writer, "label: %s\n", label);
	formattedWrite(writer, "sim_time: %.18e\n", sim_time);
	auto gmodel = myConfig.gmodel;
	auto variable_list = variable_list_for_cell(gmodel, myConfig.include_quality,
						    myConfig.MHD, myConfig.divergence_cleaning,
						    myConfig.radiation);
	formattedWrite(writer, "variables: %d\n", variable_list.length);
	// Variable list for cell on one line.
	foreach(varname; variable_list) {
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
    } // end write_solution()

    override void compute_distance_to_nearest_wall_for_all_cells(int gtl)
    // Used for the turbulence modelling.
    {
	foreach (cell; cells) {
	    double min_distance = 1.0e30; // something arbitrarily large; will be replaced
	    double cell_half_width = 1.0;
	    size_t cell_id_at_nearest_wall = 0; 
	    foreach (bndry; grid.boundaries) {
		auto nf = bndry.face_id_list.length;
		foreach (j; 0 .. nf) {
		    auto my_face = faces[bndry.face_id_list[j]];
		    auto dp = cell.pos[gtl] - my_face.pos;
		    double distance = abs(dp);
		    if (distance < min_distance) {
			min_distance =  distance;
			auto my_outsign = bndry.outsign_list[j];
			if (my_outsign == 1) {
			    auto inside0 = my_face.left_cell;
			    dp = inside0.pos[gtl] - my_face.pos;
			    cell_half_width = abs(dp);
			    cell_id_at_nearest_wall = inside0.id;
			} else {
			    auto inside0 = my_face.right_cell;
			    dp = inside0.pos[gtl] - my_face.pos;
			    cell_half_width = abs(dp);
			    cell_id_at_nearest_wall = inside0.id;
			}
		    }
		}
	    } // end foreach bndry
	    cell.distance_to_nearest_wall = min_distance;
	    cell.half_cell_width_at_wall = cell_half_width;
	    cell.cell_at_nearest_wall = cells[cell_id_at_nearest_wall];
	} // end foreach cell
    } // end compute_distance_to_nearest_wall_for_all_cells()

    override void propagate_inflow_data_west_to_east()
    {
	string msg = "propagate_inflow_data_west_to_east() " ~ 
	    "function not implemented for unstructured grid.";
	throw new FlowSolverException(msg);
    }

    override void convective_flux_phase0()
    // Compute gradients of flow quantities for higher-order reconstruction, if required.
    // To be used, later, in the convective flux calculation.
    {
	if (myConfig.interpolation_order > 1) {
	    foreach (c; cells) {
		c.gradients.compute_lsq_values(c.cell_cloud, c.ws, myConfig);
		// It is more efficient to determine limiting factor here for some usg limiters.
		final switch (myConfig.unstructured_limiter) {
		    case UnstructuredLimiter.van_albada:
			// do nothing now
		        break;
		    case UnstructuredLimiter.min_mod:
			// do nothing now
		        break;
		    case UnstructuredLimiter.barth:
		        c.gradients.barth_limit(c.cell_cloud, c.ws, myConfig);
		        break;
		    case UnstructuredLimiter.venkat:
		        c.gradients.venkat_limit(c.cell_cloud, c.ws, myConfig);
		        break;
		} // end switch
	    } // end foreach c
	} // end if interpolation_order > 1
    } // end convective_flux-phase0()

    override void convective_flux_phase1()
    // Make use of the flow gradients to actually do the high-order reconstruction
    // and then compute fluxes of conserved quantities at all faces.
    {
	if (myConfig.interpolation_order > 1) {
	    // Fill in gradients for ghost cells so that left- and right- cells at all faces,
	    // including those along block boundaries, have the latest gradient values.
	    foreach (bcond; bc) {
		bool found_mapped_cell_bc = false;
		foreach (gce; bcond.preReconAction) {
		    if (gce.type == "MappedCellExchangeCopy") {
			found_mapped_cell_bc = true;
			// There is a mapped-cell backing the ghost cell, so we can copy its gradients.
			foreach (i, f; bcond.faces) {
			    // Only FVCell objects in an unstructured-grid are expected to have
			    // precomputed gradients.  There will be an initialized reference
			    // in the FVCell object of a structured-grid block, so we need to
			    // test and avoid copying from such a reference.
			    auto mapped_cell_grad = gce.get_mapped_cell(i).gradients;
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
		    } // end if gce.type
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
	    lsq.interp_both(f, 0, Lft, Rght); // gtl assumed 0
	    f.fs.copy_average_values_from(Lft, Rght);
	    compute_interface_flux(Lft, Rght, f, myConfig, omegaz);
	} // end foreach face
    } // end convective_flux-phase1()
    
    @nogc
    override void copy_into_ghost_cells(int destination_face,
					ref Block src_blk, int src_face, int src_orientation,
					int type_of_copy, bool with_encode,
					bool reorient_vector_quantities,
					ref const(double[]) Rmatrix)
    {
	assert(false, "copy_into_ghost_cells function not implemented for unstructured grid.");
	// [TODO]
    }

} // end class UBlock
