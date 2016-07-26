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
	foreach (i, v; grid.vertices) {
	    auto new_vtx = new FVVertex(myConfig, false);
	    // should never need grad_ws storage at vertices
	    new_vtx.pos[0] = v;
	    new_vtx.id = i;
	    vertices ~= new_vtx;
	}
	foreach (i, f; grid.faces) {
	    auto new_face = new FVInterface(myConfig,
					    myConfig.retain_least_squares_work_data, 
					    myConfig.spatial_deriv_retain_lsq_work_data,
					    i);
	    faces ~= new_face;
	}
	foreach (i, c; grid.cells) {
	    auto new_cell = new FVCell(myConfig);
	    new_cell.id = i;
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
		    my_face.left_cell = c; // [TODO] check that we assign only once
		    my_face.left_cells ~= c; // [TODO] eliminate this list
		} else {
		    my_face.right_cell = c;
		    my_face.right_cells ~= c; // [TODO] eliminate this list
		}
	    }
	} // end foreach cells
	// Presently, no face should have more than one cell on its left or right side.
	foreach (f; faces) {
	    if (f.left_cells.length > 1 || f.right_cells.length > 1) {
		string msg = format("Face id= %d too many attached cells: left_cells= ", f.id);
		foreach (c; f.left_cells) { msg ~= format(" %d", c.id); }
		msg ~= " right_cells= ";
		foreach (c; f.right_cells) { msg ~= format(" %d", c.id); }
		throw new FlowSolverException(msg);
	    }
	}
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
		FVCell ghost0 = new FVCell(myConfig);
		ghost0.will_have_valid_flow = bc[i].ghost_cell_data_available;
		// Make ghost-cell id values distinct from FVCell ids so that
		// the warning/error messages are somewhat informative. 
		ghost0.id = 100000 + ghost_cell_count++;
		FVCell ghost1 = new FVCell(myConfig);
		ghost1.will_have_valid_flow = bc[i].ghost_cell_data_available;
		ghost1.id = 100000 + ghost_cell_count++;
		bc[i].faces ~= my_face;
		bc[i].outsigns ~= my_outsign;
		bc[i].ghostcells ~= ghost0;
		bc[i].ghostcells ~= ghost1;
		if (my_outsign == 1) {
		    my_face.right_cell = ghost0;
		    my_face.right_cells ~= ghost0; // [TODO] eliminate these lists
		    my_face.right_cells ~= ghost1;
		} else {
		    my_face.left_cell = ghost0;
		    my_face.left_cells ~= ghost0; // [TODO] eliminate these lists
		    my_face.left_cells ~= ghost1;
		}
	    }
	}
	// At this point, all faces should have either one finite-volume cell
	// or two ghost cells attached to a side -- check that this is true.
	// [TODO] change to not having the lists.
	foreach (f; faces) {
	    bool ok = true;
	    string msg = "";
	    if (f.is_on_boundary) {
		if ((f.left_cells.length == 2 && f.right_cells.length == 1) ||
		    (f.left_cells.length == 1 && f.right_cells.length == 2)) {
		    ok = true;
		} else {
		    ok = false;
		    msg ~= "Boundary face does not have correct number of cells per side.";
		}
	    } else {
		// not on a boundary, should have one cell only per side.
		if (f.left_cells.length != 1 || f.right_cells.length != 1) {
		    ok = false;
		    msg ~= "Non-boundary face does not have correct number of cells per side.";
		}
	    }
	    if (!ok) {
		msg = format("After adding ghost cells to face %d: ", f.id) ~ msg;
		msg ~= " left_cells= ";
		foreach (c; f.left_cells) { msg ~= format(" %d", c.id); }
		msg ~= " right_cells= ";
		foreach (c; f.right_cells) { msg ~= format(" %d", c.id); }
		throw new FlowSolverException(msg);
	    }
	} // end foreach f
	// For each side of a face with a single at
	// work around the faces of the attached cell to accumulate 
	// a cloud of cells for field reconstruction prior to computing
	// the convective fluxes, if we want high-order reconstruction.
	foreach (f; faces) {
	  // we will first add the third cell to the cloud for all boundary clouds
	  if (f.left_cells.length == 2) f.left_cells ~= f.right_cells[0];
	  if (f.right_cells.length == 2) f.right_cells ~= f.left_cells[0];
	  // now fill out the internal clouds
	  if (f.left_cells.length == 1) {
		auto cell = f.left_cells[0];
		foreach (i, other_face; cell.iface) {
		    //if (other_face.id == f.id) continue;  // REMOVED 2016-05-03 KD
		    auto other_cell = (cell.outsign[i] > 0) ? other_face.right_cells[0] : other_face.left_cells[0];
		    if (other_cell.will_have_valid_flow) { f.left_cells ~= other_cell; }
		}
	    }
	    if (f.right_cells.length == 1) {
		auto cell = f.right_cells[0];
		foreach (i, other_face; cell.iface) {
		    //if (other_face.id == f.id) continue;   // REMOVED 2016-05-03 KD
		    auto other_cell = (cell.outsign[i] > 0) ? other_face.right_cells[0] : other_face.left_cells[0];
		    if (other_cell.will_have_valid_flow) { f.right_cells ~= other_cell; }
		}
	    }
	} // end foreach f
	// Finally, all faces should have either a cloud of finite-volume cells
	// or two ghost cells attached to a side -- check that this is true.
	// For 2D triangular cells, expect 3 in a cloud. Quads will have 4
	// For 2D tetrahedral cells, expect 4 in a cloud. Hex cells will have 6.
	size_t min_count = (myConfig.dimensions == 2) ? 3 : 4;
	foreach (f; faces) {
	    bool ok = true;
	    string msg = "";
	    if (f.is_on_boundary) {
		if ((f.left_cells.length == 3 && f.right_cells.length >= min_count) ||
		    (f.left_cells.length >= min_count && f.right_cells.length == 3)) {
		    ok = true;
		} else {
		    ok = false;
		    msg ~= "Boundary face does not have correct number of cells per side.";
		}
	    } else {
		// not on a boundary, should have a cloud on both sides.
		if (f.left_cells.length < min_count || f.right_cells.length < min_count) {
		    ok = false;
		    msg ~= "Non-boundary face does not have correct number of cells per side.";
		}
	    }
	    if (!ok) {
		msg = format("After adding clouds of cells to face %d: ", f.id) ~ msg;
		msg ~= " left_cells= ";
		foreach (c; f.left_cells) { msg ~= format(" %d", c.id); }
		msg ~= " right_cells= ";
		foreach (c; f.right_cells) { msg ~= format(" %d", c.id); }
		throw new FlowSolverException(msg);
	    }
	} // end foreach f
	
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
				    f.cloud_pos ~= &(f.left_cells[0].pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.left_cells[0].fs;
				    cell = f.left_cells[0];
				} else {
				    // store "neighbour cell" information
				    f.cloud_pos ~= &(f.right_cells[0].pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.right_cells[0].fs;
				    cell = f.right_cells[0];
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
				    f.cloud_pos ~= &(f.left_cells[0].pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.left_cells[0].fs;
				    cell_list~= f.left_cells;
				    // store ghost0 information
				    if (f.right_cells[0].will_have_valid_flow) {
					f.cloud_pos ~= &(f.right_cells[0].pos[0]);
					f.cloud_fs ~= f.right_cells[0].fs;
				    }
				} else {
				    // store "neighbour cell" information
				    f.cloud_pos ~= &(f.right_cells[0].pos[0]); // assume gtl = 0
				    f.cloud_fs ~= f.right_cells[0].fs;
				    cell_list ~= f.right_cells;
				    // store ghost0 information
				    if (f.left_cells[0].will_have_valid_flow) { 
					f.cloud_pos ~= &(f.left_cells[0].pos[0]);
					f.cloud_fs ~= f.left_cells[0].fs;
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
						if (f.right_cells[0].will_have_valid_flow) {
						    f.cloud_pos ~= &(other_face.right_cells[0].pos[0]); // assume gtl = 0
						    f.cloud_fs ~= other_face.right_cells[0].fs;
						}
					    } else {
						if (f.left_cells[0].will_have_valid_flow) {
						    f.cloud_pos ~= &(other_face.left_cells[0].pos[0]); // assume gtl = 0
						    f.cloud_fs ~= other_face.left_cells[0].fs;
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
			cell_list ~= f.left_cells[0];
			cell_list ~= f.right_cells[0];
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
				f.cloud_pos ~= &(f.left_cells[0].pos[0]); // assume gtl = 0
				f.cloud_fs ~= f.left_cells[0].fs;
				cell = f.left_cells[0];
			    } else {
				// store "neighbour cell" information
				f.cloud_pos ~= &(f.right_cells[0].pos[0]); // assume gtl = 0
				f.cloud_fs ~= f.right_cells[0].fs;
				cell = f.right_cells[0];
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
			cell_list ~= f.left_cells[0];
			cell_list ~= f.right_cells[0];
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
		    if ( f.right_cells[0].will_have_valid_flow) {
			f.cloud_pos ~= &(f.right_cells[0].pos[0]); // assume gtl = 0
			f.cloud_fs ~= f.right_cells[0].fs;
		    }
		    if ( f.left_cells[0].will_have_valid_flow) {
			f.cloud_pos ~= &(f.left_cells[0].pos[0]); // assume gtl = 0
			f.cloud_fs ~= f.left_cells[0].fs;
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
		default:
		    string msg = "compute_primary_cell_geometric_data() cells: ";
		    msg ~= format("Unhandled number of vertices: %d", cell.vtx.length);
		    throw new FlowSolverException(msg);
		} // end switch
	    } // end foreach cell
	    // Face geometry in 3D.
	    foreach (f; faces) {
		switch (f.vtx.length) {
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
		    auto inside0 = my_face.left_cells[0];
		    Vector3 delta = my_face.pos - inside0.pos[gtl];
		    auto ghost0 = my_face.right_cells[0];
		    ghost0.pos[gtl] = my_face.pos + delta;
		    ghost0.iLength = inside0.iLength;
		    ghost0.jLength = inside0.jLength;
		    ghost0.kLength = inside0.kLength;
		    ghost0.L_min = inside0.L_min;
		    auto ghost1 = my_face.right_cells[1];
		    ghost1.pos[gtl] = my_face.pos + 3.0*delta;
		    ghost1.iLength = inside0.iLength;
		    ghost1.jLength = inside0.jLength;
		    ghost1.kLength = inside0.kLength;
		    ghost1.L_min = inside0.L_min;
		} else {
		    auto inside0 = my_face.right_cells[0];
		    Vector3 delta = my_face.pos - inside0.pos[gtl];
		    auto ghost0 = my_face.left_cells[0];
		    ghost0.pos[gtl] = my_face.pos + delta;
		    ghost0.iLength = inside0.iLength;
		    ghost0.jLength = inside0.jLength;
		    ghost0.kLength = inside0.kLength;
		    ghost0.L_min = inside0.L_min;
		    auto ghost1 = my_face.left_cells[1];
		    ghost1.pos[gtl] = my_face.pos + 3.0*delta;
		    ghost1.iLength = inside0.iLength;
		    ghost1.jLength = inside0.jLength;
		    ghost1.kLength = inside0.kLength;
		    ghost1.L_min = inside0.L_min;
		} // end if my_outsign
	    } // end foreach j
	} // end foreach bndry
	if (myConfig.viscous) {
	    // LSQ weights are used in the calculation of flow gradients
	    // for the viscous terms.
	    compute_leastsq_geometric_weights(gtl);
	}
	if (myConfig.interpolation_order > 1 &&
	    myConfig.retain_least_squares_work_data) {
	    // The LSQ linear model for the flow field is fitted using 
	    // information on the locations of the points. 
	    foreach (f; faces) {
		lsq.assemble_and_invert_normal_matrix(f, 0, f.left_cells, f.wsL);
		lsq.assemble_and_invert_normal_matrix(f, 0, f.right_cells, f.wsR);
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
	auto variable_list = variable_list_for_cell(gmodel);
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
		    double distance = abs(cell.pos[gtl] - my_face.pos);
		    if (distance < min_distance) {
			min_distance =  distance;
			auto my_outsign = bndry.outsign_list[j];
			if (my_outsign == 1) {
			    auto inside0 = my_face.left_cells[0];
			    cell_half_width = abs(inside0.pos[gtl] - my_face.pos);
			    cell_id_at_nearest_wall = inside0.id;
			} else {
			    auto inside0 = my_face.right_cells[0];
			    cell_half_width = abs(inside0.pos[gtl] - my_face.pos);
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

    override void convective_flux()
    {
	foreach (f; faces) {
	    lsq.interp_both(f, 0, Lft, Rght); // gtl assumed 0
	    f.fs.copy_average_values_from(Lft, Rght);
	    compute_interface_flux(Lft, Rght, f, myConfig.gmodel, omegaz);
	} // end foreach face
    } // end convective_flux()

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
