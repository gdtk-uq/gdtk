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
import sgrid;
import usgrid;
import gas;
import kinetics;
import globalconfig;
import globaldata;
import flowstate;
import fluxcalc;
import viscousflux;
import fvcore;
import fvvertex;
import fvinterface;
import fvcell;
import onedinterp;
import block;
import bc;

class UBlock: Block {
public:
    size_t ncells;
    size_t nvertices;
    size_t nfaces;
    size_t nboundaries;
    UnstructuredGrid grid;

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

    override void init_grid_and_flow_arrays(string gridFileName)
    {
	grid = new UnstructuredGrid(gridFileName, "gziptext");
	if (grid.nvertices != nvertices) {
	    throw new Error(format("UnstructuredGrid: incoming grid has %d vertices " ~
				   "but expected %d vertices.", grid.nvertices, nvertices));
	}
	if (grid.nfaces != nfaces) {
	    throw new Error(format("UnstructuredGrid: incoming grid has %d faces " ~
				   "but expected %d faces.", grid.nfaces, nfaces));
	}
	if (grid.ncells != ncells) {
	    throw new Error(format("UnstructuredGrid: incoming grid has %d cells " ~
				   "but expected %d cells.", grid.ncells, ncells));
	}
	// Assemble array storage for finite-volume cells, etc.
	foreach (i, v; grid.vertices) {
	    auto new_vtx = new FVVertex(myConfig.gmodel);
	    new_vtx.pos[0] = v;
	    new_vtx.id = i;
	    vertices ~= new_vtx;
	}
	foreach (i, f; grid.faces) {
	    auto new_face = new FVInterface(myConfig.gmodel);
	    new_face.id = i;
	    faces ~= new_face;
	}
	foreach (i, c; grid.cells) {
	    auto new_cell = new FVCell(myConfig);
	    new_cell.id = i;
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
		throw new Error(format("Mismatch in face_id_list, outsign_list lengths: %d %d",
				       grid.cells[i].face_id_list.length,
				       grid.cells[i].outsign_list.length));
	    }
	    foreach (j; 0 .. nf) {
		auto my_face = faces[grid.cells[i].face_id_list[j]];
		auto my_outsign = grid.cells[i].outsign_list[j];
		c.iface ~= my_face;
		c.outsign ~= to!double(my_outsign);
		if (my_outsign == 1) {
		    my_face.left_cells ~= c;
		} else {
		    my_face.right_cells ~= c;
		}
	    }
	} // end foreach cells
	// Presently, no face should have more than one cell on its left or right side.
	foreach (f; faces) {
	    if (f.left_cells.length > 1 || f.right_cells.length > 1) {
		string msg = format("Face id= %d too many attached cells: left_cells= ", f.id);
		foreach (c; f.left_cells) { msg ~= to!string(c.id); }
		msg ~= " right_cells= ";
		foreach (c; f.right_cells) { msg ~= to!string(c.id); }
		throw new Error(msg);
	    }
	}
	// Work through the faces on the boundaries and add ghost cells.
	if (nboundaries != grid.nboundaries) {
	    throw new Error(format("Mismatch in number of boundaries: %d %d",
				   nboundaries, grid.nboundaries));
	}
	foreach (i, bndry; grid.boundaries) {
	    auto nf = bndry.face_id_list.length;
	    if (nf != bndry.outsign_list.length) {
		throw new Error(format("Mismatch in face_id_list, outsign_list lengths: %d %d",
				       bndry.face_id_list.length,
				       bndry.outsign_list.length));
	    }
	    foreach (j; 0 .. nf) {
		FVInterface my_face = faces[bndry.face_id_list[j]];
		int my_outsign = bndry.outsign_list[j];
		BasicCell ghost0 = new BasicCell(myConfig);
		BasicCell ghost1 = new BasicCell(myConfig);
		bc[i].faces ~= my_face;
		bc[i].outsigns ~= my_outsign;
		bc[i].ghostcells ~= ghost0;
		bc[i].ghostcells ~= ghost1;
		if (my_outsign == 1) {
		    my_face.right_cells ~= ghost0;
		    my_face.right_cells ~= ghost1;
		} else {
		    my_face.left_cells ~= ghost0;
		    my_face.left_cells ~= ghost1;
		}
	    }
	}
	//
	// [TODO] for each face, work arround the faces in the attached cell to
	// accumulate a cloud of cells for field reconstruction prior to
	// computing the convective fluxes. (If we want high-order reconstruction.)
	// 
	// [TODO] store references into the FVVertex objects for derivative calc.
	//
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
		    throw new Error(msg);
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
		    throw new Error(text("Negative cell volume: Block ", id,
					 " vol for cell[", i, "]= ", vol));
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
		    throw new Error(msg);
		} // end switch
	    } // end foreach cell
	    // Face geometry in 3D.
	    foreach (f; faces) {
		switch (f.vtx.length) {
		case 4:
		    quad_properties(f.vtx[0].pos[gtl], f.vtx[1].pos[gtl],
				    f.vtx[2].pos[gtl], f.vtx[2].pos[gtl],
				    f.pos, f.n, f.t1, f.t2, f.area[gtl]);
		    break;
		default:
		    string msg = "compute_primary_cell_geometric_data(), faces: ";
		    msg ~= format("Unhandled number of vertices: %d", f.vtx.length);
		    throw new Error(msg);
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
		    ghost0.iLength = inside0.jLength;
		    ghost0.iLength = inside0.kLength;
		    ghost0.iLength = inside0.L_min;
		    auto ghost1 = my_face.right_cells[1];
		    ghost1.pos[gtl] = my_face.pos + 3.0*delta;
		    ghost1.iLength = inside0.iLength;
		    ghost1.iLength = inside0.jLength;
		    ghost1.iLength = inside0.kLength;
		    ghost1.iLength = inside0.L_min;
		} else {
		    auto inside0 = my_face.right_cells[0];
		    Vector3 delta = my_face.pos - inside0.pos[gtl];
		    auto ghost0 = my_face.left_cells[0];
		    ghost0.pos[gtl] = my_face.pos + delta;
		    ghost0.iLength = inside0.iLength;
		    ghost0.iLength = inside0.jLength;
		    ghost0.iLength = inside0.kLength;
		    ghost0.iLength = inside0.L_min;
		    auto ghost1 = my_face.left_cells[1];
		    ghost1.pos[gtl] = my_face.pos + 3.0*delta;
		    ghost1.iLength = inside0.iLength;
		    ghost1.iLength = inside0.jLength;
		    ghost1.iLength = inside0.kLength;
		    ghost1.iLength = inside0.L_min;
		} // end if my_outsign
	    } // end foreach j
	} // end foreach bndry
    } // end compute_primary_cell_geometric_data()

    override void read_grid(string filename, size_t gtl=0)
    {
	throw new Error("read_grid function NOT implemented for unstructured grid.");
    }

    override void write_grid(string filename, double sim_time, size_t gtl=0)
    // Note that we reuse the StructuredGrid object that was created on the
    // use of read_grid().
    {
	throw new Error("write_grid function not yet implemented for unstructured grid.");
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
	    throw new Error("UBlock.read_solution(): " ~
			    "format version found: " ~ format_version); 
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
	    throw new Error(text("For block[", id, "] we have a mismatch in solution size.",
				 " Have read nc=", nc, " ncells=", ncells));
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
	formattedWrite(writer, "sim_time: %20.12e\n", sim_time);
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
	throw new Error("propagate_inflow_data_west_to_east() " ~ 
			"function not implemented for unstructured grid.");
    }

    override void convective_flux()
    {
	foreach (f; faces) {
	    // For now, copy flow states without high-order (any) reconstruction.
	    Lft.copy_values_from(f.left_cells[0].fs);
	    Rght.copy_values_from(f.right_cells[0].fs);
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
