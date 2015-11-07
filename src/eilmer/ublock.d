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
	super(id, label);
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
	    string json_key = format("face_%d", boundary);
	    auto bc_json_data = json_data[json_key];
	    // [TODO] bc ~= make_BC_from_json(bc_json_data, id, boundary);
	}
    } // end init_boundary_conditions()

    override string toString() const
    {
	char[] repr;
	repr ~= "UBlock(";
	repr ~= "id=" ~ to!string(id);
	repr ~= " label=\"" ~ label ~ "\"";
	repr ~= ", active=" ~ to!string(active);
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

    override void init_grid_and_flow_arrays(string gridFileName)
    {
	grid = new UnstructuredGrid(gridFileName, "gziptext");
	throw new Error("init_grid_and_flow_arrays not yet implemented for unstructured grid.");
	// [TODO] assemble_arrays();
	// [TODO] bind_interfaces_vertices_and_cells();
	// [TODO] compute_primary_cell_geometric_data(0);
	// [TODO] compute ghost-cell details for boundaries
	// [TODO] assign_flow_locations_for_derivative_calc(0);
    }

    override double read_solution(string filename, bool overwrite_geometry_data)
    // Note that the position data is read into grid-time-level 0
    // by scan_values_from_string(). 
    // Returns sim_time from file.
    {
	size_t nc;
	double sim_time;
	if (myConfig.verbosity_level >= 1) {
	    writeln("read_solution(): Start block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, " %g", &sim_time);
	line = byLine.front; byLine.popFront();
	// ignore second line; it should be just the names of the variables
	// [TODO] We should test the incoming strings against the current variable names.
	line = byLine.front; byLine.popFront();
	formattedRead(line, "%d", &nc);
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
    // This is almost Tecplot POINT format.
    {
	if (myConfig.verbosity_level >= 1) {
	    writeln("write_solution(): Start block ", id);
	}
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	formattedWrite(writer, "%20.12e\n", sim_time);
	outfile.compress(writer.data);
	writer = appender!string();
	foreach(varname; variable_list_for_cell(myConfig.gmodel)) {
	    formattedWrite(writer, " \"%s\"", varname);
	}
	formattedWrite(writer, "\n");
	outfile.compress(writer.data);
	writer = appender!string();
	formattedWrite(writer, "%d\n", ncells);
	outfile.compress(writer.data);
	foreach(i; 0 .. ncells) {
	    outfile.compress(" " ~ cells[i].write_values_to_string() ~ "\n");
	}
	outfile.finish();
    } // end write_solution()

    override void convective_flux()
    {
	throw new Error("convective_flux function not yet implemented for unstructured grid.");
	// [TODO]
    } // end convective_flux()

} // end class UBlock
