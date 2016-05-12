/**
 * grid.d -- (abstract) grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-05-13 factored out of sgrid.d
 */

module grid;

// import std.algorithm;
// import std.string;
// import std.array;
// import std.conv;
// import std.stdio;
// import std.format;
// import std.math;
// import gzip;

import geom;
// import gpath;
// import surface;
// import volume;
// import univariatefunctions;

//-----------------------------------------------------------------

enum Grid_t {structured_grid, unstructured_grid}

string gridTypeName(Grid_t gt)
{
    final switch (gt) {
    case Grid_t.structured_grid: return "structured_grid";
    case Grid_t.unstructured_grid: return "unstructured_grid";
    }
}

Grid_t gridTypeFromName(string name)
{
    switch (name) {
    case "structured_grid": return Grid_t.structured_grid;
    case "unstructured_grid": return Grid_t.unstructured_grid;
    default: throw new Error("Unknown type of grid: " ~ name);
    }
}

class Grid {
    Grid_t grid_type;
    int dimensions; // 1, 2 or 3
    string label;
    size_t ncells;
    size_t nvertices;
    Vector3[] vertices;
    size_t[] vtx_id;
    
    this(Grid_t grid_type, int dimensions, string label="")
    {
	this.grid_type = grid_type;
	this.dimensions = dimensions;
	this.label = label;
    }

    abstract Vector3* opIndex(size_t i, size_t j, size_t k=0);
    abstract Vector3* opIndex(size_t indx);
    abstract size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const; 
    abstract size_t[] get_vtx_id_list_for_cell(size_t indx) const;
    abstract void read_from_gzip_file(string fileName);
    abstract void write_to_gzip_file(string fileName);
    abstract void write_to_vtk_file(string fileName);
}

//-----------------------------------------------------------------
