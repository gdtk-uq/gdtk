/**
 * grid.d -- (abstract) grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-05-13 factored out of sgrid.d
 */

module grid;

import std.math;
import std.stdio;
import geom;

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

    void find_enclosing_cell(double x, double y, double z, ref size_t indx, ref bool found)
    {
	found = false;
	indx = 0;
	auto p = Vector3(x, y, z);
	foreach (i; 0 .. ncells) {
	    bool inside_cell = false;
	    auto vtx_id = get_vtx_id_list_for_cell(i);
	    switch (dimensions) {
	    case 1: throw new Exception("cell search not implemented for 1D grids");
	    case 2:
		// In 2D, assume quad cells.
		inside_cell = inside_xy_quad(vertices[vtx_id[0]], vertices[vtx_id[1]],
					     vertices[vtx_id[2]], vertices[vtx_id[3]], p); 
		break;
	    case 3:
		// In 3D, assume hex cells with 8 vertices.
		inside_cell = inside_hexahedron(vertices[vtx_id[0]], vertices[vtx_id[1]],
						vertices[vtx_id[2]], vertices[vtx_id[3]],
						vertices[vtx_id[4]], vertices[vtx_id[5]],
						vertices[vtx_id[6]], vertices[vtx_id[7]], p); 
		break;
	    default: assert(0);
	    } // end switch (dimensions)
	    if (inside_cell) { found = true; indx = i; return; }
	} // foreach i
	return;
    } // end find_enclosing_cell()

    Vector3 cell_barycentre(size_t indx)
    // Returns the "centre-of-mass" of the vertices defining the cell.
    {
	auto cbc = Vector3(0.0, 0.0, 0.0);
	auto vtx_ids = get_vtx_id_list_for_cell(indx);
	foreach(vtx_id; vtx_ids) { cbc += vertices[vtx_id]; }
	double one_over_n_vtx = 1.0 / vtx_ids.length;
	cbc *= one_over_n_vtx;
	return cbc;
    } // end cell_barycentre()

    void find_nearest_cell_centre(double x, double y, double z,
				  ref size_t nearestCell, ref double minDist)
    {
	nearestCell = 0;
	auto p = cell_barycentre(0);
	double dx = x - p.x; double dy = y - p.y; double dz = z - p.z;
	minDist = sqrt(dx*dx + dy*dy + dz*dz);
	foreach (i; 1 .. ncells) {
	    p = cell_barycentre(i);
	    dx = x - p.x; dy = y - p.y; dz = z - p.z;
	    double d = sqrt(dx*dx + dy*dy + dz*dz);
	    if (d < minDist) {
		minDist = d;
		nearestCell = i;
	    }
	} // end foreach i
    } // end find_nearest_cell_centre

} // end class grid
