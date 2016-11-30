/**
 * grid.d -- (abstract) grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-05-13 factored out of sgrid.d
 */

module grid;

import std.math;
import std.stdio;
import std.conv;
import std.algorithm;
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
    // For both structured and unstructured grids, we can index vertices with i,j,k indices.
    // For structured grids, the numbers have an obvious significance.
    // For unstructured grids, niv==vertices.length, njv==1, nkv==1
    size_t niv, njv, nkv;
    // The following array can be used to get a list of faces attached to a vertex.
    // We have only filled in its data for unstructured grids, however,
    // we need to access the data via a this base class.
    // I think that my abstractions are leaking.
    size_t[][] faceIndexListPerVertex;
    //
    // For fast searches, we should sort the cells according
    // to location, so that we may limit the number of cells that
    // we need to check for any given location.
    Vector3 bb0, bb1; // bounding box for the entire grid.
    size_t[][][][] bins; // Array of bins for holding lists of cell indices.
    size_t nbx, nby, nbz; // number of bins in each coordinate direction
    double deltax, deltay, deltaz; // size of bin in each coordinate direction
    bool cells_are_sorted_into_bins = false;
    
    this(Grid_t grid_type, int dimensions, string label="")
    {
	this.grid_type = grid_type;
	this.dimensions = dimensions;
	this.label = label;
    }

    // Unified indexing.
    size_t single_index(size_t i, size_t j, size_t k=0) const
    in {
	assert (i < niv, text("index i=", i, " is invalid, niv=", niv));
	assert (j < njv, text("index j=", j, " is invalid, njv=", njv));
	assert (k < nkv, text("index k=", k, " is invalid, nkv=", nkv));
    }
    body {
	return i + niv*(j + njv*k);
    }

    size_t[] ijk_indices(size_t indx) const
    in {
	assert ( indx < vertices.length );
    }
    body {
	size_t k = indx / (niv*njv);
	indx -= k * (niv * njv);
	size_t j = indx / niv;
	indx -= j * niv;
	return [indx, j, k];
    }
	    
    abstract Vector3* opIndex(size_t i, size_t j, size_t k=0);
    abstract Vector3* opIndex(size_t indx);
    abstract size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const; 
    abstract size_t[] get_vtx_id_list_for_cell(size_t indx) const;
    abstract void read_from_gzip_file(string fileName, double scale=1.0);
    abstract void write_to_gzip_file(string fileName);
    abstract void write_to_vtk_file(string fileName);
    abstract size_t number_of_vertices_for_cell(size_t i);
    abstract int vtk_element_type_for_cell(size_t i);
    abstract Grid get_boundary_grid(size_t boundary_indx);
    abstract size_t[] get_list_of_boundary_cells(size_t boundary_indx);
    abstract double cell_volume(size_t indx);

    bool point_is_inside_cell(ref const(Vector3) p, size_t i)
    {
	bool inside_cell = false;
	size_t[] vtx_id = get_vtx_id_list_for_cell(i);
	switch (dimensions) {
	case 1: throw new Exception("cell search not implemented for 1D grids");
	case 2:
	    switch (vtx_id.length) {
	    case 3:
		inside_cell = inside_xy_triangle(vertices[vtx_id[0]], vertices[vtx_id[1]],
						 vertices[vtx_id[2]], p);
		break;
	    case 4:
		inside_cell = inside_xy_quad(vertices[vtx_id[0]], vertices[vtx_id[1]],
					     vertices[vtx_id[2]], vertices[vtx_id[3]], p);
		break;
	    default:
		assert(0);
	    } // end switch (vtx_id.length)
	    break;
	case 3:
	    switch (vtx_id.length) {
	    case 4:
		inside_cell = inside_tetrahedron(vertices[vtx_id[0]], vertices[vtx_id[1]],
						 vertices[vtx_id[2]], vertices[vtx_id[3]], p);
		break;
	    case 8:
		inside_cell = inside_hexahedron(vertices[vtx_id[0]], vertices[vtx_id[1]],
						vertices[vtx_id[2]], vertices[vtx_id[3]],
						vertices[vtx_id[4]], vertices[vtx_id[5]],
						vertices[vtx_id[6]], vertices[vtx_id[7]], p); 
		break;
	    case 5:
		throw new Exception("need to implement inside pyramid cell");
	    case 6:
		throw new Exception("need to implement inside wedge cell");
	    default:
		assert(0);
	    } // end switch (vtx_id.length)
	    break;
	default: assert(0);
	} // end switch (dimensions)
	return inside_cell;
    } // end point_is_inside_cell()

    Vector3 cell_barycentre(size_t indx)
    // Returns the "centre-of-mass" of the vertices defining the cell.
    {
	Vector3 cbc = Vector3(0.0, 0.0, 0.0);
	size_t[] vtx_ids = get_vtx_id_list_for_cell(indx);
	foreach(vtx_id; vtx_ids) { cbc += vertices[vtx_id]; }
	double one_over_n_vtx = 1.0 / vtx_ids.length;
	cbc *= one_over_n_vtx;
	return cbc;
    } // end cell_barycentre()

    // Bin-sorting methods to speed up the search for the enclosing cell.
    
    // Since we are cannot control what position we are given,
    // determine the bin index cautiously.
    int get_bin_ix(ref const(Vector3) p)
    {
	double dx = p.x - bb0.x;
	int ix = to!int(dx/deltax);
	ix = max(0, min(nbx-1, ix));
	return ix;
    }
    
    int get_bin_iy(ref const(Vector3) p)
    {
	double dy = p.y - bb0.y;
	int iy = to!int(dy/deltay);
	iy = max(0, min(nby-1, iy));
	return iy;
    }
    
    int get_bin_iz(ref const(Vector3) p)
    {
	int iz = 0;
	if (dimensions == 3) {
	    double dz = p.z - bb0.z;
	    iz = to!int(dz/deltaz);
	    iz = max(0, min(nbz-1, iz));
	}
	return iz;
    }
    
    void sort_cells_into_bins(size_t nbinx=5, size_t nbiny=5, size_t nbinz=5)
    // We have this as method separate from the constructor because
    // it may need to be called more than once for a moving grid.
    {
	if (dimensions == 2) { nbinz = 1; }
	// Numbers of bins in each coordinate direction.
	if (nbinx < 1 || nbiny < 1 || nbinz < 1) {
	    throw new Exception("Need to have nbinx, nbiny and nbinz >= 1");
	}
	nbx = nbinx; nby = nbiny; nbz = nbinz;
	// Determine the bounding box for the grid.
	Vector3 p = vertices[0];
	bb0 = p; bb1 = p;
	foreach (i; 1 .. nvertices) {
	    p = vertices[i];
	    bb0.refx = fmin(p.x, bb0.x); bb0.refy = fmin(p.y, bb0.y); bb0.refz = fmin(p.z, bb0.z);
	    bb1.refx = fmax(p.x, bb1.x); bb1.refy = fmax(p.y, bb1.y); bb1.refz = fmax(p.z, bb1.z);
	}
	// Sizes of the bins.
	deltax = (bb1.x - bb0.x)/nbx; deltay = (bb1.y - bb0.y)/nby; deltaz = (bb1.z - bb0.z)/nbz;
	// Now, set up the array of bins and sort the cells into those bins.
	bins.length = nbx;
	foreach (ix; 0 .. nbx) {
	    bins[ix].length = nby;
	    foreach (iy; 0 .. nby) { bins[ix][iy].length = nbz; }
	}
	foreach (i; 0 .. ncells) {
	    p = cell_barycentre(i);
	    int ix = get_bin_ix(p);
	    int iy = get_bin_iy(p);
	    int iz = get_bin_iz(p);
	    bins[ix][iy][iz] ~= i;
	}
	cells_are_sorted_into_bins = true;
    } // end sort_cells_into_bins()

    void find_enclosing_cell_fast(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
	// Search for the cell only in the near-by bins.
	found = false; indx = 0;
	if (!inside_bounding_box(p, bb0, bb1, dimensions)) return;
	//
	// So, the point is inside or on the grid's bounding box; look more closely.
	// Pick the most-likely bin, based on position.
	int ix0 = get_bin_ix(p);
	int iy0 = get_bin_iy(p);
	int iz0 = get_bin_iz(p);
	// Visit that bin first.
	foreach (i; bins[ix0][iy0][iz0]) {
	    if (point_is_inside_cell(p, i)) { found = true; indx = i; return; }
	}
        // If we reach this point, extend the search one bin further in each direction.
	int ix_start = max(ix0-1, 0); int ix_end = min(ix0+1, nbx-1);
	int iy_start = max(iy0-1, 0); int iy_end = min(iy0+1, nby-1);
	int iz_start = max(iz0-1, 0); int iz_end = min(iz0+1, nbz-1);
	for (auto iz = iz_start; iz <= iz_end; ++iz) {
	    for (auto iy = iy_start; iy <= iy_end; ++iy) {
		for (auto ix = ix_start; ix <= ix_end; ++ix) {
		    // Search only bins we have not already visited.
		    if (ix == ix0 && iy == iy0 && iz == iz0) continue;
		    foreach (i; bins[ix][iy][iz]) {
			if (point_is_inside_cell(p, i)) { found = true; indx = i; return; }
		    }
		}
	    }
	}
	// If we arrive here, give up anyway and look no further.
	return;
    } // end find_enclosing_cell_fast()

    void find_enclosing_cell_slow(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
	// Search every cell in the block.
	found = false; indx = 0;
	foreach (i; 0 .. ncells) {
	    if (point_is_inside_cell(p, i)) { found = true; indx = i; break; }
	}
	return;
    } // end find_enclosing_cell_slow()

    void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
	if (cells_are_sorted_into_bins) {
	    find_enclosing_cell_fast(p, indx, found);
	} else {
	    find_enclosing_cell_slow(p, indx, found);
	}
    } // end find_enclosing_cell()

    void find_nearest_cell_centre(ref const(Vector3) p, ref size_t nearestCell, ref double minDist)
    {
	Vector3 dp = cell_barycentre(0); dp -= p;
	double d = abs(dp);
	minDist = d; nearestCell = 0;
	foreach (i; 1 .. ncells) {
	    dp = cell_barycentre(i); dp -= p;
	    d = abs(dp);
	    if (d < minDist) { minDist = d; nearestCell = i; }
	}
    } // end find_nearest_cell_centre
    
} // end class grid
