/**
 * grid.d -- (abstract) grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-05-13 factored out of sgrid.d
 */

module geom.grid.grid;

import std.math;
import std.stdio;
import std.conv;
import std.algorithm;
import std.format;
import ntypes.complex;
import nm.number;

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
    @nogc
    size_t single_index(size_t i, size_t j, size_t k=0) const
    in {
        assert (i < niv, "index i is invalid");
        assert (j < njv, "index j is invalid");
        assert (k < nkv, "index k is invalid");
    }
    do {
        return i + niv*(j + njv*k);
    }

    size_t[] ijk_indices(size_t indx) const
    in {
        assert ( indx < vertices.length );
    }
    do {
        size_t k = indx / (niv*njv);
        indx -= k * (niv * njv);
        size_t j = indx / niv;
        indx -= j * niv;
        return [indx, j, k];
    }

    @nogc abstract Vector3* opIndex(size_t i, size_t j, size_t k=0);
    @nogc abstract Vector3* opIndex(size_t indx);
    abstract size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const;
    @nogc abstract void copy_vtx_id_list_for_cell(ref size_t[8] vtx_list_copy, ref size_t nvtx,
                                                  size_t i, size_t j, size_t k=0) const;
    abstract size_t[] get_vtx_id_list_for_cell(size_t indx) const;
    @nogc abstract void copy_vtx_id_list_for_cell(ref size_t[8] vtx_list_copy, ref size_t nvtx,
                                                  size_t indx) const;

    abstract void read_from_gzip_file(string fileName, double scale=1.0);
    abstract void read_from_raw_binary_file(string fileName, double scale=1.0);
    void write(string fileName, string fmt)
    {
        switch (fmt) {
        case "gziptext": write_to_gzip_file(fileName); break;
        case "rawbinary": write_to_raw_binary_file(fileName); break;
        default: write_to_gzip_file(fileName);
        }
    }
    abstract void write_to_gzip_file(string fileName);
    abstract void write_to_raw_binary_file(string fileName);
    abstract void write_to_vtk_file(string fileName);
    abstract void write_to_su2_file(string fileName, double scale=1.0,
                                    bool use_gmsh_order_for_wedges=true);

    abstract size_t number_of_vertices_for_cell(size_t i);
    abstract int vtk_element_type_for_cell(size_t i);
    abstract int get_cell_type(size_t i);
    abstract Grid get_boundary_grid(size_t boundary_indx);
    abstract size_t[] get_list_of_boundary_cells(size_t boundary_indx);

    @nogc
    void compute_cell_properties(size_t indx, bool true_centroid, ref Vector3 centroid, ref number volume)
    {
        number iLen, jLen, kLen, Lmin;
        size_t nvtx; size_t[8] vtx_ids;
        copy_vtx_id_list_for_cell(vtx_ids, nvtx, indx);

        if (dimensions == 2) {
            switch (nvtx) {
            case 3:
                xyplane_triangle_cell_properties(vertices[vtx_ids[0]], vertices[vtx_ids[1]],
                                                 vertices[vtx_ids[2]],
                                                 centroid, volume, iLen, jLen, Lmin);
                return;
            case 4:
                xyplane_quad_cell_properties(vertices[vtx_ids[0]], vertices[vtx_ids[1]],
                                             vertices[vtx_ids[2]], vertices[vtx_ids[3]],
                                             centroid, volume, iLen, jLen, Lmin);
                return;
            default:
                throw new Error("Unhandled number of vertices.");
            } // end switch
        } // end: if 2D
        // else, 3D
        switch (nvtx) {
        case 4:
            tetrahedron_properties(vertices[vtx_ids[0]], vertices[vtx_ids[1]],
                                   vertices[vtx_ids[2]], vertices[vtx_ids[3]],
                                   centroid, volume);
            return;
        case 8:
            hex_cell_properties(vertices[vtx_ids[0]], vertices[vtx_ids[1]],
                                vertices[vtx_ids[2]], vertices[vtx_ids[3]],
                                vertices[vtx_ids[4]], vertices[vtx_ids[5]],
                                vertices[vtx_ids[6]], vertices[vtx_ids[7]],
                                true_centroid, centroid, volume, iLen, jLen, kLen);
            return;
        case 5:
            pyramid_properties(vertices[vtx_ids[0]], vertices[vtx_ids[1]],
                               vertices[vtx_ids[2]], vertices[vtx_ids[3]],
                               vertices[vtx_ids[4]],
                               true_centroid, centroid, volume);
            return;
        case 6:
            wedge_properties(vertices[vtx_ids[0]], vertices[vtx_ids[1]],
                             vertices[vtx_ids[2]], vertices[vtx_ids[3]],
                             vertices[vtx_ids[4]], vertices[vtx_ids[5]],
                             true_centroid, centroid, volume);
            return;
        default:
            throw new Error("Unhandled number of vertices.");
        } // end switch
    } // end compute_cell_properties()

    @nogc
    bool point_is_inside_cell(ref const(Vector3) p, size_t i)
    {
        bool inside_cell = false;
        size_t nvtx; size_t[8] vtx_id;
        copy_vtx_id_list_for_cell(vtx_id, nvtx, i);
        switch (dimensions) {
        case 1: throw new GeometryException("cell search not implemented for 1D grids");
        case 2:
            switch (nvtx) {
            case 3:
                inside_cell = inside_xy_triangle(vertices[vtx_id[0]], vertices[vtx_id[1]],
                                                 vertices[vtx_id[2]], p);
                break;
            case 4:
                inside_cell = inside_xy_quad(vertices[vtx_id[0]], vertices[vtx_id[1]],
                                             vertices[vtx_id[2]], vertices[vtx_id[3]], p);
                break;
            default:
                throw new GeometryException("invalid cell type in 2D");
            } // end switch (vtx_id.length)
            break;
        case 3:
            switch (nvtx) {
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
                inside_cell = inside_pyramid(vertices[vtx_id[0]], vertices[vtx_id[1]],
                                             vertices[vtx_id[2]], vertices[vtx_id[3]],
                                             vertices[vtx_id[4]], p);
                break;
            case 6:
                inside_cell = inside_wedge(vertices[vtx_id[0]], vertices[vtx_id[1]],
                                           vertices[vtx_id[2]], vertices[vtx_id[3]],
                                           vertices[vtx_id[4]], vertices[vtx_id[5]], p);
                break;
            default:
                throw new GeometryException("invalid cell type in 3D");
            } // end switch (vtx_id.length)
            break;
        default: throw new GeometryException("invalid number of dimensions");
        } // end switch (dimensions)
        return inside_cell;
    } // end point_is_inside_cell()

    @nogc
    Vector3 cell_barycentre(size_t indx)
    // Returns the "centre-of-mass" of the vertices defining the cell.
    {
        size_t[8] vtx_ids; size_t nvtx;
        copy_vtx_id_list_for_cell(vtx_ids, nvtx, indx);
        Vector3 cbc; cbc.set(vertices[vtx_ids[0]]);
        foreach(i; 1 .. nvtx) { cbc.add(vertices[vtx_ids[i]]); }
        cbc.scale(1.0/nvtx);
        return cbc;
    } // end cell_barycentre()

    // Bin-sorting methods to speed up the search for the enclosing cell.

    // Since we are cannot control what position we are given,
    // determine the bin index cautiously.
    @nogc
    int get_bin_ix(ref const(Vector3) p)
    {
        double dx = p.x.re - bb0.x.re;
        int ix = cast(int) (dx/deltax);
        ix = max(0, min(nbx-1, ix));
        return ix;
    }

    @nogc
    int get_bin_iy(ref const(Vector3) p)
    {
        double dy = p.y.re - bb0.y.re;
        int iy = cast(int) (dy/deltay);
        iy = max(0, min(nby-1, iy));
        return iy;
    }

    @nogc
    int get_bin_iz(ref const(Vector3) p)
    {
        int iz = 0;
        if (dimensions == 3) {
            double dz = p.z.re - bb0.z.re;
            iz = cast(int) (dz/deltaz);
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
            throw new GeometryException("Need to have nbinx, nbiny and nbinz >= 1");
        }
        nbx = nbinx; nby = nbiny; nbz = nbinz;
        // Determine the bounding box for the grid.
        Vector3 p; p = vertices[0];
        bb0 = p; bb1 = p;
        foreach (i; 1 .. nvertices) {
            p = vertices[i];
            bb0.set(fmin(p.x.re, bb0.x.re), fmin(p.y.re, bb0.y.re), fmin(p.z.re, bb0.z.re));
            bb1.set(fmax(p.x.re, bb1.x.re), fmax(p.y.re, bb1.y.re), fmax(p.z.re, bb1.z.re));
        }
        // Sizes of the bins.
        deltax = (bb1.x.re - bb0.x.re)/nbx;
        deltay = (bb1.y.re - bb0.y.re)/nby;
        deltaz = (bb1.z.re - bb0.z.re)/nbz;
        immutable very_small_delta = 1.0e-6;
        deltax = fmax(deltax, very_small_delta);
        deltay = fmax(deltay, very_small_delta);
        deltaz = fmax(deltaz, very_small_delta);
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

    @nogc
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

    @nogc
    void find_enclosing_cell_slow(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
        // Search every cell in the block.
        found = false; indx = 0;
        foreach (i; 0 .. ncells) {
            if (point_is_inside_cell(p, i)) { found = true; indx = i; break; }
        }
        return;
    } // end find_enclosing_cell_slow()

    @nogc
    void find_enclosing_cell(ref const(Vector3) p, ref size_t indx, ref bool found)
    {
        if (cells_are_sorted_into_bins) {
            find_enclosing_cell_fast(p, indx, found);
        } else {
            find_enclosing_cell_slow(p, indx, found);
        }
    } // end find_enclosing_cell()

    @nogc
    void find_nearest_cell_centre(ref const(Vector3) p, ref size_t nearestCell, ref double minDist)
    {
        Vector3 c = cell_barycentre(0);
        double d = distance_between(c, p);
        minDist = d; nearestCell = 0;
        foreach (i; 1 .. ncells) {
            c = cell_barycentre(i);
            d = distance_between(c, p);
            if (d < minDist) { minDist = d; nearestCell = i; }
        }
    } // end find_nearest_cell_centre

    void write_to_stl_file(string fileName, double scale=1.0)
    // Write the grid (unstructured or structured) to a file in STL (ASCII) format.
    // The coordinates can be scaled arbitrarily since the STL format does not carry units.
    // Also, the STL file knows nothing about connectivity, so all facets are independent.
    // PJ 2021-06-29
    {
        if (dimensions != 2) {
            throw new Exception("write_to_stl_file is only for 2D grids.");
        }
        auto f = File(fileName, "w");
        f.writefln("solid %s", label);
        void write_STL_triangle(ref const(Vector3) p0, ref const(Vector3) p1, ref const(Vector3) p2)
        {
            Vector3 p01 = p1-p0; Vector3 p02 = p2-p0;
            Vector3 n; cross(n, p01, p02); n.normalize();
            f.writefln("facet normal %.6e %.6e %.6e", n.x, n.y, n.z);
            f.writeln("  outer loop");
            f.writefln("    vertex %.6e %.6e %.6e", p0.x*scale, p0.y*scale, p0.z*scale);
            f.writefln("    vertex %.6e %.6e %.6e", p1.x*scale, p1.y*scale, p1.z*scale);
            f.writefln("    vertex %.6e %.6e %.6e", p2.x*scale, p2.y*scale, p2.z*scale);
            f.writeln("  endloop");
            f.writeln("endfacet");
        }
        size_t[8] vtx_ids; size_t nvtx; // to suite copy_vtx_id_list_for_cell()
        foreach (i; 0 .. ncells) {
            copy_vtx_id_list_for_cell(vtx_ids, nvtx, i);
            switch (nvtx) {
            case 3:
                // Single triangle.
                write_STL_triangle(vertices[vtx_ids[0]], vertices[vtx_ids[1]], vertices[vtx_ids[2]]);
                break;
            case 4:
                // Quad, split into 4 triangles.
                Vector3 mid = cell_barycentre(i);
                write_STL_triangle(vertices[vtx_ids[0]], vertices[vtx_ids[1]], mid);
                write_STL_triangle(vertices[vtx_ids[1]], vertices[vtx_ids[2]], mid);
                write_STL_triangle(vertices[vtx_ids[2]], vertices[vtx_ids[3]], mid);
                write_STL_triangle(vertices[vtx_ids[3]], vertices[vtx_ids[0]], mid);
                break;
            default:
                throw new Error(format("Unexpected number of vertices, nvtx=%d.", nvtx));
            } // end switch nvtx
        } // end foreach cell
        f.writefln("endsolid %s", label);
        f.close();
    } // end write_to_STL_file()

    /**
     * Rotate the vertices of a grid given a unit quaternion and centre of rotation.
     *
     * The quaternion should be a unit quaternion of the form:
     * q = q0 + q1 \hat{i} + q2 \hat{j} + q3 \hat{k}
     *
     * Authors: Flynn Hack and RJG
     */
    void rotate(Quaternion q, Vector3 rotation_ctr=Vector3(0.0, 0.0, 0.0))
    {
        auto norm = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
        if (fabs(norm - 1.0) > 1.0e-8) {
            throw new Exception("grid.rotate: Invalid quaternion provided. Should be unit quaternion with norm = 1");
        }

        double[] Rmat;
        Rmat.length = 9;
        Rmat[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
        Rmat[1] = 2*q[1]*q[2] - 2*q[0]*q[3];
        Rmat[2] = 2*q[1]*q[3] + 2*q[0]*q[2];
        Rmat[3] = 2*q[1]*q[2] + 2*q[0]*q[3];
        Rmat[4] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
        Rmat[5] = 2*q[2]*q[3] - 2*q[0]*q[1];
        Rmat[6] = 2*q[1]*q[3] - 2*q[0]*q[2];
        Rmat[7] = 2*q[2]*q[3] + 2*q[0]*q[1];
        Rmat[8] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]; // Implementation checked against
                                                                 // scipy rotation package.

        foreach (ref v; vertices) {
            // Translate all points so that centre is at origin.
            v = v - rotation_ctr;
            // Apply rotation
            v.apply_matrix_transform(Rmat);
            // Translate all points so that centre is where we started.
            v = v + rotation_ctr;
        }
    }

    /**
     * Rotate the vertices of a grid given a rotation angle, an axis of rotation and a centre of rotation.
     *
     * This method converts to quaternion representation and then delegates to the method
     * that accepts a quaternion.
     *
     * Authors: Flynn Hack and RJG
     */
    void rotate(double theta, Vector3 rotation_axis, Vector3 rotation_ctr=Vector3(0.0, 0.0, 0.0))
    {
        auto C = cos(theta/2.0);
        auto S = sin(theta/2.0);
        rotation_axis.normalize();
        Quaternion q = [C, S*rotation_axis.x.re, S*rotation_axis.y.re, S*rotation_axis.z.re];
        rotate(q, rotation_ctr);
    }

} // end class grid
