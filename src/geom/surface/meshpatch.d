// meshpatch.d

module geom.surface.meshpatch;

import std.conv;
import std.algorithm;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;
import geom.grid;


class MeshPatch : ParametricSurface {
public:
    size_t niv;
    size_t njv;
    Vector3[][] mesh;
    Vector3 p00, p10, p11, p01;    // corners

    this(ref StructuredGrid grid)
    {
        if (grid.nkv != 1) {
            throw new Error("MeshPatch expected a grid with nkv == 1");
        }
        // Keep our own copy of the underlying grid.
        niv = grid.niv;
        njv = grid.njv;
        mesh.length = niv;
        foreach(i; 0 .. niv) {
            foreach (j; 0 .. njv) { mesh[i] ~= *grid[i,j]; }
        }
        p00 = mesh[0][0];
        p10 = mesh[niv-1][0];
        p01 = mesh[0][njv-1];
        p11 = mesh[niv-1][njv-1];
    }

    this(ref const(Vector3[][]) grid)
    {
        niv = grid.length;
        njv = grid[0].length;
        mesh.length = niv;
        foreach(i; 0 .. niv) {
            foreach (j; 0 .. njv) { mesh[i] ~= grid[i][j]; }
        }
        p00 = grid[0][0];
        p10 = grid[niv-1][0];
        p01 = grid[0][njv-1];
        p11 = grid[niv-1][njv-1];
    }

    this(ref const(MeshPatch) other)
    {
        niv = other.niv;
        njv = other.njv;
        mesh.length = niv;
        foreach(i; 0 .. niv) {
            foreach (j; 0 .. njv) { mesh[i] ~= other.mesh[i][j]; }
        }
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
    }

    override MeshPatch dup() const
    {
        return new MeshPatch(this.mesh);
    }

    override Vector3 opCall(double r, double s) const
    {
        // Interpolate within the background mesh.
        // This involves finding the relevant cell and
        // using a bilinear interpolation within that cell.
        double dr = 1.0 / (niv-1);
        double ds = 1.0 / (njv-1);
        int i = to!int(r / dr); i = max(0, i); i = min(i, niv-2);
        int j = to!int(s / ds); j = max(0, j); j = min(j, njv-2);
        // Parametric coordinate within the coarse cell.
        double local_r = (r - dr * i) / dr;
        double local_s = (s - ds * j) / ds;
        // BiLinear interpolation.
        Vector3 p = (1.0-local_r)*(1.0-local_s) * mesh[i][j] +
            (1.0-local_r)*local_s * mesh[i][j+1] +
            local_r*(1.0-local_s) * mesh[i+1][j] +
            local_r*local_s * mesh[i+1][j+1];
        return p;
    }

    override string toString() const
    {
        string repr = "MeshPatch(";
        repr ~= ")";
        return repr;
    }
} // end class MeshPatch

// [TODO] version(meshpatch_test)
