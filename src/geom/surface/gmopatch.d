// gmopatch.d
// A grid-metric-optimized surface patch that internally redistributes
// the r,s-interpolation parameters to get the "best" test grid.
// PJ 2020-10-15

module geom.surface.gmopatch;

import std.conv;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;
import geom.surface.coonspatch;
import geom.grid.sgrid;
import geom.misc.univariatefunctions;


class GMOPatch : ParametricSurface {
public:
    Path north, east, south, west; // bounding paths
    Vector3 p00, p10, p11, p01;    // corners
    // The following initial values for the interior parameters
    // will be adjusted during construction of a new patch.
    double[4][4] r_grid = [[0.0, 1.0/3, 2.0/3, 1.0],
                           [0.0, 1.0/3, 2.0/3, 1.0],
                           [0.0, 1.0/3, 2.0/3, 1.0],
                           [0.0, 1.0/3, 2.0/3, 1.0]];
    double[4][4] s_grid = [[0.0, 0.0, 0.0, 0.0],
                           [1.0/3, 1.0/3, 1.0/3, 1.0/3],
                           [2.0/3, 2.0/3, 2.0/3, 2.0/3],
                           [1.0, 1.0, 1.0, 1.0]];

    this(in Vector3 p00, in Vector3 p10, in Vector3 p11, in Vector3 p01,
         int niv=11, int njv=11)
    {
        north = new Line(p01, p11);
        east = new Line(p10, p11);
        south = new Line(p00, p10);
        west = new Line(p00, p01);
        this(south, north, west, east, niv, njv);
    }

    this(in Path south, in Path north, in Path west, in Path east,
         int niv=11, int njv=11)
    // The particular order for the boundary surfaces goes way back
    // to the original grid generation paper, so it doesn't match the
    // default order of NESW in the rest of the flow code.
    // niv and njv are the number of vertices for the test grid.
    {
        this.north = north.dup();
        this.east = east.dup();
        this.south = south.dup();
        this.west = west.dup();
        p00 = south(0.0);
        p10 = south(1.0);
        p01 = north(0.0);
        p11 = north(1.0);
        // Check alternate evaluation of corners for consistency.
        Vector3 p00_alt = west(0.0);
        Vector3 p10_alt = east(0.0);
        Vector3 p01_alt = west(1.0);
        Vector3 p11_alt = east(1.0);
        if (!approxEqualVectors(p00, p00_alt)) {
            throw new Error(text("GMOPatch open corner p00= ", p00,
                                 " p00_alt= ", p00_alt));
        }
        if (!approxEqualVectors(p10, p10_alt)) {
            throw new Error(text("GMOPatch open corner p10= ", p10,
                                 " p10_alt= ", p10_alt));
        }
        if (!approxEqualVectors(p11, p11_alt)) {
            throw new Error(text("GMOPatch open corner p11= ", p11,
                                 " p11_alt= ", p11_alt));
        }
        if (!approxEqualVectors(p01, p01_alt)) {
            throw new Error(text("GMOPatch open corner p01= ", p01,
                                 " p01_alt= ", p01_alt));
        }
        // Now we need to set up a temporary grid and adjust the internal
        // distribution of the interpolation parameters to get a nice grid.
        // We hang onto the redistributed parameters.
        auto my_patch = new CoonsPatch(south, north, west, east);
        auto cf = [new LinearFunction(), new LinearFunction(),
                   new LinearFunction(), new LinearFunction()];
        auto my_grid = new StructuredGrid(my_patch, niv, njv, cf, r_grid, s_grid);
        my_grid.determine_rs_grids(my_patch, cf, r_grid, s_grid);
    }

    this(ref const(GMOPatch) other)
    {
        this.north = other.north.dup();
        this.east = other.east.dup();
        this.south = other.south.dup();
        this.west = other.west.dup();
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
        r_grid[][] = other.r_grid[][];
        s_grid[][] = other.s_grid[][];
    }

    override GMOPatch dup() const
    {
        return new GMOPatch(this.south, this.north, this.west, this.east);
    }

    override Vector3 opCall(double r, double s) const
    {
        // Set up redistribution function to r,s.
        void remap(double r, double s, out double r_star, out double s_star)
        {
            // Tensor-product, cubic-Bezier interpolation of the gridded parameter values.
            double[4] rr; double[4] sr;
            foreach (j; 0 .. 4) {
                rr[j] = (1.0-r)^^3 * r_grid[0][j] + 3.0*r*(1.0-r)^^2 * r_grid[1][j] +
                    3.0*(1.0-r)*r*r * r_grid[2][j] + r^^3 * r_grid[3][j];
                sr[j] = (1.0-r)^^3 * s_grid[0][j] + 3.0*r*(1.0-r)^^2 * s_grid[1][j] +
                    3.0*(1.0-r)*r*r * s_grid[2][j] + r^^3 * s_grid[3][j];
            }
            r_star = (1.0-s)^^3 * rr[0] + 3.0*s*(1.0-s)^^2 * rr[1] + 3.0*(1.0-s)*s*s * rr[2] + s^^3 * rr[3];
            s_star = (1.0-s)^^3 * sr[0] + 3.0*s*(1.0-s)^^2 * sr[1] + 3.0*(1.0-s)*s*s * sr[2] + s^^3 * sr[3];
        }
        double rstar; double sstar;
        remap(r, s, rstar, sstar);
        Vector3 south_r = south(rstar);
        Vector3 north_r = north(rstar);
        Vector3 west_s = west(sstar);
        Vector3 east_s = east(sstar);
        Vector3 p = (1.0-sstar)*south_r + sstar*north_r + (1.0-rstar)*west_s + rstar*east_s -
            ((1.0-rstar)*(1.0-sstar)*p00 + (1.0-rstar)*sstar*p01 + rstar*(1.0-sstar)*p10 + rstar*sstar*p11);
        return p;
    }

    override string toString() const
    {
        return "GMOPatch(south=" ~ to!string(south) ~
            ", north=" ~ to!string(north) ~
            ", west=" ~ to!string(west) ~
            ", east=" ~ to!string(east) ~
            ", r_grid=" ~ to!string(r_grid) ~
            ", s_grid=" ~ to!string(s_grid) ~
            ")";
    }
} // end class GMOPatch


version(gmopatch_test) {
    import util.msg_service;
    int main() {
        auto p00 = Vector3([0.0, 0.1, 3.0]);
        auto p10 = Vector3(1.0, 0.1, 3.0);
        auto p11 = Vector3(1.0, 1.1, 3.0);
        auto p01 = Vector3(0.0, 1.1, 3.0);
        auto my_patch = new GMOPatch(p00, p10, p11, p01, 11, 11);
        auto c = my_patch(0.5, 0.5);
        // debug { import std.stdio; writeln("c=", c); }
        assert(approxEqualVectors(c, Vector3(0.477, 0.646, 3.0), 1.0e-1), failedUnitTest());
        c = my_patch(0.1, 0.1);
        // debug { import std.stdio; writeln("c=", c); }
        assert(approxEqualVectors(c, Vector3(0.117, 0.235, 3.0), 1.0e-1), failedUnitTest());
        return 0;
    }
}
