// nozzleexpansionpatch.d

// This is Wilson Chan's specialized surface for nozzle expansion regions.
// It was extracted from the nozzle.input.template in the nenzfr application.
// 2014-nov-18: Peter J.
// 2020-Jan-22: Kyle Lynch, ported to Eilmer4

module geom.surface.nozzleexpansionpatch;

import std.conv;
import nm.complex;
import nm.number;
import std.math;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;

class NozzleExpansionPatch : ParametricSurface {
public:
    Path cA; // south edge
    Path cB; // north edge
    Path cC; // west edge
    Path cD; // east edge

    this(const Path cA, const Path cB, const Path cC, const Path cD)
    {
        this.cA = cA.dup(); // south
        this.cB = cB.dup(); // north
        this.cC = cC.dup(); // west
        this.cD = cD.dup(); // east
    }

    this(ref const(NozzleExpansionPatch) other)
    {
        cA = other.cA.dup();
        cB = other.cB.dup();
        cC = other.cC.dup();
        cD = other.cD.dup();
    }

    override NozzleExpansionPatch dup() const
    {
        return new NozzleExpansionPatch(this.cA, this.cB, this.cC, this.cD);
    }

    override Vector3 opCall(double r, double s) const
    {
        // Define the expansion_region of a nozzle using a specialised Surface Function.
        // For viscous simulations, it is necessary to keep the cells near
        // the non-slip walls as orthogonal to the walls as possible. However,
        // because the "AO" option in make_patch() does not give a grid that is
        // good enough for the nozzle geometry, a specialised surface function
        // has to be used. Points in the grid along the north, east and west
        // edges follow that specified by n_north, n_east and n_west. The rest
        // of the other points in the grid are built by creating strips of
        // quadratic Bezier curves that run from the nozzle wall to the axis.
        // The use of quadratic Bezier curves allows the generated points to
        // be orthogonal to the wall near the nozzle wall and orthogonal to
        // the axis near the axis.
        // 2012-Apr-20: Wilson Chan

        Vector3 p;
        p.refz = 0.0;
        if (r == 0.0) {

            p.refx = cB(r).x; // north
            p.refy = cC(s).y; // west

        } else if (r == 1.0) {

            p.refx = cB(r).x; // north
            p.refy = cD(s).y; // east

        } else if (s == 1.0) {

            p.refx = cB(r).x; // north
            p.refy = cB(r).y; // north

        } else {

            // Wall point (Bezier control point 1)
            auto wall_pt_x = cB(r).x;
            auto wall_pt_y = cB(r).y;

            // Angle perpendicular to the wall at wall point
            auto wall_angle = atan((cB(r+0.0001).y - cB(r-0.0001).y) /
                              (cB(r+0.0001).x - cB(r-0.0001).x));

            // If the expansion region starts sharply from the throat, then we
            // need some way of transitioning from the grid in the throat region
            // to the expansion region. To do so, we tweak the wall_angle in a
            // small starting region in the nozzle (say, 2% of the length of the
            // expansion region). The wall_angle starts from 0 degrees at the
            // start of this small region and smooths it out to the actual
            // wall_angle at the end of this small region.
            auto north_length = cB(1.0).x - cB(0.0).x;
            if ((cB(r).x - cB(0.0).x) <= (0.02 * north_length)) {
                wall_angle = (cB(r).x - cB(0.0).x) / (0.02 * north_length) * wall_angle;
            }

            // Do the same for a small region at the end of the nozzle. This is
            // to accommodate to nozzles that have been a non-zero gradient at
            // at the nozzle exit (which is brought about by nozzle truncation).
            if ((cB(r).x - cB(0.0).x) >= (0.98 * north_length)) {
                wall_angle = (cB(1.0).x - cB(r).x) / (0.02 * north_length) * wall_angle;
            }

            // Mid point (Bezier control point 2).
            auto mid_pt_y = cB(r).y / 2.0;
            auto mid_pt_x = wall_pt_x + (mid_pt_y * tan(wall_angle));

            // Axis point (Bezier control point 3).
            auto axis_pt_x = mid_pt_x;
            auto axis_pt_y = 0.0;

            // Generate t for quadratic Bezier curve equation.
            auto t = (1.0 - s);

            // Generate point on quadratic Bezier curve.
            p.refx = (1-t)*(1-t)*wall_pt_x + 2*t*(1-t)*mid_pt_x + t*t*axis_pt_x;
            p.refy = (1-t)*(1-t)*wall_pt_y + 2*t*(1-t)*mid_pt_y + t*t*axis_pt_y;
        }

        return p;
    }

    override string toString() const
    {
        return "Nozzleexpansionpatch(cA=" ~ to!string(cA) ~ ", cB=" ~ to!string(cB) ~
            ", cC=" ~ to!string(cC) ~ ", cD=" ~ to!string(cD) ~ ")";
    }

} // end class NozzleExpansionPatch

// [TODO] version(nozzleexpansionpatch_test)
