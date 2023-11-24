// nozzleexpansionpatch.d
//
// This is Wilson Chan's specialized surface for nozzle expansion regions.
// It was extracted from the nozzle.input.template in the nenzfr application.
//
// Define the expansion_region of a nozzle using a specialised Surface Function.
// For viscous simulations, it is necessary to keep the cells near the non-slip
// walls as orthogonal to the walls as possible.
// However, because the "AO" option in make_patch() does not give a grid that is
// good enough for the nozzle geometry, a specialised surface has to be used.
// The controlling input is the north edge.  The south edge is the x-axis and
// the west and east edges are assumed straight and normal to the south edge.
// The interior of the surface is interpolated by creating a quadratic path
// from the south-edge point up to the north-edge point.
// The use of a quadratic Bezier curve allows the generated points to be
// orthogonal to the wall near the nozzle wall and also to the x-axis.
//
// 2012-Apr-20: Wilson Chan
// 2014-nov-18: Peter J.
// 2020-Jan-22: Kyle Lynch, ported to Eilmer4
// 2020-May-09: PJ refactor to make definition clearer.

module geom.surface.nozzleexpansionpatch;

import std.conv;
import ntypes.complex;
import nm.number;
import std.math;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;

class NozzleExpansionPatch : ParametricSurface {
public:
    Path cNorth; // the north edge is the only user-supplied data

    this(const Path cNorth)
    {
        this.cNorth = cNorth.dup();
    }

    this(ref const(NozzleExpansionPatch) other)
    {
        cNorth = other.cNorth.dup();
    }

    override NozzleExpansionPatch dup() const
    {
        return new NozzleExpansionPatch(this.cNorth);
    }

    override Vector3 opCall(double r, double s) const
    {
        Vector3 p;
        p.z = 0.0;
        if (r == 0.0) {
            // East edge
            p.x = cNorth(0.0).x;
            p.y = s*cNorth(0.0).y;
        } else if (r == 1.0) {
            // West edge
            p.x = cNorth(1.0).x;
            p.y = s*cNorth(1.0).y;
        } else if (s == 1.0) {
            // North edge
            p.x = cNorth(r).x;
            p.y = cNorth(r).y;
        } else {
            // Any other point, including the south edge.
            //
            auto pNW = cNorth(0.0);
            auto pNE = cNorth(1.0);
            auto north_length = pNE.x - pNW.x;
            // Wall point (Bezier control point at s=1)
            auto wall_pt_x = cNorth(r).x;
            auto wall_pt_y = cNorth(r).y;
            // Angle perpendicular to the wall at wall point
            double eps = 0.0001;
            number wall_angle = 0.0;
            if (r < eps) {
                // Near west boundary, use one-sided forward difference.
                wall_angle = atan((cNorth(r+eps).y - cNorth(r).y) /
                                  (cNorth(r+eps).x - cNorth(r).x));
            } else if (1.0-r < eps) {
                // Near east bounary, use one-sided backward difference.
                wall_angle = atan((cNorth(r).y - cNorth(r-eps).y) /
                                  (cNorth(r).x - cNorth(r-eps).x));
            } else {
                // Away from either boundary, use central difference.
                wall_angle = atan((cNorth(r+eps).y - cNorth(r-eps).y) /
                                  (cNorth(r+eps).x - cNorth(r-eps).x));
            }
            // If the expansion region starts sharply from the throat, then we
            // need some way of transitioning from the grid in the throat region
            // to the expansion region. To do so, we tweak the wall_angle in a
            // small starting region in the nozzle (say, 2% of the length of the
            // expansion region). The wall_angle starts from 0 degrees at the
            // start of this small region and smooths it out to the actual
            // wall_angle at the end of this small region.
            if ((cNorth(r).x - pNW.x) <= (0.02 * north_length)) {
                auto scale = (cNorth(r).x - pNW.x) / (0.02 * north_length);
                wall_angle *= scale;
            }
            // Do the same for a small region at the east end of the nozzle.
            // This is to accommodate to nozzles that have been a non-zero gradient
            // at the nozzle exit, which is often brought about by nozzle truncation.
            if ((cNorth(r).x - pNW.x) >= (0.98 * north_length)) {
                auto scale = (pNE.x - cNorth(r).x) / (0.02 * north_length);
                wall_angle *= scale;
            }
            // Mid point (Bezier control point at s=0.5).
            auto mid_pt_y = cNorth(r).y / 2.0;
            auto mid_pt_x = wall_pt_x + (mid_pt_y * tan(wall_angle));
            // Axis point (Bezier control point at s=0).
            auto axis_pt_x = mid_pt_x;
            auto axis_pt_y = 0.0;
            if (s == 0.0) {
                p.x = axis_pt_x;
                p.y = axis_pt_y;
            } else {
                // Generate point on quadratic Bezier curve.
                p.x = (1-s)*(1-s)*axis_pt_x + 2*s*(1-s)*mid_pt_x + s*s*wall_pt_x;
                p.y = (1-s)*(1-s)*axis_pt_y + 2*s*(1-s)*mid_pt_y + s*s*wall_pt_y;
            }
        }
        return p;
    } // end opCall()

    override string toString() const
    {
        return "NozzleExpansionPatch(cNorth=" ~ to!string(cNorth) ~ ")";
    }

} // end class NozzleExpansionPatch

// [TODO] version(nozzleexpansionpatch_test)
