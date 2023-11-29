// polyline.d
// Build Paths out of segments.
// This includes the building of splines with Bezier3 segments.
//
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.polyline;

import std.conv;
import std.math;
import std.stdio: File;
import std.string;

import nm.bbla;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;
import geom.gpath.line;
import geom.gpath.arc;
import geom.gpath.bezier;


class Polyline : Path {
public:
    Path[] segments; // collection of Path segments
    double[] t_values; // collection of segment break-points (in parameter t)
    bool closed = false; // assume an open path, unless told otherwise

    // Construct from the Path segments.
    this(in Path[] segments, bool isClosed=false, double tolerance=1.0e-10)
    {
        if (segments.length == 0) {
            throw new Error(text("Polyline() No segments present."));
        }
        foreach (myseg; segments) { this.segments ~= myseg.dup(); }
        if (isClosed) {
            // Test that the end-points of the path are close enough and,
            // if not, add a line segment to close the path.
            Vector3 p0 = segments[0](0.0);
            Vector3 p1 = segments[$-1](1.0);
            if (distance_between(p1, p0) > tolerance) {
                this.segments ~= new Line(p1, p0);
            }
        }
        closed = isClosed;
        t_values.length = this.segments.length;
        reset_breakpoints();
    }

    // Construct as a spline through specified points.
    this(const Vector3[] p_orig, bool isClosed=false, double tolerance=1.0e-10)
    {
        Vector3[] p; foreach(pnt; p_orig) { p ~= pnt; }
        if (isClosed && distance_between(p[0], p[$-1]) > tolerance) {
            // Add one more point to close the path.
            p ~= p[0];
        }
        auto m = p.length - 1;
        // Given m+1 interpolation points p, determine the m-segment
        // Bezier polyline that interpolates these points as a spline.
        // This is done by first determining the array of weight points
        // which define the spline and then evaluating the cubic
        // Bezier segments.
        // Reference:
        //     G. Engelin & F. Uhlig (1996)
        //     Numerical Algorithms with C
        //     Springer, Berlin
        //     Section 12.3.1

        Vector3[] d; d.length = m+1;  // weight points
        // For a natural spline, the first and last weight points
        // are also the first and last interpolation points.
        d[0] = p[0];
        d[m] = p[m];

        // For the initial guess at the remaining weight points,
        // just use the supplied data points.
        foreach (i; 1 .. m) { d[i] = p[i]; }
        // Apply Gauss-Seidel iteration until
        // the internal weight points converge.
        Vector3 old_p;
        number max_diff;
        foreach (j; 1 .. 50) {
            max_diff = 0.0;
            foreach (i; 1 .. m) {
                old_p = d[i];
                d[i] = 0.25 * (6.0*p[i] - d[i-1] - d[i+1]);
                Vector3 diff = d[i] - old_p;
                max_diff = fmax(max_diff, geom.abs(diff));
            } // end foreach i
            if ( max_diff < tolerance ) break;
        } // end foreach j

        // Final stage; calculate the Bezier segments
        Vector3[4] p03;
        Path[] seg;
        foreach (i; 0 ..  m) {
            p03[0] = p[i];
            p03[1] = (2.0*d[i] + d[i+1]) / 3.0;
            p03[2] = (d[i] + 2.0*d[i+1]) / 3.0;
            p03[3] = p[i+1];
            seg ~= new Bezier(p03);
        }
        // and pack them away.
        this(seg, isClosed);
    } // end spline constructor

    // Contructs a spline from a file containing x(,y(,z)) coordinates.
    this(string fileName, bool isClosed=false)
    {
        // This function takes a filename and processes it assuming that each
        // line contains (x,y,z) triples (space-delimited).  If any y- or z-values are
        // missing on a given line, they are assumed to be 0.0.  The x,y,z-triples
        // are gathered and used to create the Spline.
        // Ported Python code Spline2 from libgeom2.i 2015-10-05 by PJ
        Vector3[] points;
        auto f = File(fileName, "r");
        foreach (line; f.byLine) {
            auto tokens = line.strip().split();
            if (tokens.length == 0) continue; // ignore blank lines
            if (tokens[0] == "#") continue; // ignore comment lines
            number x = to!double(tokens[0]);
            number y = 0.0; if (tokens.length > 1) { y = to!double(tokens[1]); }
            number z = 0.0; if (tokens.length > 2) { z = to!double(tokens[2]); }
            points ~= Vector3(x, y, z);
        }
        this(points, isClosed);
    } // end spline constructor

    this(ref const(Polyline) other, bool isClosed=false)
    {
        this(other.segments, isClosed);
    }

    override Polyline dup() const
    {
        return new Polyline(segments, closed);
    }

    override Vector3 opCall(double t) const
    {
        // Evaluate B(t) without considering arc_length parameterization flag
        // or subrange.
        auto n = segments.length;
        if ( n == 1 ) return segments[0](t);
        size_t i;
        for ( i = 0; i < n; ++i ) {
            if ( t <= t_values[i] ) break;
        }
        if ( i >= n ) i = n - 1;  // last segment
        // At this point, t_values[i-1] < t <= t_values[i] (we hope)
        // Have assumed that the t breakpoints are well behaved.
        // i.e. There are no zero-length segments.
        double t_local;
        if ( i == 0 ) {
            t_local = t / t_values[i];
        } else {
            t_local = (t - t_values[i-1]) / (t_values[i] - t_values[i-1]);
        }
        return segments[i](t_local);
    } // end opCall()

    override string toString() const
    {
        return "Polyline(segments=" ~ to!string(segments) ~
            ", closed=" ~ to!string(closed) ~ ")";
    }
    override string classString() const
    {
        return "Polyline";
    }

    string toGmshString(ref int pointTag, ref int curveTag, ref int loopTag,
                        string label="", double len=1.0e-2)
    // We use int tags for the points, curves and loops written to the string.
    // Values passed in will be used as the starting values for these tags.
    // Subsequent calls should be careful to not double-up on tag values.
    {
        string str = "// "~label~"\n";
        str ~= format("len = %g;\n", len);
        int[] segmentTags;
        Vector3 startPoint = segments[0](0.0);
        Vector3 previousPoint = startPoint;
        pointTag += 1;
        str ~= format("Point(%d) = {%g, %g, %g, len};\n", pointTag,
                      startPoint.x, startPoint.y, startPoint.z);
        int startPointTag = pointTag;
        double tol = 1.0e-6;

        foreach (i, seg; segments) {
            switch (seg.classString()) {
            case "Line":
                // Put down the first point of the segment
                // only if it significantly different to the previous point.
                Vector3 p0 = seg(0.0);
                if (distance_between(p0, previousPoint) > tol) {
                    pointTag += 1;
                    str ~= format("Point(%d) = {%g, %g, %g, len};\n", pointTag, p0.x, p0.y, p0.z);
                    previousPoint.set(p0);
                }
                int tag0 = pointTag;
                int tag1;
                if ((i < segments.length-1) || (!closed)) {
                    Vector3 p1 = seg(1.0);
                    if (distance_between(p1, previousPoint) > tol) {
                        pointTag += 1;
                        str ~= format("Point(%d) = {%g, %g, %g, len};\n", pointTag, p1.x, p1.y, p1.z);
                        tag1 = pointTag;
                        previousPoint.set(p1);
                    }
                } else {
                    tag1 = startPointTag;
                }
                curveTag += 1;
                str ~= format("Line(%d) = {%d, %d};\n", curveTag, tag0, tag1);
                segmentTags ~= curveTag;
                break;
            case "Arc":
                str ~= "Arc not yet implemented\n";
                break;
            case "Bezier":
                // Put down the first point of the segment
                // only if it significantly different to the previous point.
                auto mySeg = cast(Bezier) seg;
                Vector3[] p; p ~= mySeg.B[0];
                if (distance_between(p[0], previousPoint) > tol) {
                    pointTag += 1;
                    str ~= format("Point(%d) = {%g, %g, %g, len};\n", pointTag, p[0].x, p[0].y, p[0].z);
                    previousPoint.set(p[0]);
                }
                int[] tags; tags ~= pointTag;
                // Put down intermediate control points.
                foreach (j; 1 .. mySeg.B.length-1) {
                    pointTag += 1;
                    str ~= format("Point(%d) = {%g, %g, %g, len};\n",
                                  pointTag, mySeg.B[j].x, mySeg.B[j].y, mySeg.B[j].z);
                    tags ~= pointTag;
                }
                // Conditionally put down final point of Bezier.
                if ((i < segments.length-1) || (!closed)) {
                    p ~= mySeg.B[$-1];
                    if (distance_between(p[$-1], previousPoint) > tol) {
                        pointTag += 1;
                        str ~= format("Point(%d) = {%g, %g, %g, len};\n",
                                      pointTag, p[$-1].x, p[$-1].y, p[$-1].z);
                        tags ~= pointTag;
                        previousPoint.set(p[$-1]);
                    }
                } else {
                    tags ~= startPointTag;
                }
                curveTag += 1;
                str ~= format("Bezier(%d) = {%d", curveTag, tags[0]);
                foreach (j; 1 .. tags.length-1) { str ~= format(", %d", tags[j]); }
                str ~= format(", %d};\n", tags[$-1]);
                segmentTags ~= curveTag;
                break;
            default:
                str ~= "// Segment type " ~ seg.classString() ~ " not handled.\n";
            } // end switch
        } // end foreach
        //
        loopTag += 1;
        str ~= format("Curve Loop(%d) = {", loopTag);
        foreach (i, segTag; segmentTags) {
            str ~= format("%d", segTag);
            str ~= (i+1 == segmentTags.length) ? "};\n" : ", ";
        }
        curveTag += 1;
        str ~= format("Physical Curve(\"%s\", %d) = {", label, curveTag);
        foreach (i, segTag; segmentTags) {
            str ~= format("%d", segTag);
            str ~= (i+1 == segmentTags.length) ? "};\n" : ", ";
        }
        return str;
    }

private:
    void reset_breakpoints()
    {
        // Set up the parameter breakpoints based on cumulative length.
        t_values[0] = segments[0].length().re;
        foreach (i; 1 .. segments.length) {
            t_values[i] = t_values[i-1] + segments[i].length().re;
        }
        double L_total = t_values[$-1];
        foreach (i; 0 .. segments.length) { t_values[i] /= L_total; }
    } // end reset_breakpoints()
} // end class Polyline


version(polyline_test) {
    import util.msg_service;
    int main() {
        auto a = Vector3([2.0, 2.0, 0.0]);
        auto b = Vector3([1.0, 2.0, 1.0]);
        auto c = Vector3([1.0, 2.0, 0.0]);
        auto abc = new Arc(a, b, c);
        auto polyline1 = new Polyline([abc, new Line(b, c)]);
        auto f = polyline1(0.5);
        assert(approxEqualVectors(f, Vector3(1.28154, 2, 0.95955)), failedUnitTest());
        auto polyline2 = new Polyline(polyline1, true);
        auto g = polyline2(0.999999);
        assert(approxEqualVectors(g, a), failedUnitTest());
                //
        version(complex_numbers) {
            // Try out the complex derivative evaluation.
            double h = 1.0e-20;
            number ih = complex(0,h);
            number zero = 0.0;
            number one = 1.0;
            double alpha = 1.0;
            auto pa_dash = Vector3(alpha+ih, alpha+ih);
            auto pb = Vector3(one, zero);
            auto pc = Vector3(zero, one);
            auto line0 = new Line(pc, pa_dash);
            auto line1 = new Line(pa_dash, pb);
            auto polyline3 = new Polyline([line0, line1]);
            // What we want to compute is the sensitivity
            // of the midpoint of the arc with respect to alpha.
            double dpmid_da_x = polyline3(0.5).x.im / h;
            double dpmid_da_y = polyline3(0.5).y.im / h;
            double dpmid_da_z = polyline3(0.5).z.im / h;
            assert(isClose(dpmid_da_x,1.0), failedUnitTest());
            assert(isClose(dpmid_da_y,1.0), failedUnitTest());
            assert(isClose(dpmid_da_z,0.0), failedUnitTest());
        }
        return 0;
    }
} // end polyline_test
