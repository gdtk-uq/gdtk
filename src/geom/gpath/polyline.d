// polyline.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.polyline;

import std.conv;
import std.math;
import std.stdio: File;
import std.string;

import nm.bbla;
import nm.complex;
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

    // Construct from the Path segments.
    this(in Path[] segments)
    {
        if (segments.length == 0) {
            throw new Error(text("Polyline() No segments present."));
        }
        foreach (myseg; segments) this.segments ~= myseg.dup(); 
        t_values.length = this.segments.length;
        reset_breakpoints();
    }
    
    // Construct as a spline through specified points.
    this(const Vector3[] p, double tolerance=1.0e-10)
    {
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
                d[i] = to!number(0.25) * (to!number(6.0) * p[i] - d[i-1] - d[i+1]);
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
            p03[1] = (to!number(2.0) * d[i] + d[i+1]) / to!number(3.0);
            p03[2] = (d[i] + to!number(2.0) * d[i+1]) / to!number(3.0);
            p03[3] = p[i+1];
            seg ~= new Bezier(p03);
        }
        // and pack them away.
        this(seg);
    } // end spline constructor

    // Contructs a spline from a file containing x(,y(,z)) coordinates.
    this(string fileName)
    {
        // This function takes a filename and processes it assuming that each
        // line contains (x,y,z) triples (space-delimited).  If any values are
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
        this(points);
    } // end spline constructor
    
    this(ref const(Polyline) other)
    {
        this(other.segments);
    }
    
    override Polyline dup() const
    {
        return new Polyline(segments);
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
        return "Polyline(segments=" ~ to!string(segments) ~ ")";
    }
    override string classString() const
    {
        return "Polyline";
    }

private:
    void reset_breakpoints()
    {
        // Set up the parameter breakpoints based on cumulative length.
        t_values[0] = segments[0].length();
        foreach (i; 1 .. segments.length) t_values[i] = t_values[i-1] + segments[i].length(); 
        double L_total = t_values[$-1];
        foreach (i; 0 .. segments.length) t_values[i] /= L_total; 
    } // end reset_breakpoints()
} // end class Polyline


version(polyline_test) {
    import util.msg_service;
    int main() {
        auto a = Vector3([2.0, 2.0, 0.0]);
        auto b = Vector3([1.0, 2.0, 1.0]);
        auto c = Vector3([1.0, 2.0, 0.0]);
        auto abc = new Arc(a, b, c);
        auto polyline = new Polyline([abc, new Line(b, c)]);
        auto f = polyline(0.5);
        assert(approxEqualVectors(f, Vector3(1.28154, 2, 0.95955)), failedUnitTest());
        return 0;
    }
} // end polyline_test
