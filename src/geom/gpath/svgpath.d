// svgpath.d
// An SVG-like path built as a Polyline Path.
// PJ 2019-04-06

module geom.gpath.svgpath;

import std.conv;
import std.math;
import std.stdio: File;
import std.string;
import std.algorithm;

import nm.bbla;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;
import geom.gpath.line;
import geom.gpath.arc;
import geom.gpath.bezier;
import geom.gpath.polyline;


class SVGPath : Polyline {
public:
    this(string txt,
         double xtrans=0.0, double ytrans=0.0,
         double xscale=1.0, double yscale=1.0,
         double tolerance=1.0e-10)
    {
        cp.clear(); // Start with current-point at origin.
        p0.clear();
        this.xtrans = xtrans; this.ytrans = ytrans;
        this.xscale = xscale; this.yscale = yscale;
        this.tolerance = tolerance;
        interpret_svg_path(txt);
        super(segments, closed, tolerance);
    }

    this(ref const(SVGPath) other)
    {
        this(other.txt,
             other.xtrans, other.ytrans,
             other.xscale, other.yscale,
             other.tolerance);
    }

    override SVGPath dup() const
    {
        return new SVGPath(txt, xtrans, ytrans, xscale, yscale, tolerance);
    }

    override string toString() const
    {
        return "SVGPath(path=\""~ txt~"\"" ~
            format(", xtrans=%g, ytrans=%g", xtrans, ytrans) ~
            format(", xscale=%g, yscale=%g", xscale, yscale) ~
            format(", tolerance=%g)", tolerance);
    }

    override string classString() const
    {
        return "SVGPath";
    }

private:
    string txt; // Keep a copy of the original text specification.
                // It will be useful for making copies of the path.
    Vector3 p0; // Starting point.
    Vector3 cp; // Current point.
    // Apply the translation first, to absolute coordinates.
    double xtrans, ytrans;
    // Apply the scales second, to both absolute and relative coordinates.
    double xscale, yscale;
    // A place to store the generated gpath segments.
    Path[] segments;
    bool closed = false;
    double tolerance;

    void interpret_svg_path(string txt)
    {
        txt = txt.strip();
        foreach (cmdStr; txt.split(";")) {
            cmdStr = cmdStr.strip();
            char cmd = cmdStr[0];
            string argStr = "";
            if (cmdStr.length > 1) { argStr = cmdStr[1..$]; }
            do_svg_path_command(cmd, argStr);
        }
    } // end interpret_svg_path()

    void do_svg_path_command(char cmd, string argStr)
    {
        // We will interpret a subset of the SVG path commands
        // that seems to be enough to follow glyph outlines
        // as written by Inkscape.
        //
        // My reference has been:
        // J. David Eisenberg & Amelia Bellamy-Royds
        // SVG Essentials
        // 2nd Edition O'Reilly
        //
        argStr = argStr.strip();
        // import std.stdio: writeln;
        // writeln("cmd=", cmd, " argStr=", argStr);
        //
        import std.format: formattedRead;
        double x1, y1, x2, y2, x3, y3;
        Vector3 p1, p2, p3;
        switch (cmd) {
        case 'M': // Move absolute, resets starting point.
            formattedRead(argStr, "%g,%g", &x1, &y1);
            p0.set((x1+xtrans)*xscale, (y1+ytrans)*yscale);
            cp.set(p0);
            break;
        case 'm': // Move relative to current position, resets starting point.
            formattedRead(argStr, "%g,%g", &x1, &y1);
            p0.set(cp.x+x1*xscale, cp.y+y1*yscale);
            cp.set(p0);
            break;
        case 'L':
            formattedRead(argStr, "%g,%g", &x1, &y1);
            p1.set((x1+xtrans)*xscale, (y1+ytrans)*yscale);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'l':
            formattedRead(argStr, "%g,%g", &x1, &y1);
            p1.set(cp.x+x1*xscale, cp.y+y1*yscale);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'H':
            formattedRead(argStr, "%g", &x1);
            p1.set(to!number((x1+xtrans)*xscale), cp.y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'h':
            formattedRead(argStr, "%g", &x1);
            p1.set(cp.x+x1*xscale, cp.y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'V':
            formattedRead(argStr, "%g", &y1);
            p1.set(cp.x, to!number((y1+ytrans)*yscale));
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'v':
            formattedRead(argStr, "%g", &y1);
            p1.set(cp.x, cp.y+y1*yscale);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'Q':
            formattedRead(argStr, "%g,%g %g,%g", &x1, &y1, &x2, &y2);
            p1.set((x1+xtrans)*xscale, (y1+ytrans)*yscale);
            p2.set((x2+xtrans)*xscale, (y2+ytrans)*yscale);
            segments ~= new Bezier([cp, p1, p2]);
            cp.set(p2);
            break;
        case 'q':
            formattedRead(argStr, "%g,%g %g,%g", &x1, &y1, &x2, &y2);
            p1.set(cp.x+x1*xscale, cp.y+y1*yscale);
            p2.set(cp.x+x2*xscale, cp.y+y2*yscale);
            segments ~= new Bezier([cp, p1, p2]);
            cp.set(p2);
            break;
        case 'C':
            formattedRead(argStr, "%g,%g %g,%g %g,%g", &x1, &y1, &x2, &y2, &x3, &y3);
            p1.set((x1+xtrans)*xscale, (y1+ytrans)*yscale);
            p2.set((x2+xtrans)*xscale, (y2+ytrans)*yscale);
            p3.set((x3+xtrans)*xscale, (y3+ytrans)*yscale);
            segments ~= new Bezier([cp, p1, p2, p3]);
            cp.set(p3);
            break;
        case 'c':
            formattedRead(argStr, "%g,%g %g,%g %g,%g", &x1, &y1, &x2, &y2, &x3, &y3);
            p1.set(cp.x+x1*xscale, cp.y+y1*yscale);
            p2.set(cp.x+x2*xscale, cp.y+y2*yscale);
            p3.set(cp.x+x3*xscale, cp.y+y3*yscale);
            segments ~= new Bezier([cp, p1, p2, p3]);
            cp.set(p3);
            break;
        case 'Z':
            p1.set(p0);
            if (distance_between(cp, p1) > tolerance) {
                segments ~= new Line(cp, p1);
                cp.set(p1);
            }
            closed = true;
            break;
        default:
            throw new Exception("Unknown SVG command: " ~ cmd);
        } // end switch
    } // end do_svg_path_command()

} // end class SVGPath


version(svgpath_test) {
    import util.msg_service;
    int main()
    {
        // We are going to use semicolons as separators between path commands,
        // single-character command names (MmLlZ) and
        // commas as separators between coordinate values.
        // Note that the path is restricted to the z=0 plane.
        auto pth1 = new SVGPath("M3.0,3.0;L4.0,3.0;v1.0;h-1.0;Z");
        auto p1 = pth1(0.5);
        // import std.stdio: writeln;
        // writeln("pth1=", pth1);
        // writeln("pth1(0.5)=", p1);
        assert(approxEqualVectors(p1, Vector3(4.0,4.0,0.0)), failedUnitTest());
        //
        auto pth2 = new SVGPath("M2.0,0.0;Q2.0,2.0 0.0,2.0");
        auto p2 = pth2(0.5);
        // writeln("pth2=", pth2);
        // writeln("pth2(0.5)=", p2);
        assert(approxEqualVectors(p2, Vector3(1.5,1.5,0.0)), failedUnitTest());
        //
        // Approximate a quarter circle
        auto pth3 = new SVGPath("M2.0,0.0;C2.0,1.10457 1.10457,2.0 0.0,2.0");
        auto p3 = pth3(0.5);
        // writeln("pth3=", pth3);
        // writeln("pth3(0.5)=", p3);
        assert(approxEqualVectors(p3, Vector3(1.41421,1.41421,0.0)), failedUnitTest());
        return 0;
    }
} // end svgpath_test
