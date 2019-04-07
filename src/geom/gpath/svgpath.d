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
import nm.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;
import geom.gpath.line;
import geom.gpath.arc;
import geom.gpath.bezier;
import geom.gpath.polyline;


class SVGPath : Polyline {
public:
    this(string txt, double tolerance=1.0e-10)
    {
        cp.clear(); // Start with current-point at origin.
        p0.clear();
        this.tolerance = tolerance;
        interpret_svg_path(txt);
        super(segments, closed, tolerance);
    }
    
    this(ref const(SVGPath) other)
    {
        this(other.txt, other.tolerance);
    }
    
    override SVGPath dup() const
    {
        return new SVGPath(txt);
    }

    override string toString() const
    {
        return "SVGPath(path=\""~ txt~"\")";
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
        double x, y;
        Vector3 p1, p2, p3;
        switch (cmd) {
        case 'M': // Move absolute, resets starting point.
            formattedRead(argStr, "%g,%g", &x, &y);
            p0.set(x, y);
            cp.set(p0);
            break;
        case 'm': // Move relative to current position, resets starting point.
            formattedRead(argStr, "%g,%g", &x, &y);
            p0.set(cp.x+x, cp.y+y);
            cp.set(p0);
            break;
        case 'L':
            formattedRead(argStr, "%g,%g", &x, &y);
            p1.set(x, y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'l':
            formattedRead(argStr, "%g,%g", &x, &y);
            p1.set(cp.x+x, cp.y+y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'H':
            formattedRead(argStr, "%g", &x);
            p1.set(to!number(x), cp.y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'h':
            formattedRead(argStr, "%g", &x);
            p1.set(cp.x+x, cp.y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'V':
            formattedRead(argStr, "%g", &y);
            p1.set(cp.x, to!number(y));
            segments ~= new Line(cp, p1);
            cp.set(p1);
            break;
        case 'v':
            formattedRead(argStr, "%g", &y);
            p1.set(cp.x, cp.y+y);
            segments ~= new Line(cp, p1);
            cp.set(p1);
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
        auto p = pth1(0.5);
        // import std.stdio: writeln;
        // writeln("pth1=", pth1);
        // writeln("pth1(0.5)=", p);
        assert(approxEqualVectors(p, Vector3(4.0,4.0,0.0)), failedUnitTest());
        return 0;
    }
} // end svgpath_test
