/**
 * gpath_demo.d Demonstrate some of the behaviour of the path primitives.
 *
 * Author: Peter J.
 * Version: 2015-02-19, 2015-04-21 Arc, Bezier
 */

import std.stdio;
import geom;
import gpath;

void main()
{
    writeln("Begin demonstration of the geometric Path primitives.");
    auto a = Vector3([1.0, 2.2, 3.0]);
    auto b = Vector3(1.0);
    writeln("a = ", a, ", b = ", b);
    auto ab = new Line(a, b);
    writeln("ab= ", ab);
    auto c = ab(0.5);
    writeln("ab(0.5)= ", c);

    writeln("Arc");
    a = Vector3([2.0, 2.0, 0.0]);
    b = Vector3([1.0, 2.0, 1.0]);
    c = Vector3([1.0, 2.0, 0.0]);
    auto abc = new Arc(a, b, c);
    writeln("abc= ", abc);
    auto d = abc(0.5);
    writeln("abc(0.5)= ", d);

    writeln("Arc3");
    auto adb = new Arc3(a, d, b);
    writeln("adb(0.5)= ", adb(0.5));

    writeln("Bezier");
    auto bez_adb = new Bezier([a, d, b]);
    writeln("Bezier adb= ", bez_adb);
    auto e = bez_adb(0.5);
    writeln("bez_adb(0.5)=", e);

    writeln("Polyline");
    auto polyline = new Polyline([abc, new Line(b, c)]);
    writeln("polyline= ", polyline);
    writeln("polyline(0.25)= ", polyline(0.25));
    writeln("polyline(0.50)= ", polyline(0.50));
    writeln("polyline(0.75)= ", polyline(0.75));

    writeln("ArcLengthParameterizedPath");
    auto p0 = Vector3([0.0, 0.0, 0.0]);
    auto p1 = Vector3([1.0, 1.0, 1.0]);
    auto p2 = Vector3([4.0, 4.0, 4.0]);
    auto alpp = new ArcLengthParameterizedPath(new Bezier([p0,p1,p2]));
    writeln("alpp= ", alpp);
    writeln("alpp(0.5)= ", alpp(0.5));

    writeln("SubRangedPath");
    auto srp = new SubRangedPath(polyline, 1.0, 0.0); // effectively reversed
    writeln("srp= ", srp);
    writeln("srp(0.25)= ", srp(0.25));
    writeln("srp(0.50)= ", srp(0.50));
    writeln("srp(0.75)= ", srp(0.75));

    writeln("ReversedPath");
    auto rp = new ReversedPath(polyline);
    writeln("rp= ", rp);
    writeln("rp(0.25)= ", rp(0.25));
    writeln("rp(0.50)= ", rp(0.50));
    writeln("rp(0.75)= ", rp(0.75));

    writeln("Spline (Polyline)");
    // A rough circle.
    auto pnts = [Vector3([0.0, -1.0, 0.0]),
		 Vector3([-1.0, 0.0, 0.0]),
		 Vector3([0.0, 1.0, 0.0]),
		 Vector3([1.0, 0.0, 0.0]),
		 Vector3([0.0, -1.0, 0.0])];
    auto circle = new Polyline(pnts);
    writeln("circle= ", circle);
    writeln("circle(5.0/8)= ", circle(5.0/8));  // approx 45 degrees
    
    writeln("Done gpath_demo.");
}
