/**
 * svg_demo.d  Demonstrate rendering of some Scalable Vector Graphics elements.
 *
 * Author(s):
 *     Peter J.
 * Version:
 *     2014-06-16
 */

import std.stdio;
import misc.svg;

void main()
{
    writeln("Begin demonstration of the SVG rendering primitives.");
    auto s = new SVGContext(100, 120.0, "mm", "test run");

    s.open("test.svg");
    s.line(0.0, 0.0, 100.0, 120.0);
    s.close();

    s.open("test2.svg");
    s.line(0.0, 0.0, 90.0, 120.0);
    s.setLineWidth(0.5);
    s.setFillColour("yellow");
    s.circle(25.0, 85.0, 12.3);
    s.setLineWidth(0.25);
    s.polyline([0.0, 10.0, 20.0, 30.0], [50.0, 60.0, 50.0, 60.0], true);
    s.text(25.0, 85.0, "Circle", -30.0, "middle", 10);
    s.arc(90.0, 0.0, 60.0, 30.0, 60.0, 0.0);
    s.setLineWidth(0.75);
    s.arc(30.0, 0.0, 60.0, 30.0, 60.0, 0.0, true);
    s.bezier3(80.0, 80.0, 40.0, 80.0, 80.0, 100.0, 40.0, 100.0);
    s.setLineWidth(0.01);
    s.dotlabel(70.0, 20.0, "a");
    s.setLineWidth(0.1);
    s.polygon([60.0, 95.0, 95.0, 60.0],[10.0, 10.0, 50.0, 50.0], true, true, true);
    s.close();

    writeln("Done.");
}

