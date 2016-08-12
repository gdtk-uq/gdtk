/**
 * sketch_demo.d  Demonstrate rendering of some geometric elements.
 *
 * Author(s):
 *     Peter J.
 * Version:
 *     2016-08-10 adapted from svg_demo.d
 */

import std.stdio;
import geom;
import gpath;
import surface;
import sketch;

void main()
{
    writeln("Begin demonstration of the geometric rendering.");
    auto s = new Sketch("svg", "xyortho");
    s.canvas.set(0.0,0.0,120.0,120.0);
    s.viewport.set(-2.0,-2.0,2.0,2.0);
    auto a = Vector3(2.0, 0.0);
    auto b = Vector3(1.0, 1.0);
    auto c = Vector3(1.0, 0.0);
    auto abc = new Arc(a, b, c);
    s.start("test.svg");
    s.line(Vector3(-1.0,1.0,0.0), Vector3(1.0,-1.0,0.0));
    s.render(abc);
    s.rule("x", -1.2, 1.2, 0.4, Vector3(0.0,-1.3), 0.03, "%.1f", 0.12, 8);
    s.rule("y", -1.2, 1.2, 0.4, Vector3(-1.3,0.0), 0.03, "%.1f", 0.06, 8);
    s.setFillColour("green");
    auto p00 = Vector3(0.0, 0.1);
    auto p10 = Vector3(1.0, 0.1);
    auto p11 = Vector3(1.0, 1.1);
    auto p01 = Vector3(0.0, 1.1);
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    s.dotlabel(p00, "p00"); s.dotlabel(p10, "p10");
    s.dotlabel(p11, "p11"); s.dotlabel(p01, "p01");
    s.render(my_patch);
    s.finish();
    //
    s.start("test2.svg");
    s.viewport.set(0.0,0.0,120.0,120.0);
    s.line(Vector3(0.0,0.0), Vector3(90.0,120.0));
    s.setLineWidth(0.5);
    s.setFillColour("yellow");
    s.circle(Vector3(25.0,85.0), 12.3);
    s.setLineWidth(0.25);
    s.polyline([Vector3(0,50), Vector3(10,60), Vector3(20,50), Vector3(30,60)], true);
    s.text(Vector3(25.0,85.0), "Circle", Vector3(0.866,-0.5), Vector3(0.0,0.0,1.0),
	   "middle", 10);
    s.arc(Vector3(90.0,0.0), Vector3(60.0,30.0), Vector3(60.0,0.0));
    s.setLineWidth(0.75);
    s.arc(Vector3(30.0,0.0), Vector3(60.0,30.0), Vector3(60.0,0.0), true);
    s.bezier3(Vector3(80.0,80.0), Vector3(40.0,80.0),
	      Vector3(80.0,100.0), Vector3(40.0,100.0));
    s.setLineWidth(0.01);
    s.dotlabel(Vector3(70.0,20.0), "a");
    s.setLineWidth(0.1);
    s.polygon([Vector3(60,10), Vector3(95,10), Vector3(95,50), Vector3(60,50)],
	      true, true, true);
    s.finish();
    writeln("Done.");
}

