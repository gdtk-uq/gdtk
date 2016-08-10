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
import sketch;

void main()
{
    writeln("Begin demonstration of the geometric rendering.");
    auto s = new Sketch("svg", "xyortho");
    s.canvas.set(0.0,0.0,120.0,120.0);
    s.viewport.set(-1.0,-1.0,1.0,1.0);
    s.begin("test.svg");
    s.line(Vector3(-1.0,1.0,0.0), Vector3(1.0,-1.0,0.0));
    s.end();
    //
    s.begin("test2.svg");
    s.viewport.set(0.0,0.0,120.0,120.0);
    s.line(Vector3(0.0,0.0), Vector3(90.0,120.0));
    s.setLineWidth(0.5);
    s.circle(Vector3(25.0,85.0), 12.3);
    s.setLineWidth(0.25);
    s.polyline([Vector3(0.0,50), Vector3(10.0,60),
		Vector3(20.0,50), Vector3(30.0,60)], true);
    s.text(Vector3(25.0,85.0), "Circle", Vector3(0.866,-0.5), Vector3(0.0,0.0,1.0),
	   "middle", 10);
    s.arc(Vector3(90.0,0.0), Vector3(60.0,30.0), Vector3(60.0,0.0));
    s.setLineWidth(0.75);
    s.arc(Vector3(30.0,0.0), Vector3(60.0,30.0), Vector3(60.0,0.0), true);
    s.bezier3(Vector3(80.0,80.0), Vector3(40.0,80.0),
	      Vector3(80.0,100.0), Vector3(40.0,100.0));
    s.setLineWidth(0.01);
    s.dotlabel(Vector3(70.0,20.0), "a");
    s.end();
    writeln("Done.");
}

