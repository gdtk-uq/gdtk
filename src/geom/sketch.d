/**
 * sketch.d  Sketch geometric elements using a selected renderer.
 *
 * These are the user-callable functions for sketching geometric elements
 * while preparing the data for a flow simulation.  They are modelled on 
 * the capabilities of the SVG renderer but work in the geometric space
 * of the flow simulation.
 *
 * Author(s):
 *     Peter J.
 * Version:
 *     2016-08-10 Let's start over with just the 2D SVG renderer
 *     but keep in mind that we'll want to be able to use Tcl/Tk,
 *     Cairo and OpenGL, eventually.
 *
 */

module sketch;

import std.stdio;
import std.math;
import std.string;
import std.algorithm;

import svg;
import geom;

struct Extents{
    double x0, x1, y0, y1;
    double width, height;
    void set(double x0, double y0, double x1, double y1)
    {
	this.x0 = x0; this.y0 = y0;
	this.x1 = x1; this.y1 = y1;
	width = x1 - x0;
	height = y1 - y0;
    } 
}

string[] renderers = ["svg", "tcltk", "cairo", "opengl"];
string[] projections = ["xyortho", "isometric", "oblique"];

class Sketch {
public:
    // We have some state and configuration data.
    string renderer_name;
    SVGContext svg;
    string projection_name;
    string title;
    string description;
    // Rendering canvas and viewport into user-space
    // 
    //     (x0,y1)---------(x1,y1)
    //        |               |
    //        |               |
    //        |               |
    //     (x0,y0)---------(x1,y0)
    //
    Extents canvas; // canvas-space in mm or pixels, depending on renderer
    Extents viewport; // user-space in metres

public:
    this(string renderer_name="svg", string projection_name="xyortho")
    {
	if (!find(renderers, renderer_name).empty) { 
	    this.renderer_name = renderer_name;
	} else {
	    throw new Exception("Unknown renderer: " ~ renderer_name);
	}
	if (!find(projections, projection_name).empty) {
	    this.projection_name = projection_name;
	} else {
	    throw new Exception("Unknown projection: " ~ projection_name);
	}
	assert(this.renderer_name == "svg", "other renderers unimplemented"); // FIX-ME
	assert(this.projection_name == "xyortho", "other projections unimplemented"); // FIX-ME
	canvas.set(0.0, 0.0, 120.0, 120.0); // reasonable size for a drawing in SVG
	viewport.set(0.0, 0.0, 1.0, 1.0); // unit space because we don't yet know better
    }

    void begin(string file_name="sketch.svg")
    {
	switch(renderer_name) {
	case "svg":
	    svg = new SVGContext(canvas.width, canvas.height, "mm", title, description);
	    svg.open(file_name);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
    }

    void end()
    {
	switch(renderer_name) {
	case "svg":
	    svg.close();
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
    }

    void setLineWidth(double w)
    // Sets line width, in units appropriate to the renderer.
    {
	switch(renderer_name) {
	case "svg":
	    svg.setLineWidth(w); // in mm
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void setLineColour(string colour)
    {
	switch(renderer_name) {
	case "svg":
	    svg.setLineColour(colour);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void setFillColour(string colour)
    {
	switch(renderer_name) {
	case "svg":
	    svg.setFillColour(colour);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void clearFillColour()
    {
	switch(renderer_name) {
	case "svg":
	    svg.clearFillColour();
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void setDashArray(double dashLength=2.0, double gapLength=2.0)
    // Sets length of dashes and gaps, in units appropriate to the renderer.
    {
	switch(renderer_name) {
	case "svg":
	    svg.setDashArray(dashLength, gapLength); // in mm
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    double toCanvasX(double x)
    // Map from user-space to canvas-space.
    {
	return canvas.x0 + (x - viewport.x0)/(viewport.width)*canvas.width;
    }

    double toCanvasY(double y)
    // Map from user-space to canvas-space.
    {
	return canvas.y0 + (y - viewport.y0)/(viewport.height)*canvas.height;
    }

    void line(const Vector3 p0, const Vector3 p1, bool dashed=false)
    // Render a line from point 1 to point 2.
    {
	auto x0 = toCanvasX(p0.x); auto y0 = toCanvasY(p0.y);
	auto x1 = toCanvasX(p1.x); auto y1 = toCanvasY(p1.y);
	switch(renderer_name) {
	case "svg":
	    svg.line(x0, y0, x1, y1, dashed);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void polyline(const Vector3[] pa, bool dashed=false)
    {
	double[] xlist, ylist;
	foreach(p; pa) {
	    xlist ~= toCanvasX(p.x); ylist ~= toCanvasY(p.y);
	}
	switch(renderer_name) {
	case "svg":
	    svg.polyline(xlist, ylist, dashed);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void polygon(const Vector3[] pa, bool fill=true, bool stroke=true, bool dashed=false)
    {
	double[] xlist, ylist;
	foreach(p; pa) {
	    xlist ~= toCanvasX(p.x); ylist ~= toCanvasY(p.y);
	}
	switch(renderer_name) {
	case "svg":
	    svg.polygon(xlist, ylist, fill, stroke, dashed);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void arc(const Vector3 p0, const Vector3 p1, const Vector3 pc,
	     bool dashed=false)
    {
	double x0 = toCanvasX(p0.x); double y0 = toCanvasY(p0.y);
	double x1 = toCanvasX(p1.x); double y1 = toCanvasY(p1.y);
	double xc = toCanvasX(pc.x); double yc = toCanvasY(pc.y);
	switch(renderer_name) {
	case "svg":
	    svg.arc(x0, y0, x1, y1, xc, yc, dashed);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void circle(const Vector3 pc, double r, bool fill=true,
		bool stroke=true, bool dashed=false)
    {
	double xc = toCanvasX(pc.x); double yc = toCanvasY(pc.y);
	switch(renderer_name) {
	case "svg":
	    svg.circle(xc, yc, r, fill, stroke, dashed);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void bezier3(const Vector3 p0, const Vector3 p1,
		 const Vector3 p2, const Vector3 p3,
		 bool dashed=false)
    {
	double x0 = toCanvasX(p0.x); double y0 = toCanvasY(p0.y);
	double x1 = toCanvasX(p1.x); double y1 = toCanvasY(p1.y);
	double x2 = toCanvasX(p2.x); double y2 = toCanvasY(p2.y);
	double x3 = toCanvasX(p3.x); double y3 = toCanvasY(p3.y);
	switch(renderer_name) {
	case "svg":
	    svg.bezier3(x0, y0, x1, y1, x2, y2, x3, y3, dashed);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void text(const Vector3 p, string textString,
	      const Vector3 direction=Vector3(1.0,0.0,0.0), // x-direction
	      const Vector3 normal=Vector3(0.0,0.0,1.0), // z-direction
	      string anchor="start",
	      int fontSize=10, string colour="black", string fontFamily="sanserif")
    {
	double xp = toCanvasX(p.x); double yp = toCanvasY(p.y);
	double angle = atan2(direction.y, direction.x)*180.0/PI;
	switch(renderer_name) {
	case "svg":
	    svg.text(xp, yp, textString, angle, anchor, fontSize, colour, fontFamily);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void dotlabel(const Vector3 p, string label="", 
		  string anchor="middle", double dotSize=2.0,
		  int fontSize=10, string colour="black", string fontFamily="sanserif")
    {
	double xp = toCanvasX(p.x); double yp = toCanvasY(p.y);
	switch(renderer_name) {
	case "svg":
	    svg.dotlabel(xp, yp, label, anchor, dotSize, fontSize, colour, fontFamily);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }
 } // end class Sketch
