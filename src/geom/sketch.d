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

import std.conv;
import std.stdio;
import std.math;
import std.string;
import std.algorithm;

import svg;
import geom;
import gpath;
import surface;
import volume;

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
    Extents canvas; // canvas- or device-space in mm or pixels, depending on renderer
    Extents viewport; // view coordinates in metres (following projection)
    double[4][4] proj_mat; // projection matrix, as per OpenGL description
    double[4][4] view_mat; // for when we want to view the model differently
    // We keep in mind the usual 4-component homogeneous coordinate representation
    // (x,y,z,w) so that we can accommodate perspective projection eventually.
    // Graphics pipeline processes are:
    // (1) 3D vertex data for model --> multiply by view_matrix --> 3D world coordinates
    // (2) 3D world coordinates --> multiply by projection_matrix --> 2D view coordinates
    // (3) 2D view coordinates --> viewport transform --> 2D window/canvas coordinates
    //
    // Default view is having the z-axis coming towards your eye, x-axis to the right
    // and y-axis in the up direction.
    Vector3 eye = Vector3(0.0, 0.0, 0.0);
    Vector3 centre = Vector3(0.0, 0.0, -1.0);
    Vector3 up = Vector3(0.0, 1.0, 0.0);

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
	set_to_identity(proj_mat);
	set_to_identity(view_mat); // consistent with default view
	switch (this.projection_name) {
	case "xyortho":
	    proj_mat[2][2] = 0.0; // x,y unchanged, z eliminated
	    break;
	case "isometric":
	    double cos30 = cos(30.0/180.0*PI);
	    double sin30 = sin(30.0/180.0*PI);
	    proj_mat[0][0] =  cos30; proj_mat[0][1] = 0.0; proj_mat[0][2] = -cos30; proj_mat[0][3] = 0.0;
	    proj_mat[1][0] = -sin30; proj_mat[1][1] = 1.0; proj_mat[1][2] = -sin30; proj_mat[1][3] = 0.0;
	    proj_mat[2][0] =    0.0; proj_mat[2][1] = 0.0; proj_mat[2][2] =    0.0; proj_mat[2][3] = 0.0;
	    proj_mat[3][0] =    0.0; proj_mat[3][1] = 0.0; proj_mat[3][2] =    0.0; proj_mat[3][3] = 1.0;
	    break;
	case "oblique":
	    double zfactor = 1.0/(2.0*sqrt(2.0));
	    proj_mat[0][0] =  1.0; proj_mat[0][1] = 0.0; proj_mat[0][2] = -zfactor; proj_mat[0][3] = 0.0;
	    proj_mat[1][0] =  0.0; proj_mat[1][1] = 1.0; proj_mat[1][2] = -zfactor; proj_mat[1][3] = 0.0;
	    proj_mat[2][0] =  0.0; proj_mat[2][1] = 0.0; proj_mat[2][2] =      0.0; proj_mat[2][3] = 0.0;
	    proj_mat[3][0] =  0.0; proj_mat[3][1] = 0.0; proj_mat[3][2] =      0.0; proj_mat[3][3] = 1.0;
	    break;
	default: throw new Exception("other projections unimplemented");
	}
	canvas.set(0.0, 0.0, 120.0, 120.0); // reasonable size for a drawing in SVG
	viewport.set(0.0, 0.0, 1.0, 1.0); // unit space because we don't yet know better
    } // end this

    override string toString()
    {
	string str = "Sketch(";
	str ~= "renderer=\""~renderer_name~"\"";
	str ~= ", projection=\""~projection_name~"\"";
	str ~= ", eye=" ~ to!string(eye);
	str ~= ", centre=" ~ to!string(centre);
	str ~= ", up=" ~ to!string(up);
	str ~= ")";
	return str;
    }

    void set_to_identity(ref double[4][4] mat)
    {
	foreach (i; 0 .. 4) {
	    foreach (j; 0 .. 4) { mat[i][j] = 0.0; }
	    mat[i][i] = 1.0;
	}
    }

    void matrix_vector_multiply(ref double[4][4] mat, ref double[4] vec)
    {
	double[4] result;
	foreach (i; 0 .. 4) {
	    result[i] = 0.0;
	    foreach (j; 0 .. 4) { result[i] += mat[i][j] * vec[j]; }
	}
	foreach (i; 0 .. 4) { vec[i] = result[i]; }
    }
    
    void apply_transform(string transform_name, ref Vector3 p, ref double w)
    {
	double[4] ph; // homogeneous coordinate representation
	ph[0] = p.x; ph[1] = p.y; ph[2] = p.z; ph[3] = w;
	switch (transform_name) {
	case "view":
	    matrix_vector_multiply(view_mat, ph);
	    break;
	case "projection":
	    matrix_vector_multiply(proj_mat, ph);
	    break;
	default:
	    // do nothing.
	}
	p.refx = ph[0]; p.refy = ph[1]; p.refz = ph[2]; w = ph[3];
    } // end apply_transform

    void look_at(const Vector3 eye, const Vector3 centre, const Vector3 up)
    // Set the view matrix so that we have a new view of the 3D model.
    // In the new coordinates, we are looking back along the znew-axis,
    // with ynew up and xnew toward the right. 
    {
	this.eye = eye;
	this.centre = centre;
	this.up = up;
	Vector3 znew = (eye - centre).normalize();
	Vector3 xnew = (cross(up, znew)).normalize();
	Vector3 ynew = (cross(znew, xnew)).normalize();
	view_mat[0][0] = xnew.x; view_mat[0][1] = xnew.y; view_mat[0][2] = xnew.z;
	view_mat[1][0] = ynew.x; view_mat[1][1] = ynew.y; view_mat[1][2] = ynew.z;
	view_mat[2][0] = znew.x; view_mat[2][1] = znew.y; view_mat[2][2] = znew.z;
	view_mat[0][3] = -(dot(centre,xnew));
	view_mat[1][3] = -(dot(centre,ynew));
	view_mat[2][3] = -(dot(centre,znew));
	view_mat[3][0] = 0.0; view_mat[3][1] = 0.0; view_mat[3][2] = 0.0; view_mat[3][3] = 1.0;
    } // end look_at
    
    void start(string file_name="sketch.svg")
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

    void finish()
    {
	switch(renderer_name) {
	case "svg":
	    svg.close();
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
    }

    void setLineWidth(double width)
    // Sets line width, in units appropriate to the renderer.
    {
	switch(renderer_name) {
	case "svg":
	    svg.setLineWidth(width); // in mm
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

    void begin_group(string id="", double opacity=1.0)
    // id needs to be a unique identifier string, if supplied
    // opacity is in range 0.0, 1.0
    {
	switch(renderer_name) {
	case "svg":
	    svg.begin_group(id, opacity);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void end_group()
    {
	switch(renderer_name) {
	case "svg":
	    svg.end_group();
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }

    void line(const Vector3 p0, const Vector3 p1, bool dashed=false)
    // Render a line from point 1 to point 2.
    {
	Vector3 p0tmp = Vector3(p0); Vector3 p1tmp = Vector3(p1);
	double w0 = 1.0; double w1 = 1.0;
	apply_transform("view", p0tmp, w0); apply_transform("view", p1tmp, w1);
	apply_transform("projection", p0tmp, w0); apply_transform("projection", p1tmp, w1);
	auto x0 = toCanvasX(p0tmp.x); auto y0 = toCanvasY(p0tmp.y);
	auto x1 = toCanvasX(p1tmp.x); auto y1 = toCanvasY(p1tmp.y);
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
	Vector3[] ptmp; double[] wtmp;
	foreach (p; pa) { ptmp ~= Vector3(p); wtmp ~= 1.0; }
	foreach (i; 0 .. ptmp.length) {
	    apply_transform("view", ptmp[i], wtmp[i]);
	    apply_transform("projection", ptmp[i], wtmp[i]);
	}	    
	double[] xlist, ylist;
	foreach(p; ptmp) {
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
	Vector3[] ptmp; double[] wtmp;
	foreach (p; pa) { ptmp ~= Vector3(p); wtmp ~= 1.0; }
	foreach (i; 0 .. ptmp.length) {
	    apply_transform("view", ptmp[i], wtmp[i]);
	    apply_transform("projection", ptmp[i], wtmp[i]);
	}	    
	double[] xlist, ylist;
	foreach(p; ptmp) {
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
	// This is a purely 2D xy-plane function.
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
	// This is a purely 2D xy-plane function.
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
	Vector3 p0tmp = Vector3(p0); Vector3 p1tmp = Vector3(p1);
	Vector3 p2tmp = Vector3(p0); Vector3 p3tmp = Vector3(p1);
	double w0 = 1.0; double w1 = 1.0;
	double w2 = 1.0; double w3 = 1.0;
	apply_transform("view", p0tmp, w0); apply_transform("view", p1tmp, w1);
	apply_transform("view", p2tmp, w2); apply_transform("view", p3tmp, w3);
	apply_transform("projection", p0tmp, w0); apply_transform("projection", p1tmp, w1);
	apply_transform("projection", p2tmp, w2); apply_transform("projection", p3tmp, w3);
	double x0 = toCanvasX(p0tmp.x); double y0 = toCanvasY(p0tmp.y);
	double x1 = toCanvasX(p1tmp.x); double y1 = toCanvasY(p1tmp.y);
	double x2 = toCanvasX(p2tmp.x); double y2 = toCanvasY(p2tmp.y);
	double x3 = toCanvasX(p3tmp.x); double y3 = toCanvasY(p3tmp.y);
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
	Vector3 p_tmp = Vector3(p);
	Vector3 direction_tmp = Vector3(direction);
	Vector3 normal_tmp = Vector3(normal);
	double w_p = 1.0; double w_dir = 1.0; double w_n = 1.0;
	apply_transform("view", p_tmp, w_p); apply_transform("projection", p_tmp, w_p);
	apply_transform("view", direction_tmp, w_dir); apply_transform("projection", direction_tmp, w_dir);
	apply_transform("view", normal_tmp, w_n); apply_transform("projection", normal_tmp, w_n);
	double xp = toCanvasX(p_tmp.x); double yp = toCanvasY(p_tmp.y);
	double angle = atan2(direction_tmp.y, direction_tmp.x)*180.0/PI;
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
	Vector3 p_tmp = Vector3(p); double w_p = 1.0;
	apply_transform("view", p_tmp, w_p); apply_transform("projection", p_tmp, w_p);
	double xp = toCanvasX(p_tmp.x); double yp = toCanvasY(p_tmp.y);
	switch(renderer_name) {
	case "svg":
	    svg.dotlabel(xp, yp, label, anchor, dotSize, fontSize, colour, fontFamily);
	    break;
	default:
	    throw new Exception("oops, invalid render name " ~ renderer_name);
	}
	return;
    }
    
    // ----------------------------------------------------------------------
    // The following methods should be independent of the particular renderer.
    // ----------------------------------------------------------------------
    
    void rule(string direction, double vmin, double vmax, double vtic, Vector3 anchorPoint,
	      double ticMarkSize, string numberFormat, double textOffset, int fontSize)
    // direction: one of "x", "y", "z"
    // vmin, vmax, vtic: range and increments for our rule
    // anchorPoint: rule will go through this point (if extended, maybe)
    // ticMarkSize: in world units
    // numberFormat: print specifier string, "%.1f" for example
    // textOffset: world units to anchor point for text labels
    // fontSize: in points
    {
	Vector3 p0 = Vector3(anchorPoint);
	Vector3 p1 = Vector3(anchorPoint);
	Vector3 dp = Vector3(0.0,0.0,0.0);
	Vector3 dpTic = Vector3(0.0,0.0,0.0);
	int n = to!int(1.000001*(vmax - vmin)/vtic); // be sure to get final tic mark 
	if (n > 50) {
	    writeln("You have asked for a lot of tic marks: ", n);
	}
	Vector3 dpText = Vector3(0.0,0.0,0.0);
	Vector3 textDirection;
	Vector3 textNormal;
	string textAnchor;
	switch(direction) {
	case "x":
	    p0.refx = vmin; p1.refx = vmax; dp.refx = vtic;
	    dpTic.refy = -ticMarkSize; // tic marks go in -y direction
	    dpText.refy = -textOffset;
	    textDirection = Vector3(1.0,0.0,0.0); // text horizontal
	    textNormal = Vector3(0.0,0.0,1.0);
	    textAnchor = "middle";
	    break;
	case "y":
	    p0.refy = vmin; p1.refy = vmax; dp.refy = vtic;
	    dpTic.refx = -ticMarkSize; // tic marks go in -x direction
	    dpText.refx = -textOffset;
	    textDirection = Vector3(1.0,0.0,0.0); // text horizontal
	    textNormal = Vector3(0.0,0.0,1.0);
	    textAnchor = "end";
	    break;
	case "z":
	    p0.refz = vmin; p1.refz = vmax; dp.refz = vtic;
	    dpTic.refx = -ticMarkSize; // tic marks go in -x direction
	    dpText.refx = -textOffset;
	    textDirection = Vector3(1.0,0.0,0.0);
	    textNormal = Vector3(0.0,1.0,0.0);
	    textAnchor = "end";
	    break;
	default:
	    throw new Exception("invalid direction for rule: "~direction);
	}
	begin_group(); // Keep tic marks together
	line(p0, p1);
	foreach(i; 0 .. n+1) {
	    Vector3 p = p0 + i*dp;
	    Vector3 p2 = p+dpTic;
	    line(p, p2);
	}
	end_group();
	// Label the tic marks
	foreach(i; 0 .. n+1) {
	    Vector3 p = p0 + i*dp;
	    Vector3 p3 = p+dpText;
	    double value = vmin + i*vtic;
	    string myText = format(numberFormat, value);
	    text(p3, myText, textDirection, textNormal, textAnchor, fontSize);
	}
    } // end rule()
    
    // Render functions for our geometric entities used in the flow simulation.
    void render(const Path pth, bool dashed=false, size_t n=30)
    {
	Vector3[] psample;
	double dt = 1.0/n;
	foreach(i; 0 .. n+1) { psample ~= pth(i*dt); }
	polyline(psample, dashed);
    }

    void render(const ParametricSurface surf, bool fill=true,
		bool stroke=true, bool dashed=false,
		size_t n=30, bool facets=false)
    {
	double dt = 1.0/n;
	// Sample boundaries, progressing counter-clockwise.
	Vector3[] p_boundary;
	foreach(i; 0 .. n+1) { p_boundary ~= surf(i*dt, 0.0); } // south boundary
	foreach(i; 1 .. n+1) { p_boundary ~= surf(1.0, i*dt); } // east
	foreach(i; 1 .. n+1) { p_boundary ~= surf(1.0-i*dt, 1.0); } // north
	foreach(i; 1 .. n+1) { p_boundary ~= surf(0.0, 1.0-i*dt); } // west
	if (fill) {
	    if (facets) {
		// Fill in the surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double r0 = i*dt; double r1 = r0+dt;
		    foreach(j; 0 .. n) {
			double s0 = j*dt; double s1 = s0+dt;
			Vector3[] psample;
			psample ~= surf(r0, s0); psample ~= surf(r1, s0);
			psample ~= surf(r1, s1); psample ~= surf(r0, s1);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
	    } else {
		// Fill a single polygon boundary.
		// This will allow a much smaller SVG file for 2D renderings.
		polygon(p_boundary, fill, false, false);
	    }
	} // end if fill
	if (stroke) {
	    polygon(p_boundary, false, stroke, dashed);
	}
    } // end render (ParametricSurface)

    void render(const ParametricVolume pvol, bool fill=true,
		bool stroke=true, bool dashed=false,
		size_t n=30, bool facets=false)
    {
	// Render the bounding surface and edges of the ParametricVolume.
	double dt = 1.0/n;
	// Sample boundaries, progressing counter-clockwise.
	Vector3[] p_bottom;
	foreach(i; 0 .. n+1) { p_bottom ~= pvol(i*dt, 0.0, 0.0); } // south edge
	foreach(i; 1 .. n+1) { p_bottom ~= pvol(1.0, i*dt, 0.0); } // east
	foreach(i; 1 .. n+1) { p_bottom ~= pvol(1.0-i*dt, 1.0, 0.0); } // north
	foreach(i; 1 .. n+1) { p_bottom ~= pvol(0.0, 1.0-i*dt, 0.0); } // west
	Vector3[] p_top;
	foreach(i; 0 .. n+1) { p_top ~= pvol(i*dt, 0.0, 1.0); } // south edge
	foreach(i; 1 .. n+1) { p_top ~= pvol(1.0, i*dt, 1.0); } // east
	foreach(i; 1 .. n+1) { p_top ~= pvol(1.0-i*dt, 1.0, 1.0); } // north
	foreach(i; 1 .. n+1) { p_top ~= pvol(0.0, 1.0-i*dt, 1.0); } // west
	Vector3[] p_west;
	foreach(i; 0 .. n+1) { p_west ~= pvol(0.0, 0.0, i*dt); }
	foreach(i; 1 .. n+1) { p_west ~= pvol(0.0, i*dt, 1.0); }
	foreach(i; 1 .. n+1) { p_west ~= pvol(0.0, 1.0, 1.0-i*dt); }
	foreach(i; 1 .. n+1) { p_west ~= pvol(0.0, 1.0-i*dt, 0.0); }
	Vector3[] p_east;
	foreach(i; 0 .. n+1) { p_east ~= pvol(1.0, 0.0, i*dt); }
	foreach(i; 1 .. n+1) { p_east ~= pvol(1.0, i*dt, 1.0); }
	foreach(i; 1 .. n+1) { p_east ~= pvol(1.0, 1.0, 1.0-i*dt); }
	foreach(i; 1 .. n+1) { p_east ~= pvol(1.0, 1.0-i*dt, 0.0); }
	Vector3[] p_south;
	foreach(i; 0 .. n+1) { p_south ~= pvol(i*dt, 0.0, 0.0); }
	foreach(i; 1 .. n+1) { p_south ~= pvol(1.0, 0.0, i*dt); }
	foreach(i; 1 .. n+1) { p_south ~= pvol(1.0-i*dt, 0.0, 1.0); }
	foreach(i; 1 .. n+1) { p_south ~= pvol(0.0, 0.0, 1.0-i*dt); }
	Vector3[] p_north;
	foreach(i; 0 .. n+1) { p_north ~= pvol(i*dt, 1.0, 0.0); }
	foreach(i; 1 .. n+1) { p_north ~= pvol(1.0, 1.0, i*dt); }
	foreach(i; 1 .. n+1) { p_north ~= pvol(1.0-i*dt, 1.0, 1.0); }
	foreach(i; 1 .. n+1) { p_north ~= pvol(0.0, 1.0, 1.0-i*dt); }
	if (fill) {
	    if (facets) {
		begin_group(); // hold all surfaces together
		// Fill in the bottom surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double r0 = i*dt; double r1 = r0+dt;
		    foreach(j; 0 .. n) {
			double s0 = j*dt; double s1 = s0+dt;
			Vector3[] psample;
			psample ~= pvol(r0, s0, 0.0); psample ~= pvol(r1, s0, 0.0);
			psample ~= pvol(r1, s1, 0.0); psample ~= pvol(r0, s1, 0.0);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
		// Fill in the top surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double r0 = i*dt; double r1 = r0+dt;
		    foreach(j; 0 .. n) {
			double s0 = j*dt; double s1 = s0+dt;
			Vector3[] psample;
			psample ~= pvol(r0, s0, 1.0); psample ~= pvol(r1, s0, 1.0);
			psample ~= pvol(r1, s1, 1.0); psample ~= pvol(r0, s1, 1.0);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
		// Fill in the west surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double t0 = i*dt; double t1 = t0+dt;
		    foreach(j; 0 .. n) {
			double s0 = j*dt; double s1 = s0+dt;
			Vector3[] psample;
			psample ~= pvol(0.0, s0, t0); psample ~= pvol(0.0, s0, t1);
			psample ~= pvol(0.0, s1, t1); psample ~= pvol(0.0, s1, t0);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
		// Fill in the east surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double t0 = i*dt; double t1 = t0+dt;
		    foreach(j; 0 .. n) {
			double s0 = j*dt; double s1 = s0+dt;
			Vector3[] psample;
			psample ~= pvol(1.0, s0, t0); psample ~= pvol(1.0, s0, t1);
			psample ~= pvol(1.0, s1, t1); psample ~= pvol(1.0, s1, t0);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
		// Fill in the south surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double r0 = i*dt; double r1 = r0+dt;
		    foreach(j; 0 .. n) {
			double t0 = j*dt; double t1 = t0+dt;
			Vector3[] psample;
			psample ~= pvol(r0, 0.0, t0); psample ~= pvol(r1, 0.0, t0);
			psample ~= pvol(r1, 0.0, t1); psample ~= pvol(r0, 0.0, t1);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
		// Fill in the north surface as (many) quadrilateral patches.
		begin_group();
		foreach(i; 0 .. n) {
		    double r0 = i*dt; double r1 = r0+dt;
		    foreach(j; 0 .. n) {
			double t0 = j*dt; double t1 = t0+dt;
			Vector3[] psample;
			psample ~= pvol(r0, 1.0, t0); psample ~= pvol(r1, 1.0, t0);
			psample ~= pvol(r1, 1.0, t1); psample ~= pvol(r0, 1.0, t1);
			polygon(psample, fill, false, false);
		    }
		}
		end_group();
		end_group(); // for all surfaces, together
	    } else {
		// Fill a single polygon boundary for each bounding surface.
		// This will allow a much smaller SVG file for 2D renderings.
		begin_group();
		polygon(p_bottom, fill, false, false);
		polygon(p_top, fill, false, false);
		polygon(p_west, fill, false, false);
		polygon(p_east, fill, false, false);
		polygon(p_south, fill, false, false);
		polygon(p_north, fill, false, false);
		end_group();
	    }
	} // end if fill
	if (stroke) {
	    begin_group();
	    polygon(p_bottom, false, stroke, dashed);
	    polygon(p_top, false, stroke, dashed);
	    polygon(p_west, false, stroke, dashed);
	    polygon(p_east, false, stroke, dashed);
	    polygon(p_south, false, stroke, dashed);
	    polygon(p_north, false, stroke, dashed);
	    end_group();
	}
    } // end render (ParametricVolume)
    
 } // end class Sketch
