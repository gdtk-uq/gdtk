/**
 * sketch.d  Sketch geometric elements using a selected renderer.
 *
 * These are the user-callable functions for sketching geometric elements
 * while preparing the data for a flow simulation.  They are modelled on
 * the capabilities of the SVG renderer but work in the geometric space
 * of the flow simulation, where the units of lengths are metres.
 * The "canvas" referred to below is a virtual canvas, where the units
 * of lengths are millimetres.
 *
 * Author(s):
 *     Peter J.
 * Version:
 * 2016-08-10 Let's start over with just the 2D SVG renderer
 *     but keep in mind that we'll want to be able to use Tcl/Tk,
 *     Cairo and OpenGL, eventually.
 * 2018-03-11 Add an X render, provided by GNU libplot.
 */

module geom.misc.sketch;

import std.conv;
import std.stdio;
import core.stdc.stdio: stdin, stdout, stderr;
import std.math;
import std.string;
import std.algorithm;

import ntypes.complex;
import nm.number;

import misc.svg;
version(with_libplot) {
    import libplot;
}
import geom;

const double dpi = 90.0;  // Expected screen resolution.

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

enum Renderer {svg, xplot}
enum Projection {xyortho, isometric, oblique}
string[] rendererNames = ["svg", "xplot"];
string[] projectionNames = ["xyortho", "isometric", "oblique"];

class Sketch {
public:
    // We have some state and configuration data.
    Renderer myRenderer = Renderer.svg;
    SVGContext svg;
    version(with_libplot) {
        plPlotter* myXplotter;
        double xplot_line_width = 0.25; // Line thickness in mm.
        string xplot_bgcolourname = "white";
        string xplot_pencolourname = "black";
        string xplot_fillcolourname = "yellow";
        int xplot_greylevel = 0x8000; // 50% of 0xffff
    }
    Projection myProjection = Projection.xyortho;
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
    Extents canvas; // Virtual-canvas or device-space in mm.
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
    this(string renderer_name="svg", string projection_name="xyortho",
         double[] canvas_mm = [0.0, 0.0, 120.0, 120.0])
    {
        if (!find(rendererNames, renderer_name).empty) {
            switch (renderer_name) {
            case "svg": myRenderer = Renderer.svg; break;
            case "xplot": myRenderer = Renderer.xplot; break;
            default: myRenderer = Renderer.svg;
            }
        } else {
            throw new Exception("Unknown renderer: " ~ renderer_name);
        }
        if (!find(projectionNames, projection_name).empty) {
            switch (projection_name) {
            case "xyortho": myProjection = Projection.xyortho; break;
            case "isometric": myProjection = Projection.isometric; break;
            case "oblique": myProjection = Projection.oblique; break;
            default: myProjection = Projection.xyortho;
            }
        } else {
            throw new Exception("Unknown projection: " ~ projection_name);
        }
        set_to_identity(proj_mat);
        set_to_identity(view_mat); // consistent with default view
        final switch (myProjection) {
        case Projection.xyortho:
            proj_mat[2][2] = 0.0; // x,y unchanged, z eliminated
            break;
        case Projection.isometric:
            double cos30 = cos(30.0/180.0*PI);
            double sin30 = sin(30.0/180.0*PI);
            proj_mat[0][0] =  cos30; proj_mat[0][1] = 0.0; proj_mat[0][2] = -cos30; proj_mat[0][3] = 0.0;
            proj_mat[1][0] = -sin30; proj_mat[1][1] = 1.0; proj_mat[1][2] = -sin30; proj_mat[1][3] = 0.0;
            proj_mat[2][0] =    0.0; proj_mat[2][1] = 0.0; proj_mat[2][2] =    0.0; proj_mat[2][3] = 0.0;
            proj_mat[3][0] =    0.0; proj_mat[3][1] = 0.0; proj_mat[3][2] =    0.0; proj_mat[3][3] = 1.0;
            break;
        case Projection.oblique:
            double zfactor = 1.0/(2.0*sqrt(2.0));
            proj_mat[0][0] =  1.0; proj_mat[0][1] = 0.0; proj_mat[0][2] = -zfactor; proj_mat[0][3] = 0.0;
            proj_mat[1][0] =  0.0; proj_mat[1][1] = 1.0; proj_mat[1][2] = -zfactor; proj_mat[1][3] = 0.0;
            proj_mat[2][0] =  0.0; proj_mat[2][1] = 0.0; proj_mat[2][2] =      0.0; proj_mat[2][3] = 0.0;
            proj_mat[3][0] =  0.0; proj_mat[3][1] = 0.0; proj_mat[3][2] =      0.0; proj_mat[3][3] = 1.0;
        }
        // Set corner positions (in mm) for the canvas.
        canvas.set(canvas_mm[0], canvas_mm[1], canvas_mm[2], canvas_mm[3]);
        // For the model coordinates, assume unit space because we don't yet know better.
        viewport.set(0.0, 0.0, 1.0, 1.0);
        final switch(myRenderer) {
        case Renderer.svg:
            // We're good for now, since we'll just be writing to a file
            // for each SVG plot.
            break;
        case Renderer.xplot:
            version(with_libplot) {
                // We need to get the Xplotter started up.
                plPlotterParams* params = pl_newplparams();
                string bitmapsize = format("%dx%d", to!int(canvas.width*dpi/25.4),
                                           to!int(canvas.height*dpi/25.4));
                pl_setplparam(params, "BITMAPSIZE", cast(void*)toStringz(bitmapsize));
                pl_setplparam(params, "VANISH_ON_DELETE", cast(void*)toStringz("yes"));
                pl_setplparam(params, "USE_DOUBLE_BUFFERING", cast(void*)toStringz("yes"));
                // Create an X Plotter with the specified parameters.
                if ((myXplotter = pl_newpl_r("X", stdin, stdout, stderr, params)) == null) {
                    throw new Error("Couldn't create X Plotter.");
                }
            } else {
                throw new Exception("You asked for xplot renderer but libplot is not included.");
            }
        }
    } // end this

    ~this()
    {
        // Cleanup
        final switch(myRenderer) {
        case Renderer.svg:
            // Do nothing
            break;
        case Renderer.xplot:
            version(with_libplot) {
                if (pl_deletepl_r(myXplotter) < 0) {
                    throw new Error("Couldn't delete X Plotter\n");
                }
            }
        }
    }

    override string toString()
    {
        string str = "Sketch(";
        str ~= "renderer=\""~rendererNames[myRenderer]~"\"";
        str ~= ", projection=\""~projectionNames[myProjection]~"\"";
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
        ph[0] = p.x.re; ph[1] = p.y.re; ph[2] = p.z.re; ph[3] = w;
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
        p.set(ph[0], ph[1], ph[2]); w = ph[3];
    } // end apply_transform

    void look_at(const Vector3 eye, const Vector3 centre, const Vector3 up)
    // Set the view matrix so that we have a new view of the 3D model.
    // In the new coordinates, we are looking back along the znew-axis,
    // with ynew up and xnew toward the right.
    {
        this.eye = eye;
        this.centre = centre;
        this.up = up;
        Vector3 znew = (eye - centre); znew.normalize();
        Vector3 xnew = (cross(up, znew)); xnew.normalize();
        Vector3 ynew = (cross(znew, xnew)); ynew.normalize();
        view_mat[0][0] = xnew.x.re; view_mat[0][1] = xnew.y.re; view_mat[0][2] = xnew.z.re;
        view_mat[1][0] = ynew.x.re; view_mat[1][1] = ynew.y.re; view_mat[1][2] = ynew.z.re;
        view_mat[2][0] = znew.x.re; view_mat[2][1] = znew.y.re; view_mat[2][2] = znew.z.re;
        view_mat[0][3] = -(dot(centre,xnew).re);
        view_mat[1][3] = -(dot(centre,ynew).re);
        view_mat[2][3] = -(dot(centre,znew).re);
        view_mat[3][0] = 0.0; view_mat[3][1] = 0.0; view_mat[3][2] = 0.0; view_mat[3][3] = 1.0;
    } // end look_at()

    void start(string file_name="sketch.svg")
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg = new SVGContext(canvas.width, canvas.height, "mm", title, description);
            svg.open(file_name);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                if (pl_openpl_r(myXplotter) < 0) { throw new Error("Couldn't open X Plotter."); }
                // Specify the virtual-canvas coordinate system,
                // starting at (0,0) bottom-left.
                pl_fspace_r(myXplotter, 0.0, 0.0, canvas.width, canvas.height);
                xplot_line_width = 0.5; // Line thickness in mm.
                xplot_bgcolourname = "white";
                xplot_pencolourname = "black";
                xplot_fillcolourname = "yellow";
                pl_flinewidth_r(myXplotter, xplot_line_width);
                pl_joinmod_r(myXplotter, toStringz("round"));
                pl_capmod_r(myXplotter, toStringz("round"));
                pl_pentype_r(myXplotter, 1);
                pl_linemod_r(myXplotter, toStringz("solid"));
                pl_filltype_r(myXplotter, xplot_greylevel);
                pl_bgcolorname_r(myXplotter, toStringz(xplot_bgcolourname));
                pl_pencolorname_r(myXplotter, toStringz(xplot_pencolourname));
                pl_fillcolorname_r(myXplotter, toStringz(xplot_fillcolourname));
                pl_erase_r(myXplotter);
            }
        }
    } // end start()

    void finish()
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.close();
            break;
        case Renderer.xplot:
            version(with_libplot) {
                if (pl_closepl_r(myXplotter) < 0) { throw new Error("Couldn't close X Plotter."); }
                writeln("Rendering finished; press ENTER to continue.");
                string junk_text = readln();
            }
        }
    } // end finish()

    void setLineWidth(double width)
    // Sets line width in mm on our virtual canvas.
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.setLineWidth(width); // SVG canvas was set up in mm.
            break;
        case Renderer.xplot:
            version(with_libplot) {
                xplot_line_width = width;
                pl_flinewidth_r(myXplotter, width);
            }
        }
        return;
    } // end setLineWidth()

    void setLineColour(string colour)
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.setLineColour(colour);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                xplot_pencolourname = colour;
                pl_pencolorname_r(myXplotter, toStringz(colour));
            }
        }
        return;
    } // end setLineColour()

    void setFillColour(string colour)
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.setFillColour(colour);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                xplot_fillcolourname = colour;
                pl_filltype_r(myXplotter, xplot_greylevel);
                pl_fillcolorname_r(myXplotter, toStringz(colour));
            }
        }
        return;
    } // end setFillColour()

    void clearFillColour()
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.clearFillColour();
            break;
        case Renderer.xplot:
            version(with_libplot) {
                xplot_fillcolourname = "none";
                pl_fillcolorname_r(myXplotter, "none");
            }
        }
        return;
    } // end clearFillColour()

    void setDashArray(double dashLength=2.0, double gapLength=2.0)
    // Sets length of dashes and gaps, in units appropriate to the renderer.
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.setDashArray(dashLength, gapLength); // in mm
            break;
        case Renderer.xplot:
            version(with_libplot) {
            // [TODO]
            }
        }
        return;
    } // end setDashArray()

    double toCanvasX(double x)
    // Map from user-space to canvas-space.
    {
        return canvas.x0 + (x - viewport.x0)/viewport.width*canvas.width;
    }

    double toCanvasY(double y)
    // Map from user-space to canvas-space.
    {
        return canvas.y0 + (y - viewport.y0)/viewport.height*canvas.height;
    }

    void begin_group(string id="", double opacity=1.0)
    // id needs to be a unique identifier string, if supplied
    // opacity is in range 0.0, 1.0
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.begin_group(id, opacity);
            break;
        case Renderer.xplot:
            // not applicable
        }
        return;
    }

    void end_group()
    {
        final switch(myRenderer) {
        case Renderer.svg:
            svg.end_group();
            break;
        case Renderer.xplot:
            // not applicable
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
        auto x0 = toCanvasX(p0tmp.x.re); auto y0 = toCanvasY(p0tmp.y.re);
        auto x1 = toCanvasX(p1tmp.x.re); auto y1 = toCanvasY(p1tmp.y.re);
        final switch(myRenderer) {
        case Renderer.svg:
            svg.line(x0, y0, x1, y1, dashed);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                pl_filltype_r(myXplotter, 0);
                pl_pentype_r(myXplotter, 1);
                pl_flinewidth_r(myXplotter, xplot_line_width);
                if (dashed) {
                    pl_linemod_r(myXplotter, toStringz("longdashed"));
                } else {
                    pl_linemod_r(myXplotter, toStringz("solid"));
                }
                pl_fline_r(myXplotter, x0, y0, x1, y1);
                pl_flushpl_r(myXplotter);
            }
        }
        return;
    } // end line()

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
            xlist ~= toCanvasX(p.x.re); ylist ~= toCanvasY(p.y.re);
        }
        final switch(myRenderer) {
        case Renderer.svg:
            svg.polyline(xlist, ylist, dashed);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                pl_filltype_r(myXplotter, 0);
                pl_pentype_r(myXplotter, 1);
                pl_flinewidth_r(myXplotter, xplot_line_width);
                if (dashed) {
                    pl_linemod_r(myXplotter, toStringz("longdashed"));
                } else {
                    pl_linemod_r(myXplotter, toStringz("solid"));
                }
                foreach (i; 1 .. xlist.length) {
                    pl_fline_r(myXplotter, xlist[i-1], ylist[i-1], xlist[i], ylist[i]);
                }
                pl_flushpl_r(myXplotter);
            }
        }
        return;
    } // end polyline()

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
            xlist ~= toCanvasX(p.x.re); ylist ~= toCanvasY(p.y.re);
        }
        final switch(myRenderer) {
        case Renderer.svg:
            svg.polygon(xlist, ylist, fill, stroke, dashed);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                if (fill) {
                    pl_filltype_r(myXplotter, xplot_greylevel);
                } else {
                    pl_filltype_r(myXplotter, 0);
                }
                if (stroke) {
                    pl_pentype_r(myXplotter, 1);
                    pl_flinewidth_r(myXplotter, xplot_line_width);
                } else {
                    pl_pentype_r(myXplotter, 0);
                }
                if (dashed) {
                    pl_linemod_r(myXplotter, toStringz("longdashed"));
                } else {
                    pl_linemod_r(myXplotter, toStringz("solid"));
                }
                pl_fmove_r(myXplotter, xlist[0], ylist[0]);
                foreach (i; 1 .. xlist.length) {
                    pl_fline_r(myXplotter, xlist[i-1], ylist[i-1], xlist[i], ylist[i]);
                }
                pl_fline_r(myXplotter, xlist[0], ylist[0], xlist[$-1], ylist[$-1]);
                pl_closepath_r(myXplotter);
                pl_endpath_r(myXplotter);
                pl_flushpl_r(myXplotter);
            }
        }
        return;
    } // end polyline()

    void arc(const Vector3 p0, const Vector3 p1, const Vector3 pc,
             bool dashed=false)
    {
        // This is a purely 2D xy-plane function.
        double x0 = toCanvasX(p0.x.re); double y0 = toCanvasY(p0.y.re);
        double x1 = toCanvasX(p1.x.re); double y1 = toCanvasY(p1.y.re);
        double xc = toCanvasX(pc.x.re); double yc = toCanvasY(pc.y.re);
        final switch(myRenderer) {
        case Renderer.svg:
            svg.arc(x0, y0, x1, y1, xc, yc, dashed);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                pl_filltype_r(myXplotter, 0);
                pl_pentype_r(myXplotter, 1);
                pl_flinewidth_r(myXplotter, xplot_line_width);
                if (dashed) {
                    pl_linemod_r(myXplotter, toStringz("longdashed"));
                } else {
                    pl_linemod_r(myXplotter, toStringz("solid"));
                }
                pl_farc_r(myXplotter, xc, yc, x0, y0, x1, y1);
                pl_flushpl_r(myXplotter);
            }
        }
        return;
    } // end arc()

    void circle(const Vector3 pc, double r, bool fill=true,
                bool stroke=true, bool dashed=false)
    {
        // This is a purely 2D xy-plane function.
        // r is in mm on the canvas.
        double xc = toCanvasX(pc.x.re); double yc = toCanvasY(pc.y.re);
        final switch(myRenderer) {
        case Renderer.svg:
            svg.circle(xc, yc, r, fill, stroke, dashed);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                if (fill) {
                    pl_filltype_r(myXplotter, xplot_greylevel);
                } else {
                    pl_filltype_r(myXplotter, 0);
                }
                if (stroke) {
                    pl_pentype_r(myXplotter, 1);
                    pl_flinewidth_r(myXplotter, xplot_line_width);
                } else {
                    pl_pentype_r(myXplotter, 0);
                }
                if (dashed) {
                    pl_linemod_r(myXplotter, toStringz("longdashed"));
                } else {
                    pl_linemod_r(myXplotter, toStringz("solid"));
                }
                pl_fcircle_r(myXplotter, xc, yc, r);
                pl_flushpl_r(myXplotter);
            }
        }
        return;
    } // end circle()

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
        double x0 = toCanvasX(p0tmp.x.re); double y0 = toCanvasY(p0tmp.y.re);
        double x1 = toCanvasX(p1tmp.x.re); double y1 = toCanvasY(p1tmp.y.re);
        double x2 = toCanvasX(p2tmp.x.re); double y2 = toCanvasY(p2tmp.y.re);
        double x3 = toCanvasX(p3tmp.x.re); double y3 = toCanvasY(p3tmp.y.re);
        final switch(myRenderer) {
        case Renderer.svg:
            svg.bezier3(x0, y0, x1, y1, x2, y2, x3, y3, dashed);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                pl_filltype_r(myXplotter, 0);
                pl_pentype_r(myXplotter, 1);
                pl_flinewidth_r(myXplotter, xplot_line_width);
                if (dashed) {
                    pl_linemod_r(myXplotter, toStringz("longdashed"));
                } else {
                    pl_linemod_r(myXplotter, toStringz("solid"));
                }
                pl_fbezier3_r(myXplotter, x0, y0, x1, y1, x2, y2, x3, y3);
                pl_flushpl_r(myXplotter);
            }
        }
        return;
    } // end bezier()

    void text(const Vector3 p, string textString,
              const double angle=0.0, // angle (in degrees) of text line wrt x-axis in canvas system
              // (counterclockwise is positive)
              string anchor="start",
              int fontSize=10, string colour="black", string fontFamily="sanserif")
    {
        Vector3 p_tmp = Vector3(p);
        double w_p = 1.0; double w_dir = 1.0; double w_n = 1.0;
        apply_transform("view", p_tmp, w_p); apply_transform("projection", p_tmp, w_p);
        double xp = toCanvasX(p_tmp.x.re); double yp = toCanvasY(p_tmp.y.re);
        final switch(myRenderer) {
        case Renderer.svg:
            svg.text(xp, yp, textString, angle, anchor, fontSize, colour, fontFamily);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                pl_pencolorname_r(myXplotter, toStringz(colour));
                // Convert font size from points to virtual-canvas millimetres.
                double fontHeight = fontSize/72.0*25.4;
                pl_ftextangle_r(myXplotter, angle);
                pl_fmove_r(myXplotter, xp, yp);
                pl_fontname_r(myXplotter, toStringz(fontFamily));
                pl_ffontsize_r(myXplotter, fontHeight);
                int horizJust, vertJust;
                switch (anchor) {
                case "start": horizJust = 'l'; vertJust = 'x'; break;
                case "middle": horizJust = 'c'; vertJust = 'x'; break;
                case "end": horizJust = 'r'; vertJust = 'x'; break;
                default:
                    horizJust = 'l'; vertJust = 'x';
                }
                pl_alabel_r(myXplotter, horizJust, vertJust, toStringz(textString));
                pl_flushpl_r(myXplotter);
                pl_pencolorname_r(myXplotter, toStringz(xplot_pencolourname));
            }
        }
        return;
    } // end text()

    void dotlabel(const Vector3 p, string label="",
                  string anchor="middle", double dotSize=2.0,
                  int fontSize=10, string colour="black", string fontFamily="sanserif")
    // dotSize is already in virtual-canvas mm
    // fontSize is in points
    {
        Vector3 p_tmp = Vector3(p); double w_p = 1.0;
        apply_transform("view", p_tmp, w_p); apply_transform("projection", p_tmp, w_p);
        double xp = toCanvasX(p_tmp.x.re); double yp = toCanvasY(p_tmp.y.re);
        final switch(myRenderer) {
        case Renderer.svg:
            svg.dotlabel(xp, yp, label, anchor, dotSize, fontSize, colour, fontFamily);
            break;
        case Renderer.xplot:
            version(with_libplot) {
                pl_pencolorname_r(myXplotter, toStringz(colour));
                pl_fillcolorname_r(myXplotter, toStringz(colour));
                pl_fcircle_r(myXplotter, xp, yp, dotSize/2);
                if (label.length > 0) {
                    // Convert font size from points to virtual-canvas millimetres.
                    double fontHeight = fontSize/72.0*25.4;
                    pl_ftextangle_r(myXplotter, 0.0);
                    pl_fmove_r(myXplotter, xp, yp+0.75*dotSize);
                    pl_fontname_r(myXplotter, toStringz(fontFamily));
                    pl_ffontsize_r(myXplotter, fontHeight);
                    int horizJust, vertJust;
                    switch (anchor) {
                    case "start": horizJust = 'l'; vertJust = 'x'; break;
                    case "middle": horizJust = 'c'; vertJust = 'x'; break;
                    case "end": horizJust = 'r'; vertJust = 'x'; break;
                    default:
                        horizJust = 'l'; vertJust = 'x';
                    }
                    pl_alabel_r(myXplotter, horizJust, vertJust, toStringz(label));
                }
                pl_flushpl_r(myXplotter);
                pl_pencolorname_r(myXplotter, toStringz(xplot_pencolourname));
                pl_fillcolorname_r(myXplotter, toStringz(xplot_fillcolourname));
            }
        }
        return;
    } // end dotlabel()

    // ----------------------------------------------------------------------
    // The following methods should be independent of the particular renderer.
    // ----------------------------------------------------------------------

    void rule(string direction, double vmin, double vmax, double vtic, Vector3 anchorPoint,
              double ticMarkSize, string numberFormat, double textOffset, double textAngle, int fontSize)
    // direction: one of "x", "y", "z"
    // vmin, vmax, vtic: range and increments for our rule
    // anchorPoint: rule will go through this point (if extended, maybe)
    // ticMarkSize: in world units
    // numberFormat: print specifier string, "%.1f" for example
    // textOffset: world units to anchor point for text labels
    // textAngle: angle (in degrees) of text line wrt x-axis in canvas system
    //            (counterclockwise is positive)
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
        string textAnchor;
        switch(direction) {
        case "x":
            p0.x = vmin; p1.x = vmax; dp.x = vtic;
            dpTic.y = -ticMarkSize; // tic marks go in -y direction
            dpText.y = -textOffset;
            textAnchor = "middle";
            break;
        case "y":
            p0.y = vmin; p1.y = vmax; dp.y = vtic;
            dpTic.x = -ticMarkSize; // tic marks go in -x direction
            dpText.x = -textOffset;
            textAnchor = "end";
            break;
        case "z":
            p0.z = vmin; p1.z = vmax; dp.z = vtic;
            dpTic.x = -ticMarkSize; // tic marks go in -x direction
            dpText.x = -textOffset;
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
            text(p3, myText, textAngle, textAnchor, fontSize);
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
