/**
 * svg.d  Render a drawing in Scalable Vector Graphics (SVG) format.
 *
 * This module provides a few convenient functions for rendering 
 * a 2D drawing in Scalable Vector Graphics format.
 * The main transformation is 
 * from a user-space coordinate system with (0,0) at the lower-left corner 
 * to the SVG coordinate system with (0,0) in the upper-left corner of the page.
 * Along the way, user-space units are converted to points because
 * Inkscape seems to behave better if everything is specified in points.
 *
 * Author(s):
 *     Peter J.
 * Version: 
 *     2005-10-25 first cut of Python version
 *     2016-05-18 Ported from cfcfd3/lib/cfpylib/svg_render.py.
 *
 */
module misc.svg;

import std.stdio;
import std.math;
import std.string;
import std.algorithm;

class SVGContext {
private:
    // Somewhere to keep the rendering configuration.
    double width;
    double height;
    string unitLength;
    double unitToPixels = 1.0;
    string lineColour = "black";
    string fillColour = "none";
    double lineWidth = 1.0; // points
    string lineCap = "round";
    string lineJoin = "round";
    double dashLength = 2.0*72.0/25.4; // 2mm in points
    double gapLength = 2.0*72.0/25.4;
    string title;
    string description;
    File f;

public:
    this(double width=120.0, double height=120.0, string unitLength="mm",
         string title="Untitled", string description="No description")
    // Creates a SVG graphics context with a particular canvas size.
    {
        this.width = width;
        this.height = height;
        this.unitLength = unitLength;
        // Anonymous units in SVG 1.1 are pixels with 90px=1in or 1.25px=1pt
        switch (unitLength) {
        case "in": unitToPixels = 90.0; break;
        case "mm": unitToPixels = 90.0/25.4; break;
        case "cm": unitToPixels = 90.0/2.54; break;
        case "m": unitToPixels = 90.0/0.0254; break;
        default:
            throw new Exception("SvgEnvironment: Unknown units " ~ unitLength);
        } // end switch
        setLineWidth(0.25);
        setDashArray();
        this.title = title;
        this.description = description;
        return;
    }

    double toPixelsX(double x)
    // Transforms x-coordinate from user-space to SVG space.
    // x: x-coordinate in units (with the origin in the lower-left corner)       
    // Returns points for SVG
    {
        return x * unitToPixels;
    }

    double toPixelsY(double y)
    // Transforms y-coordinate from user-space to SVG space.
    // y: y-coordinate in units (with the origin in the lower-left corner)
    // Returns points in SVG coordinate system (with the origin in the upper left)
    {
        return (height - y) * unitToPixels;
    }

    void setLineWidth(double w)
    // Sets line width (in mm).
    {
        lineWidth = w * 72.0 / 25.4;
        return;
    }

    void setLineColour(string colour)
    // 6 hexadecimal digits in form "#rrggbb"
    // aqua, black, blue, fichsia, gray, green, lime, maroon,
    // navy, olive, purple, red, silver, teal, white, yellow
    {
        lineColour = colour;
        return;
    }

    void setFillColour(string colour)
    {
        fillColour = colour;
        return;
    }

    void clearFillColour()
    {
        fillColour = "none";
        return;
    }
    
    void setDashArray(double dashLength=2.0, double gapLength=2.0)
    // Sets length of dashes and gaps (in mm).
    {
        this.dashLength = dashLength * 72.0/25.4;
        this.gapLength = gapLength * 72.0/25.4;
        return;
    }
    
    string styleStr(bool fill, bool stroke, bool dashed,
                    double fillOpacity=1.0, double strokeOpacity=1.0)
    // Assembles a suitable style specification string.
    // fill: flag to indicate that the interior should be painted
    // stroke: flag to indicate that the outline should be drawn
    // dashed: flag to indicate that the outline is dashed
    // Returns style string.
    {
        string style = "";
        if (fill) {
            style ~= format("fill:%s;fill-opacity:%g", fillColour, fillOpacity);
        } else {
            style ~= format("fill:%s", "none");
        }
        if (stroke) {
            style ~= format(";stroke:%s;stroke-width:%.2f;stroke-linecap:%s",
                            lineColour, lineWidth, lineCap);
            style ~= format(";stroke-opacity:%g", strokeOpacity);
        }
        if (dashed) {
            style ~= format(";stroke-dasharray: %.2f %.2f", dashLength, gapLength);
        }
        return style;
    }
    
    void open(string fileName="drawing.svg")
    // Opens the SVG file and writes the preamble.
    {
        f = File(fileName, "w");
        f.write("<?xml version=\"1.0\"?>\n");
        f.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"\n");
        f.write("\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n");
        f.write("<svg xmlns:svg=\"http://www.w3.org/2000/svg\"\n");
        f.write("     xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
        f.writef("     width=\"%.2f\" height=\"%.2f\"\n",
                 width*unitToPixels, height*unitToPixels);
        f.write("     version=\"1.0\">\n");
        f.writef("<title>%s</title>\n", title);
        f.writef("<desc>%s</desc>\n", description);
        f.write("<!-- drawing begins here -->\n");
        return;
    }

    void add_comment(string text)
    // Inserts a comment into the SVG file.
    {
        f.write("<!-- ");
        f.write(text);
        f.write(" -->\n");
        return;
    }
    
    void close()
    // Finishes off the SVG file and closes it.
    {
        f.write("<!-- end of drawing -->\n");
        f.write("</svg>\n");
        f.close();
        return;
    }

    void begin_group(string id="", double opacity=1.0)
    // id needs to be a unique identifier string, if supplied
    // opacity is in range 0.0, 1.0
    {
        f.write("<g ");
        if (id.length > 0) { f.writef(" id=\"%s\"", id); }
        f.writef("opacity=\"%g\" >\n", opacity);
        return;
    }

    void end_group()
    {
        f.write("</g>\n");
        return;
    }
    
    void line(double x1, double y1, double x2, double y2,
              bool dashed=false, double opacity=1.0)
    // Render a line from point 1 to point 2.
    {
        auto x1p = toPixelsX(x1); auto y1p = toPixelsY(y1);
        auto x2p = toPixelsX(x2); auto y2p = toPixelsY(y2);
        f.writef("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" ",
                 x1p, y1p, x2p, y2p);
        f.writef("style=\"%s\" />\n", styleStr(false, true, dashed, 1.0, opacity));
        return;
    }

    void polyline(double[] xlist, double[] ylist,
                  bool dashed=false, double opacity=1.0)
    // Render a polyline from lists of x and y coordinates.
    {
        auto x0 = toPixelsX(xlist[0]); auto y0 = toPixelsY(ylist[0]);
        f.writef("<path d=\"M %.2f,%.2f", x0, y0);
        foreach (i; 1 .. min(xlist.length, ylist.length)) {
            auto x = toPixelsX(xlist[i]); auto y = toPixelsY(ylist[i]);
            f.writef(" L %.2f,%.2f", x, y);
        }
        f.write("\""); // finish the coordinate string
        f.writef(" style=\"%s\"", styleStr(false, true, dashed, 1.0, opacity));
        f.write("/>\n");
        return;
    }

    void polygon(double[] xlist, double[] ylist, bool fill=true,
                 bool stroke=true, bool dashed=false,
                 string id="", double opacity=0.5)
    // Render a (closed) polygon from lists of x and y coordinates.
    // The default opacity is 0.5 so that covered-over objects are
    // still visible.
    {
        auto x0 = toPixelsX(xlist[0]); auto y0 = toPixelsY(ylist[0]);
        f.writef("<polygon points=\"%.2f,%.2f", x0, y0);
        foreach (i; 1 .. min(xlist.length, ylist.length)) {
            auto x = toPixelsX(xlist[i]); auto y = toPixelsY(ylist[i]);
            f.writef(" %.2f,%.2f", x, y);
        }
        f.write("\""); // finish the coordinate string
        f.writef(" style=\"%s\"", styleStr(fill, stroke, dashed, opacity, opacity));
        f.write("/>\n");
        return;
    }

    void arc(double x0, double y0, double x1, double y1, double xc, double yc,
             bool dashed=false, double opacity=1.0)
    // Render a circular arc from point 0 to point 1 with centre at c.
    // in user-space coordinates, comput the properties of the arc.
    {
        double dx0 = x0 - xc; double dy0 = y0 - yc;
        double dx1 = x1 - xc; double dy1 = y1 - yc;
        double r0 = sqrt(dx0*dx0 + dy0*dy0);
        double r1 = sqrt(dx1*dx1 + dy1*dy1);
        if (fabs(r0 - r1) > 1.0e-6*(r0+1.0)) {
            throw new Exception("Radii don't match");
        }
        double crossz = (dx0 * dy1 - dx1 * dy0); // z-component of vector product
        bool clockwise = (crossz < 0.0);
        // now, do the rendering in SVG space (in points)
        double x0p = toPixelsX(x0); double y0p = toPixelsY(y0);
        double x1p = toPixelsX(x1); double y1p = toPixelsY(y1);
        double rp = r0 * unitToPixels;
        double x_axis_rotation = 0.0;
        int large_arc_flag = 0;
        int sweep_flag = (clockwise)? 1 : 0;
        f.writef("<path d=\"M %.2f %.2f A %.2f %.2f, %.2f, %d, %d, %.2f %.2f\" ",
                 x0p, y0p, rp, rp, x_axis_rotation, large_arc_flag,
                 sweep_flag, x1p, y1p);
        f.writef("style=\"%s\"", styleStr(false, true, dashed, 1.0, opacity));
        f.write("/>\n");
        return;
    }
    
    void circle(double x, double y, double r, bool fill=true,
                bool stroke=true, bool dashed=false,
                string id="", double opacity=0.5)
    // Render a circle of radius r at centre (x,y).
    // The default opacity is 0.5 so that covered-over objects are
    // still visible.
    {
        double xp = toPixelsX(x); double yp = toPixelsY(y);
        double rp = r * unitToPixels;
        f.writef("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" ", xp, yp, rp);
        f.writef("style=\"%s\" />\n", styleStr(fill, stroke, dashed, opacity, opacity));
        return;
    }

    void bezier3(double x0, double y0, double x1, double y1,
                 double x2, double y2, double x3, double y3, 
                 bool dashed=false, double opacity=1.0)
    // Render a thrid-order Bezier curve.
    {
        double x0p = toPixelsX(x0); double y0p = toPixelsY(y0);
        double x1p = toPixelsX(x1); double y1p = toPixelsY(y1);
        double x2p = toPixelsX(x2); double y2p = toPixelsY(y2);
        double x3p = toPixelsX(x3); double y3p = toPixelsY(y3);
        f.writef("<path d=\"M %.2f %.2f ", x0p, y0p);
        f.writef("C %.2f %.2f %.2f %.2f %.2f %.2f\" ",
                 x1p, y1p, x2p, y2p, x3p, y3p);
        f.writef("style=\"%s\" />\n", styleStr(false, true, dashed, 1.0, opacity));
        return;
    }
    
    void text(double x, double y, string textString,
              double angle=0.0, string anchor="start",
              int fontSize=10, string colour="black", string fontFamily="sanserif")
    // Render the textString at point (x,y).
    // x: x-coordinate of anchor in user-space
    // y: y-coordinate of anchor in user-space
    // textString: string of characters to render
    // angle: angle (in degrees) of text line wrt x-axis (counterclockwise is positive)
    // anchor: one of "start", "middle" or "end"
    // fontSize: size of font in points
    // colour: of rendered text
    // fontFamily: a font name known in SVG
    {
        double xp = toPixelsX(x); double yp = toPixelsY(y);
        f.writef("<text x=\"%.2f\" y=\"%.2f\"", xp, yp);
        if (angle != 0.0) {
            // my angle is positive counterclockwise.
            f.writef(" transform=\"rotate(%.2f,%.2f,%.2f)\" ", -angle, xp, yp);
        }
        f.writef(" style=\"font-size:%dpt;text-anchor:%s;fill:%s;stroke:none",
                fontSize, anchor, colour);
        f.writef(";font-weight:normal;font-family:%s\" ", fontFamily);
        f.write(">"); // end of opening tag
        f.write(textString);
        f.write("</text>\n");
        return;
    }

    void dotlabel(double x, double y, string label="",
                  string anchor="middle", double dotSize=2.0,
                  int fontSize=10, string colour="black", string fontFamily="sanserif")
    // Render a dot with a text label.
    // x: x-coordinate in user-space
    // y: y-coordinate in user-space
    // label: label text
    // anchor: anchor location on label
    // dotSize: dot diameter in mm
    // fontSize: font size in points
    // colour: of both the label and the dot
    // fontFamily: a font name known in SVG
    {
        double xp = toPixelsX(x); double yp = toPixelsY(y);
        double rp = dotSize/2.0 * 72.0/25.4;
        f.writef("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" ", xp, yp, rp);
        f.writef("style=\"fill:%s\" />\n", colour);
        if (label.length > 0) {
            // put the label slightly above the dot.
            f.writef("<text x=\"%.2f\" y=\"%.2f\"", xp, yp-1.5*rp);
            f.writef(" style=\"font-size:%dpt;text-anchor:%s;fill:%s;stroke:none",
                     fontSize, anchor, colour);
            f.writef(";font-weight:normal;font-family:%s\" ", fontFamily);
            f.write(">"); // end of opening tag
            f.write(label);
            f.write("</text>\n");
        }
        return;
    }
} // end class SVGContext   

version(svg_test) {
    import util.msg_service;
    int main() {
        // [TODO] fill me in with something meaningful.
        return 0;
    }
} // end svg_test
