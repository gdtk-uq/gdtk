/** parametricsurface.d
 * Geometry-building elements for our 3D world -- two-parameter surfaces.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 *          2017-11-29 refactored into a number of files.
 */

module geom.surface.parametricsurface;

import std.math;
import std.stdio;
import std.conv;
import std.algorithm;

import ntypes.complex;
import nm.number;
import geom;


// Nomenclature for the parametric distances, bounding paths and corners.
//
//         north
//     p01-------p11 s=1
//      |         |
// west |         | east
//      |         |
//     p00-------p10 s=0
//         south
//     r=0       r=1
//
// We'll try to use this notation consistently in the classes below.

class ParametricSurface {
public:
    abstract Vector3 opCall(double r, double s) const;
    abstract ParametricSurface dup() const;
    abstract override string toString() const;

    Vector3 area(int nr=10, int ns=10)
    {
        Vector3 vector_area; vector_area.set(0.0, 0.0, 0.0);
        double dr = 1.0/nr;
        double ds = 1.0/ns;
        Vector3[][] p;
        p.length = nr+1;
        foreach (i; 0 .. nr+1) {
            double r = dr * i;
            foreach (j; 0 .. ns+1) {
                double s = ds * j;
                p[i] ~= opCall(r, s);
            }
        }
        number dA;
        Vector3 centroid;
        Vector3 n, t1, t2;
        foreach (i; 0 .. nr) {
            foreach (j; 0 .. ns) {
                quad_properties(p[i][j], p[i+1][j], p[i+1][j+1], p[i][j+1],
                                centroid, n, t1, t2, dA);
                n *= dA;
                vector_area += n;
            }
        }
        return vector_area;
    } // end area()

} // end class ParametricSurface

void writeSurfaceAsVtkXml(ParametricSurface surf, string fileName, int nrPts, int nsPts)
{
    double dr = 1.0/(nrPts-1);
    double ds = 1.0/(nsPts-1);
    auto f = File(fileName, "w");
    f.writeln("<VTKFile type=\"StructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
    f.writefln("  <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">", 0, nrPts-1, 0, nsPts-1);
    f.writefln("    <Piece Extent=\"%d %d %d %d 0 0\">",  0, nrPts-1, 0, nsPts-1);
    f.writeln("      <Points>");
    f.writeln("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
    foreach (j; 0 .. nsPts) {
        foreach (i; 0 .. nrPts) {
            auto r = i*dr;
            auto s = j*ds;
            auto p = surf(r, s);
            f.writefln("       %20.16e %20.16e %20.16e", p.x, p.y, p.z);
        }
    }
    f.writeln("        </DataArray>");
    f.writeln("      </Points>");
    f.writeln("    </Piece>");
    f.writeln("  </StructuredGrid>");
    f.writeln("</VTKFile>");
    f.close();
}
