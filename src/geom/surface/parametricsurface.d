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
} // end class ParametricSurface
