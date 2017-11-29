/** parametricvolume.d
 * Geometry-building elements for our 3D world -- three-parameter volumes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-04-07 first code
 */

module geom.volume.parametricvolume;

import geom.elements;

// Nomenclature for the parametric distances, bounding surfaces, paths and corners.
//
// t=1 at top surface
//         north
//    p011-------p111 s=1
//      |         |
// west |   Top   | east
//      |         |
//    p001-------p101 s=0
//         south
//     r=0       r=1
//
//
// t=0 at Bottom surface
//         north
//    p010-------p110 s=1
//      |         |
// west |  Bottom | east
//      |         |
//    p000-------p100 s=0
//         south
//     r=0       r=1
//
// Faces:
// North = face[0]; East = face[1]; South = face[2]; West = face[3];
// Top = face[4]; Bottom = face[5]
//
// Corners:
// Bottom surface: p000 == p[0]; p100 == p[1]; p110 == p[2]; p010 == p[3]
// Top surface   : p001 == p[4]; p101 == p[5]; p111 == p[6]; p011 == p[7]
//
// Edges:
// edge[0] p[0] --> p[1] around Bottom surface
//     [1] p[1] --> p[2]
//     [2] p[3] --> p[2]
//     [3] p[0] --> p[3]
//
//     [4] p[4] --> p[5] around Top surface
//     [5] p[5] --> p[6]
//     [6] p[7] --> p[6]
//     [7] p[4] --> p[7]
//
//     [8] p[0] --> p[4] connecting Bottom to Top
//     [9] p[1] --> p[5]
//    [10] p[2] --> p[6]
//    [11] p[3] --> p[7]
//
// We'll try to use this notation consistently in the classes of this package.

class ParametricVolume {
public:
    abstract Vector3 opCall(double r, double s, double t) const;
    abstract ParametricVolume dup() const;
    abstract override string toString() const;
} // end class ParametricVolume
