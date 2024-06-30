// nomenclature.d
// Names related to geometry that will be used throughout the CFD code.

module geom.elements.nomenclature;

import std.conv: text;
import std.math: PI;

// Symbolic names and indices for the cells' faces.
// The names of the faces of the structured-grid blocks will be the same.
// New face order, 2023-03-29, corresponding to the index notation used in Chicken.
enum Face {
    west = 0, iminus = 0,
    east = 1, iplus = 1,
    south = 2, jminus = 2,
    north = 3, jplus = 3,
    bottom = 4, kminus = 4,
    top = 5, kplus = 5
}
string[] face_name = [ "west", "east", "south", "north", "bottom", "top" ];
// Chicken names:       iminus, iplus,  jminus,  jplus,   kminus,   kplus

uint face_index(string name)
{
    switch ( name ) {
    case "west": case "iminus": return Face.west;
    case "east": case "iplus": return Face.east;
    case "south": case "jminus": return Face.south;
    case "north": case "jplus": return Face.north;
    case "bottom": case "kminus": return Face.bottom;
    case "top": case "kplus": return Face.top;
    default:
        throw new Error(text("Invalid face name: ", name));
    }
} // end face_index

// VTK cell types, for use when writing and reading VTK files.
// With wedge and pyramid items from SU2 Mesh File documentation.
enum VTKElement {
    vertex = 1,
    polyvertex = 2,
    line = 3,
    polyline = 4,
    triangle = 5,
    triangle_strip = 6,
    polygon = 7,
    pixel = 8,
    quad = 9,
    tetra = 10,
    voxel = 11,
    hexahedron = 12,
    wedge = 13,
    pyramid = 14
}

double radians(double degrees) { return PI*degrees/180.0; }

/**
 * Quaternion represented as four scalars.
 *
 * Mathematical notation --> D implementation:
 *
 * q = q0 + q1 \hat{i} + q2 \hat{j} + q2 \hat{k}
 *
 * Quaternion q = (q0, q1, q2, q3);
 * double rot_scalar = q[0];
 */
alias Quaternion = double[4];
