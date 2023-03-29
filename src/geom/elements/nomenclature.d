// nomenclature.d
// Names related to geometry that will be used throughout the CFD code.

module geom.elements.nomenclature;

import std.conv: text;
import std.math: PI;

// Symbolic names and indices for the cells' faces.
// The names of the faces of the structured-grid blocks will be the same.
version (OLD_FACE_ORDER) {
    enum Face {
        north = 0,
        east = 1,
        south = 2,
        west = 3,
        top = 4,
        bottom = 5
    }
    string[] face_name = [ "north", "east", "south", "west", "top", "bottom" ];
 } else {
    // New face order, 2023-03-29, corresponding to the index notation used in Chicken.
    // iminus, iplus, jminus, jplus, kminus, kplus
    enum Face {
        west = 0, iminus = 0,
        east = 1, iplus = 1,
        south = 2, jminus = 2,
        north = 3, jplus = 3,
        bottom = 4, kminus = 4,
        top = 5, kplus = 5
    }
    string[] face_name = [ "west", "east", "south", "north", "bottom", "top" ];
 }

uint face_index(string name)
{
    switch ( name ) {
    case "west":
    case "iminus": return Face.west;
    case "east":
    case "iplus": return Face.east;
    case "south":
    case "jminus": return Face.south;
    case "north":
    case "jplus": return Face.north;
    case "bottom":
    case "kminus": return Face.bottom;
    case "top":
    case "kplus": return Face.top;
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
