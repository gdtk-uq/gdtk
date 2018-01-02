// nomenclature.d
// Names related to geometry that will be used throughout the CFD code.

module geom.elements.nomenclature;

import std.conv: text;
import std.math: PI;

// Symbolic names and indices for the cells' faces.
// The names of the faces of the structured-grid blocks will be the same.
enum Face {
    north = 0,
    east = 1,
    south = 2,
    west = 3,
    top = 4,
    bottom = 5
}

string[] face_name = [ "north", "east", "south", "west", "top", "bottom" ];
uint face_index(string name)
{
    switch ( name ) {
    case "north": return Face.north;
    case "east": return Face.east;
    case "south": return Face.south;
    case "west": return Face.west;
    case "top": return Face.top;
    case "bottom": return Face.bottom;
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
