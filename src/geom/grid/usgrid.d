/**
 * usgrid.d -- unstructured-grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-11-02 First code
 */

module geom.grid.usgrid;

import std.string;
import std.array;
import std.conv;
import std.stdio;
import std.file;
import std.algorithm;
import std.array;
import std.format;
import std.math;
import std.datetime;
import gzip;
import nm.complex;
import nm.number;

import geom;
import paver: PavedGrid, POINT_LIST, FACE_LIST, CELL_LIST;
import paver2d;

//-----------------------------------------------------------------
// For the USGCell types, we will have only the linear elements,
// with straight-line edges joining the vertices.
// We will use VTK names and vertex ordering to align with su2 file format,
// but we will have different integer tags and thus need a translation array
// between the tag values.  This is vtk_element_types below.
// 2016-dec-01: Kyle has found that Pointwise writes the su2 file with GMSH
// vertex ordering for the wedge elements.  We retain the VTK ordering.
// 2022-oct-10: Pointwise now writes the wedge/prism element out with VTK order
// for the nodes.

enum USGCell_type {
    none = 0,
    // For 2D simulations.
    triangle = 1,
    quad = 2,
    polygon = 3,
    // For 3D simulations.
    tetra = 4,
    wedge = 5,
    hexahedron = 6,
    pyramid = 7,
    // For subsetting the grid, need lines
    line = 8
}
string[] cell_names = ["none", "triangle", "quad", "polygon",
                       "tetra", "wedge", "hexahedron", "pyramid",
                       "line"];

int[] vtk_element_types = [0, VTKElement.triangle, VTKElement.quad, VTKElement.polygon,
                           VTKElement.tetra, VTKElement.wedge, VTKElement.hexahedron,
                           VTKElement.pyramid, VTKElement.line];

@nogc
USGCell_type convert_cell_type(int vtk_element_type)
{
    switch (vtk_element_type) {
    case VTKElement.triangle: return USGCell_type.triangle;
    case VTKElement.quad: return USGCell_type.quad;
    case VTKElement.tetra: return USGCell_type.tetra;
    case VTKElement.wedge: return USGCell_type.wedge;
    case VTKElement.hexahedron: return USGCell_type.hexahedron;
    case VTKElement.pyramid: return USGCell_type.pyramid;
    case VTKElement.line: return USGCell_type.line;
    default:
        return USGCell_type.none;
    }
}

@nogc
int n_vertices_for_vtk_element(int vtk_element_type)
{
    switch (vtk_element_type) {
    case VTKElement.triangle: return 3;
    case VTKElement.quad: return 4;
    case VTKElement.tetra: return 4;
    case VTKElement.wedge: return 6;
    case VTKElement.hexahedron: return 8;
    case VTKElement.pyramid: return 5;
    case VTKElement.line: return 2;
    default:
        return 0;
    }
}

@nogc
USGCell_type cell_type_from_name(string name)
{
    switch (name) {
    case "none": return USGCell_type.none;
    case "triangle": return USGCell_type.triangle;
    case "quad": return USGCell_type.quad;
    // case "polygon": return USGCell_type.polygon; // don't allow polygon cells
    case "tetra": return USGCell_type.tetra;
    case "wedge": return USGCell_type.wedge;
    case "hexahedron": return USGCell_type.hexahedron;
    case "pyramid": return USGCell_type.pyramid;
    case "line": return USGCell_type.line;
    default:
        return USGCell_type.none;
    }
}

string makeFaceTag(const size_t[] vtx_id_list)
{
    // We make a tag for this face out of the vertex id numbers
    // sorted in ascending order, to be used as a key in an
    // associative array of indices.  Any cycle of the same vertices
    // should define the same face, so we will use this tag as
    // a way to check if we have made a particular face before.
    // Of course, permutations of the same values that are not
    // correct cycles will produce the same tag, however,
    // we'll not worry about that for now because we should only
    // ever be providing correct cycles.
    size_t[] my_id_list = vtx_id_list.dup();
    sort(my_id_list);
    string tag = "";
    size_t n = my_id_list.length;
    foreach(i; 0 .. n) {
        if (i > 0) { tag ~= "-"; }
        tag ~= format("%d", my_id_list[i]);
    }
    return tag;
} // end makeFaceTag()

bool sameOrientation(const size_t[] vtx_id_list_a, const size_t[] vtx_id_list_b)
{
    assert(vtx_id_list_a.length == vtx_id_list_b.length, "Mismatch in vtx list lengths.");
    assert(isPermutation(vtx_id_list_a, vtx_id_list_b), "Lists should be permutations of vtx ids");
    if (vtx_id_list_a.length == 2) {
        // In 2D, faces have only 2 points, so they must be in order.
        return equal(vtx_id_list_a, vtx_id_list_b);
    } else {
        // In 3D, the cycle of points may start with any of them,
        // so we have to check them all.
        size_t[] test_list = vtx_id_list_b.dup();
        foreach (j; 0 .. test_list.length) {
            if (equal(vtx_id_list_a, test_list)) { return true; }
            bringToFront(test_list[1..$], test_list[$-1..$]); // for the next cycle
        }
        // At this point, we have checked all cycles and none match.
        return false;
    }
} // end sameOrientation()

class USGFace {
public:
    size_t[] vtx_id_list;
    string tag;
    // Note that, when construction a grid from a set of vertices and lists
    // of vertex indices, we will be able to tell if a face in a boundary list
    // is facing in or out by whether it has a right or left cell.
    // This information is filled in when constructing a new mesh but is not copied
    // in the USGFace copy constructor.  Beware.
    USGCell left_cell;
    USGCell right_cell;
    bool is_on_boundary = false; // presume not on boundary until marked as such

    this(const size_t[] vtx_id_list)
    {
        this.vtx_id_list = vtx_id_list.dup();
    }

    this(const USGFace other)
    {
        vtx_id_list = other.vtx_id_list.dup();
        // The const constraint means that we cannot copy the left_cell and
        // right_cell references.  The compiler cannot be sure that we will
        // honour the const promise once it has allowed us to copy references.
    }

    this(string str)
    // Defines the format for input from a text stream.
    {
        auto tokens = str.strip().split();
        vtx_id_list.length = to!int(tokens[0]);
        foreach(i; 0 .. vtx_id_list.length) { vtx_id_list[i] = to!int(tokens[i+1]); }
    }

    this(ref File fin)
    // constructor from raw binary file content.
    {
        int[1] buf1;
        fin.rawRead(buf1); vtx_id_list.length = buf1[0];
        foreach(i; 0 .. vtx_id_list.length) { fin.rawRead(buf1); vtx_id_list[i] = buf1[0]; }
    }

    override string toString()
    {
        return format("USGFace(vtx_id_list=%s)", vtx_id_list);
    }

    string toIOString()
    // Defines the format for output to a text stream.
    {
        string str = to!string(vtx_id_list.length);
        foreach (vtx_id; vtx_id_list) { str ~= " " ~ to!string(vtx_id); }
        return str;
    }

    void rawWrite(ref File fout)
    {
        int[1] buf1;
        buf1 = to!int(vtx_id_list.length); fout.rawWrite(buf1);
        foreach (vtx_id; vtx_id_list) { buf1 = to!int(vtx_id); fout.rawWrite(buf1); }
        return;
    }
} // end class USGFace

class USGCell {
public:
    USGCell_type cell_type;
    size_t id;
    size_t[] vtx_id_list;
    size_t[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this(const USGCell_type cell_type, size_t id, const size_t[] vtx_id_list,
         const size_t[] face_id_list, const int[] outsign_list)
    {
        this.id = id;
        this.cell_type = cell_type;
        this.vtx_id_list = vtx_id_list.dup();
        this.face_id_list = face_id_list.dup();
        this.outsign_list = outsign_list.dup();
    }

    this(const USGCell other, size_t id)
    {
        id = other.id;
        cell_type = other.cell_type;
        vtx_id_list = other.vtx_id_list.dup();
        face_id_list = other.face_id_list.dup();
        outsign_list = other.outsign_list.dup();
    }

    this(string str, size_t id)
    // Defines the format for input from a text stream.
    {
        this.id = id;
        auto tokens = str.strip().split();
        int itok = 0;
        cell_type = cell_type_from_name(tokens[itok++]);
        if (cell_type == USGCell_type.none) {
            throw new GeometryException("Unexpected cell type from line: " ~ str);
        }
        auto junk = tokens[itok++]; // "vtx"
        vtx_id_list.length = to!int(tokens[itok++]);
        foreach(i; 0 .. vtx_id_list.length) {
            vtx_id_list[i] = to!int(tokens[itok++]);
        }
        junk = tokens[itok++]; // "faces"
        face_id_list.length = to!int(tokens[itok++]);
        foreach(i; 0 .. face_id_list.length) {
            face_id_list[i] = to!int(tokens[itok++]);
        }
        junk = tokens[itok++]; // "outsigns"
        outsign_list.length = to!int(tokens[itok++]);
        foreach(i; 0 .. outsign_list.length) {
            outsign_list[i] = to!int(tokens[itok++]);
        }
    } // end constructor from string

    this(ref File fin, size_t id)
    // constructor from raw binary file content
    {
        int[1] buf1;
        fin.rawRead(buf1);
        cell_type = to!USGCell_type(buf1[0]);
        if (cell_type == USGCell_type.none) {
            throw new GeometryException("Unexpected cell type: " ~ to!string(cell_type));
        }
        fin.rawRead(buf1); vtx_id_list.length = buf1[0];
        foreach(i; 0 .. vtx_id_list.length) { fin.rawRead(buf1); vtx_id_list[i] = buf1[0]; }
        fin.rawRead(buf1); face_id_list.length = buf1[0];
        foreach(i; 0 .. face_id_list.length) { fin.rawRead(buf1); face_id_list[i] = buf1[0]; }
        fin.rawRead(buf1); outsign_list.length = buf1[0];
        foreach(i; 0 .. outsign_list.length) { fin.rawRead(buf1);  outsign_list[i] = buf1[0]; }
    } // end constructor from raw binary file content

    override string toString()
    {
        string repr = "USGCell(";
        repr ~= format("id=%d cell_type=%s", id, cell_type);
        repr ~= format(", vtx_id_list=%s", vtx_id_list);
        repr ~= format(", face_id_list=%s", face_id_list);
        repr ~= format(", outsign_list=%s", outsign_list);
        repr ~= ")";
        return repr;
    }

    string toIOString()
    // Defines the format for output to a text stream.
    {
        string str = cell_names[cell_type];
        str ~= " vtx " ~ to!string(vtx_id_list.length);
        foreach (vtx_id; vtx_id_list) { str ~= " " ~ to!string(vtx_id); }
        str ~= " faces " ~ to!string(face_id_list.length);
        foreach (face_id; face_id_list) { str ~= " " ~ to!string(face_id); }
        str ~= " outsigns " ~ to!string(outsign_list.length);
        foreach (outsign; outsign_list) { str ~= " " ~ to!string(outsign); }
        return str;
    }

    void rawWrite(ref File fout)
    {
        int[1] buf1;
        buf1 = to!int(cell_type); fout.rawWrite(buf1);
        buf1 = to!int(vtx_id_list.length); fout.rawWrite(buf1);
        foreach (vtx_id; vtx_id_list) { buf1 = to!int(vtx_id); fout.rawWrite(buf1); }
        buf1 = to!int(face_id_list.length); fout.rawWrite(buf1);
        foreach (face_id; face_id_list) { buf1 = to!int(face_id); fout.rawWrite(buf1); }
        buf1 = to!int(outsign_list.length); fout.rawWrite(buf1);
        foreach (outsign; outsign_list) { buf1 = to!int(outsign); fout.rawWrite(buf1); }
        return;
    }
} // end class USGCell

class BoundaryFaceSet {
public:
    string tag;
    size_t[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this(string tag, size_t[] face_id_list, int[] outsign_list)
    {
        this.tag = tag;
        this.face_id_list = face_id_list.dup();
        this.outsign_list = outsign_list.dup();
    } // end constructor

    this(string str)
    // Defines the format for input from a text stream.
    {
        auto tokens = str.strip().split();
        if (tokens.length == 1) {
            // Special case; just the string tag.
            this.tag = tokens[0];
        } else {
            // We have more information so try to use it.
            int itok = 0;
            this.tag = tokens[itok++];
            auto junk = tokens[itok++]; // "faces"
            face_id_list.length = to!int(tokens[itok++]);
            foreach(i; 0 .. face_id_list.length) {
                face_id_list[i] = to!int(tokens[itok++]);
            }
            junk = tokens[itok++]; // "outsigns"
            outsign_list.length = to!int(tokens[itok++]);
            foreach(i; 0 .. outsign_list.length) {
                outsign_list[i] = to!int(tokens[itok++]);
            }
            assert(face_id_list.length == outsign_list.length,
                   "Mismatch in numbers of faces and outsigns");
        }
    } // end constructor from a string

    this(ref File fin)
    {
        int[1] buf1; fin.rawRead(buf1);
        int tag_length = buf1[0];
        if (tag_length > 0) {
            char[] found_tag = new char[tag_length];
            fin.rawRead(found_tag);
            tag = to!string(found_tag);
        }
        fin.rawRead(buf1); face_id_list.length = buf1[0];
        foreach(i; 0 .. face_id_list.length) { fin.rawRead(buf1); face_id_list[i] = buf1[0]; }
        fin.rawRead(buf1); outsign_list.length = buf1[0];
        foreach(i; 0 .. outsign_list.length) { fin.rawRead(buf1);  outsign_list[i] = buf1[0]; }
    } // end constructor from raw file content

    this(const BoundaryFaceSet other)
    {
        tag = other.tag;
        face_id_list = other.face_id_list.dup();
        outsign_list = other.outsign_list.dup();
    } // end constructor from another

    override string toString()
    {
        string repr = "BoundaryFaceSet(";
        repr ~= format("tag=\"%s\"", tag);
        repr ~= format(", face_id_list=[%d ... %d]", face_id_list[0], face_id_list[$-1]);
        repr ~= format(", outsign_list=[%d ... %d]", outsign_list[0], outsign_list[$-1]);
        repr ~= ")";
        return repr;
    }

    string toIOString()
    // Defines the format for output to a text stream.
    {
        string str = tag;
        str ~= " faces " ~ to!string(face_id_list.length);
        foreach (face_id; face_id_list) { str ~= " " ~ to!string(face_id); }
        str ~= " outsigns " ~ to!string(outsign_list.length);
        foreach (outsign; outsign_list) { str ~= " " ~ to!string(outsign); }
        return str;
    }

    void rawWrite(ref File fout)
    {
        int[1] buf1;
        buf1 = to!int(tag.length); fout.rawWrite(buf1);
        fout.rawWrite(to!(char[])(tag));
        buf1 = to!int(face_id_list.length); fout.rawWrite(buf1);
        foreach (face_id; face_id_list) { buf1 = to!int(face_id); fout.rawWrite(buf1); }
        buf1 = to!int(outsign_list.length); fout.rawWrite(buf1);
        foreach (outsign; outsign_list) { buf1 = to!int(outsign); fout.rawWrite(buf1); }
        return;
    }
} // end class BoundaryFaceSet


class UnstructuredGrid : Grid {
public:
    size_t nfaces, nboundaries, ninteriorfaces;
    USGFace[] faces;
    // The following dictionary of faces, identified by their defining vertices,
    // will enable us to quickly search for an already defined face.
    size_t[string] faceIndices;
    USGCell[] cells;
    BoundaryFaceSet[] boundaries;

    this(int dimensions, string label="")
    // A new empty grid that will have its details filled in later.
    {
        super(Grid_t.unstructured_grid, dimensions, label);
    }

    this(const UnstructuredGrid other, const string new_label="")
    {
        super(Grid_t.unstructured_grid, other.dimensions,
              ((new_label == "")? other.label : new_label));
        nvertices = other.nvertices;
        ncells = other.ncells;
        nfaces = other.nfaces;
        nboundaries = other.nboundaries;
        foreach(v; other.vertices) { vertices ~= Vector3(v); }
        foreach(f; other.faces) {
            faceIndices[makeFaceTag(f.vtx_id_list)] = faces.length;
            faces ~= new USGFace(f);
        }
        foreach(id, c; other.cells) { cells ~= new USGCell(c, id); }
        foreach(b; other.boundaries) { boundaries ~= new BoundaryFaceSet(b); }
        niv = vertices.length; njv = 1; nkv = 1;
    } // end constructor from another

    this(const Vector3[][] bPoints, const string new_label="", double tol=1.0e-9)
    // Construction via paving with (mostly) quadrilateral cells.
    // Input is one or more arrays of points that define a closed region.
    {
        super(Grid_t.unstructured_grid, 2, new_label);
        // Copy the boundary boints into the grid's vertex list,
        // making sure that the start-point of each new boundary-vertex list
        // coincides with the end-point of the previous boundary-vertex list.
        if (distance_between(bPoints[0][0], bPoints[$-1][$-1]) > tol) {
            string msg = "Beginning point and end point do not coincide: " ~
                to!string(bPoints[0][0]) ~ " " ~ to!string(bPoints[$-1][$-1]);
            throw new Error(msg);
        }
        size_t[] bndry_vtx_ids;
        foreach (i, bpa; bPoints) {
            size_t[] face_id_list;
            int[] outsign_list;
            if (i == 0) {
                vertices ~= Vector3(bpa[0]);
                bndry_vtx_ids ~= vertices.length - 1;
            }
            foreach (j; 1 .. bpa.length) {
                vertices ~= Vector3(bpa[j]);
                bndry_vtx_ids ~= vertices.length - 1;
                // Also, build a face object on the boundary.
                size_t[] vtx_id_list = [vertices.length-2, vertices.length-1];
                string faceTag = makeFaceTag(vtx_id_list);
                faceIndices[faceTag] = faces.length; // before new face
                faces ~= new USGFace(vtx_id_list);
                face_id_list ~= faces.length-1;
                outsign_list ~= 1; // always pointing out for a counter-clockwise vertex order
            }
            boundaries ~= new BoundaryFaceSet("bndry-"~to!string(i), face_id_list, outsign_list);
            if (i > 0) {
                if (distance_between(bPoints[i][0], bPoints[i-1][$-1]) > tol) {
                    string msg = "Beginning point for segment[" ~ to!string(i) ~ "] and end point " ~
                        "for segment [" ~ to!string(i-1) ~ "] do not coincide: " ~
                        to!string(bPoints[i][0]) ~ " " ~ to!string(bPoints[i-1][$-1]);
                    throw new Error(msg);
                }
            }
        }
        // At this point, we should have a closed region defined by its bounding polygon.
        // Need to fill in interior points, faces and cells.
        paver2d.fill_interior(vertices, faces, faceIndices, cells, bndry_vtx_ids);
        // Should also be able to use any other algorithm, such as an advancing-front method.
        //
        // Fill in some other properties that the flow-code expects.
        niv = vertices.length; njv = 1; nkv = 1;
        nvertices = vertices.length;
        nfaces = faces.length;
        ncells = cells.length;
        nboundaries = boundaries.length;
        sort_faces_in_boundary_order();
    } // end construction via paving in 2D

    this(const Vector3[] boundary, BoundaryFaceSet[] in_boundaries, const string new_label="")
    // Paved Grid Constructor by Heather Muir, 2016.
    // Deprecated in 2018. We will use the constructor above.
    {
        double[][] boundary_points;
        foreach(p; boundary){
            boundary_points ~= [p.x.re, p.y.re, p.z.re];
        }
        PavedGrid grid = new PavedGrid(boundary_points);
        super(Grid_t.unstructured_grid, grid.dimensions,
              ((new_label == "")? grid.label : new_label));

        foreach(b; in_boundaries) { boundaries ~= new BoundaryFaceSet(b); }
        this.nvertices = grid.nvertices;
        this.ncells = grid.ncells;
        this.nfaces = grid.nfaces;
        this.nboundaries = grid.nboundaries;
        foreach(p; POINT_LIST){
            double[] v = [p.x, p.y, p.z];
            this.vertices ~= Vector3(v);
        }
        foreach(f; FACE_LIST){
            size_t[] vtx_id_list = f.point_IDs;
            faceIndices[makeFaceTag(vtx_id_list)] = faces.length;
            this.faces ~= new USGFace(vtx_id_list);
        }
        foreach(id, c; CELL_LIST){
            c.auto_cell_type();
            if(c.cell_type != "quad"){
                throw new Error(text("paver generated a non-quad cell"));
            } else {
                USGCell_type cell_type = cell_type_from_name(c.cell_type);
                this.cells ~= new USGCell(cell_type, id, c.point_IDs, c.face_IDs, c.outsigns);
            }
        }
        niv = vertices.length; njv = 1; nkv = 1;
        sort_faces_in_boundary_order();
    } //end Heather's paved grid constructor

    this(const StructuredGrid sg, const string new_label="")
    {
        string my_block_label = (new_label == "")? sg.label : new_label;
        super(Grid_t.unstructured_grid, sg.dimensions, my_block_label);
        vertices.length = 0;
        faces.length = 0;
        cells.length = 0;
        boundaries.length = 0;
        if (dimensions == 2) {
            nvertices = sg.niv * sg.njv;
            nfaces = (sg.niv)*(sg.njv-1) + (sg.niv-1)*(sg.njv);
            ncells = (sg.niv-1)*(sg.njv-1);
            nboundaries = 4;
            foreach(ib; 0 .. nboundaries) {
                string boundary_tag;
                if (sg.bcTags[ib].length > 0) {
                    // Copy the original tag string.
                    boundary_tag = sg.bcTags[ib];
                } else {
                    // Fabricate a tag string.
                    if (my_block_label.length > 0) { boundary_tag = my_block_label ~ "-"; }
                    boundary_tag ~= face_name[ib]; // 0=west, 1=east, 2=south, 3=north
                }
                boundaries ~= new BoundaryFaceSet(boundary_tag);
            }
            // vertex index array
            size_t[][] vtx_id;
            vtx_id.length = sg.niv;
            foreach (i; 0 .. sg.niv) {
                vtx_id[i].length = sg.njv;
            }
            // vertex list in standard order, indexed
            foreach (j; 0 .. sg.njv) {
                foreach (i; 0 .. sg.niv) {
                    vertices ~= Vector3(sg.vertices[sg.single_index(i,j)]);
                    vtx_id[i][j] = vertices.length - 1;
                }
            }
            // i-faces
            // As we move from the first vertex to the second vertex,
            // the unit normal will be pointing to the right.
            size_t[][] iface_id;
            iface_id.length = sg.niv;
            foreach (i; 0 .. sg.niv) {
                iface_id[i].length = sg.njv-1;
            }
            foreach (j; 0 .. sg.njv-1) {
                foreach (i; 0 .. sg.niv) {
                    size_t[] vtx_id_list = [vtx_id[i][j], vtx_id[i][j+1]];
                    faceIndices[makeFaceTag(vtx_id_list)] = faces.length;
                    faces ~= new USGFace(vtx_id_list);
                    iface_id[i][j] = faces.length - 1;
                    if (i == 0) {
                        boundaries[Face.west].face_id_list ~= iface_id[i][j];
                        boundaries[Face.west].outsign_list ~= -1;
                    }
                    if (i == sg.niv-1) {
                        boundaries[Face.east].face_id_list ~= iface_id[i][j];
                        boundaries[Face.east].outsign_list ~= +1;
                    }
                }
            }
            // j-faces
            size_t[][] jface_id;
            jface_id.length = sg.niv - 1;
            foreach (i; 0 .. sg.niv-1) {
                jface_id[i].length = sg.njv;
            }
            foreach (j; 0 .. sg.njv) {
                foreach (i; 0 .. sg.niv-1) {
                    size_t[] vtx_id_list = [vtx_id[i+1][j], vtx_id[i][j]];
                    faceIndices[makeFaceTag(vtx_id_list)] = faces.length;
                    faces ~= new USGFace(vtx_id_list);
                    jface_id[i][j] = faces.length - 1;
                    if (j == 0) {
                        boundaries[Face.south].face_id_list ~= jface_id[i][j];
                        boundaries[Face.south].outsign_list ~= -1;
                    }
                    if (j == sg.njv-1) {
                        boundaries[Face.north].face_id_list ~= jface_id[i][j];
                        boundaries[Face.north].outsign_list ~= +1;
                    }
                }
            }
            // cells
            foreach (j; 0 .. sg.njv-1) {
                foreach (i; 0 .. sg.niv-1) {
                    auto cell_vertices = [vtx_id[i][j], vtx_id[i+1][j],
                                          vtx_id[i+1][j+1], vtx_id[i][j+1]];
                    // For faces, work clockwise around the cell.
                    // This is likely to give us nearly orthogonal faces 0 and 1,
                    // something that is assumed by the CFL signal_speed calcaulation.
                    auto cell_faces = [jface_id[i][j+1], // north
                                       iface_id[i+1][j], // east
                                       jface_id[i][j], // south
                                       iface_id[i][j]]; // west
                    auto outsigns = [+1, +1, -1, -1];
                    size_t id = cells.length;
                    auto my_cell = new USGCell(USGCell_type.quad, id, cell_vertices,
                                               cell_faces, outsigns);
                    cells ~= my_cell;
                    // Now that we have the new cell, we can make the connections
                    // from the faces back to the cell.
                    faces[iface_id[i][j]].right_cell = my_cell;
                    faces[iface_id[i+1][j]].left_cell = my_cell;
                    faces[jface_id[i][j]].right_cell = my_cell;
                    faces[jface_id[i][j+1]].left_cell = my_cell;
                }
            }
        } else {
            // Assume dimensions == 3
            nvertices = sg.niv * sg.njv * sg.nkv;
            nfaces = (sg.niv)*(sg.njv-1)*(sg.nkv-1) +
                (sg.niv-1)*(sg.njv)*(sg.nkv-1) +
                (sg.niv-1)*(sg.njv-1)*(sg.nkv);
            ncells = (sg.niv-1)*(sg.njv-1)*(sg.nkv-1);
            nboundaries = 6;
            foreach(ib; 0 .. nboundaries) {
                string boundary_tag;
                if (sg.bcTags[ib].length > 0) {
                    // Copy the original tag string.
                    boundary_tag = sg.bcTags[ib];
                } else {
                    // Fabricate a tag string.
                    if (my_block_label.length > 0) { boundary_tag = my_block_label ~ "-"; }
                    boundary_tag ~= face_name[ib]; // 0=west, 1=east, 2=south, 3=north, 4=bottom, 5=top
                }
                boundaries ~= new BoundaryFaceSet(boundary_tag);
            }
            // vertex index array
            size_t[][][] vtx_id;
            vtx_id.length = sg.niv;
            foreach (i; 0 .. sg.niv) {
                vtx_id[i].length = sg.njv;
                foreach (j; 0 .. sg.njv) {
                    vtx_id[i][j].length = sg.nkv;
                }
            }
            // vertex list in standard order, indexed
            foreach (k; 0 .. sg.nkv) {
                foreach (j; 0 .. sg.njv) {
                    foreach (i; 0 .. sg.niv) {
                        vertices ~= Vector3(sg.vertices[sg.single_index(i,j,k)]);
                        vtx_id[i][j][k] = vertices.length - 1;
                    }
                }
            }
            // i-faces
            // Using the right-hand rule:
            // as we progress around the vertices, following the curl of the fingers,
            // the right thumb is pointing in the direction of the unit normal.
            size_t[][][] iface_id;
            iface_id.length = sg.niv;
            foreach (i; 0 .. sg.niv) {
                iface_id[i].length = sg.njv-1;
                foreach (j; 0 .. sg.njv-1) {
                    iface_id[i][j].length = sg.nkv-1;
                }
            }
            foreach (k; 0 .. sg.nkv-1) {
                foreach (j; 0 .. sg.njv-1) {
                    foreach (i; 0 .. sg.niv) {
                        size_t[] vtx_id_list = [vtx_id[i][j][k], vtx_id[i][j+1][k],
                                                vtx_id[i][j+1][k+1], vtx_id[i][j][k+1]];
                        faceIndices[makeFaceTag(vtx_id_list)] = faces.length;
                        faces ~= new USGFace(vtx_id_list);
                        iface_id[i][j][k] = faces.length - 1;
                        if (i == 0) {
                            boundaries[Face.west].face_id_list ~= iface_id[i][j][k];
                            boundaries[Face.west].outsign_list ~= -1;
                        }
                        if (i == sg.niv-1) {
                            boundaries[Face.east].face_id_list ~= iface_id[i][j][k];
                            boundaries[Face.east].outsign_list ~= +1;
                        }
                    }
                }
            }
            // j-faces
            size_t[][][] jface_id;
            jface_id.length = sg.niv-1;
            foreach (i; 0 .. sg.niv-1) {
                jface_id[i].length = sg.njv;
                foreach (j; 0 .. sg.njv) {
                    jface_id[i][j].length = sg.nkv-1;
                }
            }
            foreach (k; 0 .. sg.nkv-1) {
                foreach (j; 0 .. sg.njv) {
                    foreach (i; 0 .. sg.niv-1) {
                        size_t[] vtx_id_list = [vtx_id[i+1][j][k], vtx_id[i][j][k],
                                                vtx_id[i][j][k+1], vtx_id[i+1][j][k+1]];
                        faceIndices[makeFaceTag(vtx_id_list)] = faces.length;
                        faces ~= new USGFace(vtx_id_list);
                        jface_id[i][j][k] = faces.length - 1;
                        if (j == 0) {
                            boundaries[Face.south].face_id_list ~= jface_id[i][j][k];
                            boundaries[Face.south].outsign_list ~= -1;
                        }
                        if (j == sg.njv-1) {
                            boundaries[Face.north].face_id_list ~= jface_id[i][j][k];
                            boundaries[Face.north].outsign_list ~= +1;
                        }
                    }
                }
            }
            // k-faces
            size_t[][][] kface_id;
            kface_id.length = sg.niv-1;
            foreach (i; 0 .. sg.niv-1) {
                kface_id[i].length = sg.njv-1;
                foreach (j; 0 .. sg.njv-1) {
                    kface_id[i][j].length = sg.nkv;
                }
            }
            foreach (k; 0 .. sg.nkv) {
                foreach (j; 0 .. sg.njv-1) {
                    foreach (i; 0 .. sg.niv-1) {
                        size_t[] vtx_id_list = [vtx_id[i][j][k], vtx_id[i+1][j][k],
                                                vtx_id[i+1][j+1][k], vtx_id[i][j+1][k]];
                        faceIndices[makeFaceTag(vtx_id_list)] = faces.length;
                        faces ~= new USGFace(vtx_id_list);
                        kface_id[i][j][k] = faces.length - 1;
                        if (k == 0) {
                            boundaries[Face.bottom].face_id_list ~= kface_id[i][j][k];
                            boundaries[Face.bottom].outsign_list ~= -1;
                        }
                        if (k == sg.nkv-1) {
                            boundaries[Face.top].face_id_list ~= kface_id[i][j][k];
                            boundaries[Face.top].outsign_list ~= +1;
                        }
                    }
                }
            }
            // cells
            foreach (k; 0 .. sg.nkv-1) {
                foreach (j; 0 .. sg.njv-1) {
                    foreach (i; 0 .. sg.niv-1) {
                        auto cell_vertices = [vtx_id[i][j][k], vtx_id[i+1][j][k],
                                              vtx_id[i+1][j+1][k], vtx_id[i][j+1][k],
                                              vtx_id[i][j][k+1], vtx_id[i+1][j][k+1],
                                              vtx_id[i+1][j+1][k+1], vtx_id[i][j+1][k+1]];
                        // Keep the old face order,
                        // so that the CFL signal_speed calculation picks
                        // a good set of faces, at least for quad or hex meshes.
                        auto cell_faces = [jface_id[i][j+1][k], // north
                                           iface_id[i+1][j][k], // east
                                           jface_id[i][j][k], // south
                                           iface_id[i][j][k], // west
                                           kface_id[i][j][k+1], // top
                                           kface_id[i][j][k]]; // bottom
                        auto outsigns = [+1, +1, -1, -1, +1, -1];
                        size_t id = cells.length;
                        auto my_cell = new USGCell(USGCell_type.hexahedron, id, cell_vertices,
                                                   cell_faces, outsigns);
                        cells ~= my_cell;
                        // Now that we have the new cell, we can make the connections
                        // from the faces back to the cell.
                        faces[iface_id[i][j][k]].right_cell = my_cell;
                        faces[iface_id[i+1][j][k]].left_cell = my_cell;
                        faces[jface_id[i][j][k]].right_cell = my_cell;
                        faces[jface_id[i][j+1][k]].left_cell = my_cell;
                        faces[kface_id[i][j][k]].right_cell = my_cell;
                        faces[kface_id[i][j][k+1]].left_cell = my_cell;
                    }
                }
            }
        } // end if dimensions
        assert(nvertices == vertices.length, "mismatch in number of vertices");
        assert(nfaces == faces.length, "mismatch in number of faces");
        assert(ncells == cells.length, "mismatch in number of cells");
        niv = vertices.length; njv = 1; nkv = 1;

        sort_faces_in_boundary_order();
    } // end constructor from StructuredGrid object

    // Imported grid.
    this(string fileName, string fmt, bool true_centroids, double scale=1.0,
         bool expect_gmsh_order_for_wedges=false, string new_label="")
    {
        super(Grid_t.unstructured_grid, 0, new_label);
        // dimensions will be reset on reading grid
        switch (fmt) {
        case "gziptext": read_from_gzip_file(fileName, scale); break;
        case "rawbinary": read_from_raw_binary_file(fileName, scale); break;
        case "su2text": read_from_su2_file(fileName, true_centroids, scale, expect_gmsh_order_for_wedges); break;
        case "vtktext": read_from_vtk_text_file(fileName, scale); break;
        case "vtkxml": throw new Error("Reading from VTK XML format not implemented.");
        default: throw new Error("Import an UnstructuredGrid, unknown format: " ~ fmt);
        }
        if (new_label != "") { label = new_label; }
        niv = vertices.length; njv = 1; nkv = 1;
        sort_faces_in_boundary_order();
    } // end constructor from data imported from a file

    UnstructuredGrid dup() const
    {
        return new UnstructuredGrid(this);
    }


    // -----------------------------
    // Indexing and location methods.
    // -----------------------------
    @nogc
    override Vector3* opIndex(size_t i, size_t j, size_t k=0)
    in {
        assert (i < nvertices, "index i is invalid");
        assert (j == 0, "index j is invalid");
        assert (k == 0, "index k is invalid");
    }
    do {
        return &(vertices[i]);
    }

    @nogc
    override Vector3* opIndex(size_t indx)
    in { assert (indx < nvertices, "index indx is invalid"); }
    do {
        return &(vertices[indx]);
    }

    override string toString()
    {
        string repr = "UnstructuredGrid(";
        repr ~= format("dimensions=%d, nvertices=%d, nfaces=%d, ncells=%d",
                       dimensions, nvertices, nfaces, ncells);
        repr ~= format(", vertices=[%s ... %s]", vertices[0], vertices[$-1]);
        repr ~= format(", faces=[%s ... %s]", faces[0], faces[$-1]);
        repr ~= format(", cells=[%s ... %s]", cells[0], cells[$-1]);
        repr ~= format(", boundaries=%s", boundaries);
        repr ~= ")";
        return repr;
    }

    void sort_faces_in_boundary_order()
    /*
        This routine enables some optimisation trickery later on in the flow solver.
        Specifically, when looping through the faces, we want to be able to know where
        the interior ones are without having to manually check interiorness every time.

        We do this by sorting the face references so that all the interior ones are first,
        then the boundary faces grouped together in the order their boundaries appear in
        in "boundaries" the array of BoundaryFaceSet objects.

        @author: Nick Gibbons (13/12/22)
    */
    {

        //writefln("Called sort_faces_in_boundary_order... ");
        //SysTime start = Clock.currTime();

        // First thing we do is build look up tables for finding which boundary set and position
        // each boundary face belongs to. We use a pair of temporary associative arrays for this.
        size_t[size_t] faceidx_to_boundary_table;
        size_t[size_t] faceidx_to_boundary_idx_table;


        size_t nboundaryfaces = 0;
        foreach(i, bfs; boundaries){
            foreach(j, faceid; bfs.face_id_list){
                faces[faceid].is_on_boundary = true;
                faceidx_to_boundary_table[faceid] = i;
                faceidx_to_boundary_idx_table[faceid] = j;
                nboundaryfaces += 1;
            }
        }

        ninteriorfaces = nfaces-nboundaryfaces;

        size_t nswaps = 0;
        size_t j = nfaces-1;
        // Now we do the sort, in place, by swapping references (see 14/12/22 notes for derivation)
        foreach(i; 0 .. ninteriorfaces){
            // If face at i is internal, that's good. We can skip an iteration.
            if (!faces[i].is_on_boundary) continue;

            // At this point we know that the i face is on a boundary. We want to swap
            // it with an internal face, so scan along until we find one:
            // In v2 of this algorithm, we actually start at the far end.
            while (faces[j].is_on_boundary) {
                j -= 1;
            }
            // We've found an internal face in position j. This needs to be swapped
            // in the faces array, but also anywhere a face index is recorded. There
            // are four places in the UnstructuredGrid object where this occurs.

            // 1.) the faceIndices arrays, which map vertex ids to face ids
            string ifaceTag = makeFaceTag(faces[i].vtx_id_list);
            string jfaceTag = makeFaceTag(faces[j].vtx_id_list);
            faceIndices[ifaceTag] = j;
            faceIndices[jfaceTag] = i;

            // 2.) The USGCell objects, which have a list of face ids attached to them
            // Note that i's cells may not be initialised, which we need to test for.
            size_t il, ir, jl, jr;
            if (faces[i].left_cell)  il = countUntil(faces[i].left_cell.face_id_list, i);
            if (faces[i].right_cell) ir = countUntil(faces[i].right_cell.face_id_list, i);
            if (faces[j].left_cell)  jl = countUntil(faces[j].left_cell.face_id_list, j);
            if (faces[j].right_cell) jr = countUntil(faces[j].right_cell.face_id_list, j);

            if (faces[i].left_cell)  faces[i].left_cell.face_id_list[il] = j;
            if (faces[i].right_cell) faces[i].right_cell.face_id_list[ir] = j;
            if (faces[j].left_cell)  faces[j].left_cell.face_id_list[jl] = i;
            if (faces[j].right_cell) faces[j].right_cell.face_id_list[jr] = i;

            // 3.) The BoundaryFaceSet collection called boundaries. We use the tables
            // we made eariler to figure out which one i is sitting inside.
            size_t ibfs = faceidx_to_boundary_table[i];
            size_t ibfsidx = faceidx_to_boundary_idx_table[i];
            assert(boundaries[ibfs].face_id_list[ibfsidx] == i);
            boundaries[ibfs].face_id_list[ibfsidx] = j;

            // 4.) This face might come up again later. That means we need to update the
            // BFS lookup tables to reflect the new index. The old key at i is no longer valid
            // but we should never be looking that up in future iterations, because i is always
            // increasing.
            faceidx_to_boundary_table[j] = ibfs;
            faceidx_to_boundary_idx_table[j] = ibfsidx;

            // 5.) Finally, swap the face references in memory
            swap(faces[i], faces[j]);
            nswaps += 1;
        } // loop through i interior worth of faces

        // Check that this was done properly...
        //foreach(i, face; faces) {
        //    writefln(" %d %s %s", i, face.vtx_id_list, face.is_on_boundary);
        //    writefln(" %s %d ", makeFaceTag(face.vtx_id_list), faceIndices[makeFaceTag(face.vtx_id_list)]);
        //    if (face.left_cell)  writefln("   lcell %s" , face.left_cell.face_id_list); // checked
        //    if (face.right_cell) writefln("   rcell %s ", face.right_cell.face_id_list); // checked
        //}
        //foreach(i, bfs; boundaries){
        //    writefln(" %d %s %s", i, bfs.tag, bfs.face_id_list);
        //}

        //double elapsed = to!double((Clock.currTime() - start).total!"msecs"());
        //writefln("Done %d swaps in %f msec", nswaps, elapsed);
    }

    @nogc
    override size_t number_of_vertices_for_cell(size_t i)
    {
        return cells[i].vtx_id_list.length;
    }

    @nogc
    override int vtk_element_type_for_cell(size_t i)
    {
        return vtk_element_types[cells[i].cell_type];
    }

    @nogc
    override int get_cell_type(size_t i)
    {
        return cells[i].cell_type;
    }

    override size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const
    in {
        assert (i < ncells, text("index i=", i, " is invalid, ncells=", ncells));
        assert (j == 0, text("index j=", j, " is invalid for unstructured grid"));
        assert (k == 0, text("index k=", k, " is invalid for unstructured grid"));
    }
    do {
        return cells[i].vtx_id_list.dup();
    }

    @nogc
    override void copy_vtx_id_list_for_cell(ref size_t[8] vtx_list_copy, ref size_t nvtx,
                                            size_t indx, size_t j, size_t k=0) const
    {
        copy_vtx_id_list_for_cell(vtx_list_copy, nvtx, indx);
    }

    override size_t[] get_vtx_id_list_for_cell(size_t indx) const
    in {
        assert (indx < ncells,
                text("index indx=", indx, " is invalid, ncells=", ncells));
    }
    do {
        return cells[indx].vtx_id_list.dup();
    }

    @nogc
    override void copy_vtx_id_list_for_cell(ref size_t[8] vtx_list_copy, ref size_t nvtx,
                                            size_t indx) const
    {
        nvtx = cells[indx].vtx_id_list.length;
        foreach(i; 0 .. nvtx) { vtx_list_copy[i] = cells[indx].vtx_id_list[i]; }
    }

    size_t[] get_connected_cell_ids_for_cell(size_t i)
    {
        bool[size_t] connectedCells;
        foreach (iFace; cells[i].face_id_list) {
            if (faces[iFace].left_cell)
                connectedCells[faces[iFace].left_cell.id] = true;
            if (faces[iFace].right_cell)
                connectedCells[faces[iFace].right_cell.id] = true;
        }
        return connectedCells.keys.dup;
    }

    override Grid get_boundary_grid(size_t boundary_indx)
    // Returns the grid defining a particular boundary of the original grid.
    // For an 3D block, a 2D surface grid will be returned, with index directions
    // as defined on the debugging cube.
    // This new grid is intended for use just in the output of a subset
    // of the simulation data and not for actual simulation.
    // It does not carry all of the information required for flow simulation.
    {
        int new_dimensions = dimensions - 1;
        if (new_dimensions < 1 || new_dimensions > 2) {
            throw new GeometryException(format("Invalid new_dimensions=%d", new_dimensions));
        }
        BoundaryFaceSet bfs = boundaries[boundary_indx];
        UnstructuredGrid new_grid = new UnstructuredGrid(new_dimensions, bfs.tag);
        // We'll make a full copy of the original grid's vertices
        // so that we don't have to renumber them.
        new_grid.nvertices = nvertices;
        new_grid.vertices.length = nvertices;
        foreach (i; 0 .. nvertices) { new_grid.vertices[i].set(vertices[i]); }
        foreach (fid; bfs.face_id_list) {
            // fid is the face in the original grid.
            // It will become the cell in the new grid.
            USGCell_type new_cell_type = USGCell_type.none;
            switch (faces[fid].vtx_id_list.length) {
            case 2:
                assert(dimensions == 2, "Assumed 2D but it was not so.");
                new_cell_type = USGCell_type.line;
                break;
            case 3:
                assert(dimensions == 3, "Assumed 3D but it was not so.");
                new_cell_type = USGCell_type.triangle;
                break;
            case 4:
                assert(dimensions == 3, "Assumed 3D but it was not so.");
                new_cell_type = USGCell_type.quad;
                break;
            default:
                throw new GeometryException("Did not know what to do with this number of vertices.");
            }
            size_t[] newCell_face_id_list; // empty
            int[] newCell_outsign_list; // empty
            size_t id = new_grid.cells.length;
            new_grid.cells ~= new USGCell(new_cell_type, id, faces[fid].vtx_id_list,
                                          newCell_face_id_list, newCell_outsign_list);
        }
        new_grid.ncells = new_grid.cells.length;
        return new_grid;
    } // end get_boundary_grid()

    override size_t[] get_list_of_boundary_cells(size_t boundary_indx)
    // Prepares list of cells indicies that match the boundary grid selected by
    // the function get_boundary_grid().  See above.
    {
        size_t new_nic, new_njc, new_nkc;
        size_t[] cellList;
        // Make a dictionary of cell indices, so that we may look them up
        // when building the list of indices.
        size_t[size_t] cell_indx;
        foreach (i, c; cells) { cell_indx[c.id] = i; }
        // Work through the set of boundary faces and identify the inside cell
        // as being the one we want to associate with the boundary face.
        BoundaryFaceSet bfs = boundaries[boundary_indx];
        foreach (i, fid; bfs.face_id_list) {
            USGCell insideCell = (bfs.outsign_list[i] == 1) ? faces[fid].left_cell : faces[fid].right_cell;
            if (insideCell) {
                cellList ~= cell_indx[insideCell.id];
            } else {
                string msg = format("Oops, while processing boundary faces, " ~
                                    "i=%d, insideCell is null for face.id=%d", i, fid);
                throw new Exception(msg);
            }
        }
        return cellList;
    } // end get_list_of_boundary_cells()

    // Functions for checking if defined regions are fully closed.
    Vector3 vectorAreaOfFacet(size_t[] vtx_id_list)
    {
        size_t nv = vtx_id_list.length;
        assert(nv > 1, "do not have a valid list of vertex ids");
        Vector3 varea = Vector3(0,0,0);
        switch (nv) {
        case 2:
            auto p0 = vertices[vtx_id_list[0]];
            auto p1 = vertices[vtx_id_list[1]];
            varea = cross(p1-p0, Vector3(0,0,1.0)); // unit depth in z for this 2D face
            break;
        case 3:
            auto p0 = vertices[vtx_id_list[0]];
            auto p1 = vertices[vtx_id_list[1]];
            auto p2 = vertices[vtx_id_list[2]];
            varea = 0.5 * cross(p1-p0, p2-p0);
            break;
        default:
            Vector3 pmid = Vector3(0,0,0);
            foreach (vid; vtx_id_list) { pmid += vertices[vid]; }
            pmid /= nv;
            // Compute the vector area is triangle sections.
            foreach (i; 0 .. nv-1) {
                auto p0 = vertices[vtx_id_list[i]];
                auto p1 = vertices[vtx_id_list[i+1]];
                varea += 0.5 * cross(p1-p0, pmid-p0);
            }
            auto p0 = vertices[vtx_id_list[nv-1]];
            auto p1 = vertices[vtx_id_list[0]];
            varea += 0.5 * cross(p1-p0, pmid-p0);
        }
        return varea;
    } // end vectorAreaOfFacet

    Vector3 vectorAreaOfCell(const USGCell c)
    {
        Vector3 varea = Vector3(0,0,0);
        foreach (j, fid; c.face_id_list) {
            varea += c.outsign_list[j] * vectorAreaOfFacet(faces[fid].vtx_id_list);
        }
        return varea;
    } // end vectorAreaOfCell

    // -----------------------
    // Cell insertion methods
    // -----------------------
    void splitTriangleCell(size_t id)
    {
        assert(id < cells.length);
        if (dimensions != 2) {
            string excMsg = "Exception in UnstructureGrid:splitTriangleCell():\n";
            excMsg ~= "This method can only be used on 2D grids.\n";
            throw new GeometryException(excMsg);
        }
        if (cells[id].cell_type != USGCell_type.triangle) {
            string excMsg = "Exception in UnstructureGrid:splitTriangleCell():\n";
            excMsg ~= "This method can only be used on TRIANGLE cells.\n";
            throw new GeometryException(excMsg);
        }

        auto vtxList = get_vtx_id_list_for_cell(id, 0, 0);
        auto p0 = vertices[vtxList[0]];
        auto p1 = vertices[vtxList[1]];
        auto p2 = vertices[vtxList[2]];
        Vector3 centroid = Vector3(0.0, 0.0, 0.0);
        centroid.set((p0.x + p1.x + p2.x)/3.0,
                     (p0.y + p1.y + p2.y)/3.0,
                     (p0.z + p1.z + p2.z)/3.0);
        // Add centroid as new vertex in grid.
        vertices ~= centroid;
        nvertices = vertices.length;
        auto ctrId = vertices.length-1;
        // Add new faces
        faces ~= new USGFace([vtxList[0], ctrId]);
        auto if0_id = faces.length-1;
        faceIndices[makeFaceTag([vtxList[0], ctrId])] = if0_id;

        faces ~= new USGFace([vtxList[1], ctrId]);
        auto if1_id = faces.length-1;
        faceIndices[makeFaceTag([vtxList[1], ctrId])] = if1_id;

        faces ~= new USGFace([vtxList[2], ctrId]);
        auto if2_id = faces.length-1;
        faceIndices[makeFaceTag([vtxList[2], ctrId])] = if2_id;

        nfaces = faces.length;
        // We'll replace the present cell with sub-triangle 0,
        // then add sub-triangle 1 and 2 as new cells.
        cells[id].vtx_id_list[0] = ctrId;
        cells[id].face_id_list[0] = if1_id;
        cells[id].face_id_list[2] = if2_id;
        cells[id].outsign_list[0] = -1;

        auto nc1_id = cells.length;
        cells ~= new USGCell(USGCell_type.triangle, nc1_id,
                             [vtxList[0], ctrId, vtxList[2]],
                             [if0_id, if2_id, cells[id].face_id_list[2]],
                             [1, -1, 1]);

        auto nc2_id = cells.length;
        cells ~= new USGCell(USGCell_type.triangle, nc2_id,
                             [vtxList[0], vtxList[1], ctrId],
                             [cells[id].face_id_list[0], if1_id, if0_id],
                             [1, 1, -1]);


        ncells = cells.length;

        // Set left and right cells on new faces.
        faces[if0_id].left_cell = cells[nc1_id];
        faces[if0_id].right_cell = cells[nc2_id];

        faces[if1_id].left_cell = cells[nc2_id];
        faces[if1_id].right_cell = cells[id];

        faces[if2_id].left_cell = cells[id];
        faces[if2_id].right_cell = cells[nc1_id];
    }



    // ------------------------
    // Import-from-file methods.
    // ------------------------

    void link_faces_back_to_cells()
    {
        // Make links from faces back to cells.
        foreach (c; cells) {
            size_t nf = c.face_id_list.length;
            assert (nf == c.outsign_list.length, "Mismatch in number of faces and number of outsigns");
            foreach (j; 0 .. nf) {
                auto f = faces[c.face_id_list[j]];
                if (c.outsign_list[j] == 1) {
                    f.left_cell = c;
                } else {
                    f.right_cell = c;
                }
            }
        }
    } // end link_faces_back_to_cells()

    override void read_from_gzip_file(string fileName, double scale=1.0)
    // This function, together with the constructors (from strings) for
    // the classes USGFace, USGCell and BoundaryFaceSet (above),
    // define the Eilmer4 native format.
    // For output, there is the matching write_to_gzip_file() below.
    // scale = unit length in metres
    {
        auto byLine = new GzipByLine(fileName);
        auto line = byLine.front; byLine.popFront();
        string format_version;
        formattedRead(line, "unstructured_grid %s", &format_version);
        if (format_version != "1.0") {
            throw new Error("UnstructuredGrid.read_from_gzip_file(): " ~
                            "format version found: " ~ format_version);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "label: %s", &label);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "dimensions: %d", &dimensions);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "vertices: %d", &nvertices);
        double x, y, z;
        vertices.length = 0;
        foreach (i; 0 .. nvertices) {
            line = byLine.front; byLine.popFront();
            // Note that the line starts with whitespace.
            formattedRead(line, " %g %g %g", &x, &y, &z);
            vertices ~= Vector3(scale*x, scale*y, scale*z);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "faces: %d", &nfaces);
        faces.length = 0;
        foreach (i; 0 .. nfaces) {
            line = byLine.front; byLine.popFront();
            auto myFace = new USGFace(line);
            faceIndices[makeFaceTag(myFace.vtx_id_list)] = faces.length;
            faces ~= myFace;
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "cells: %d", &ncells);
        cells.length = 0;
        foreach (i; 0 .. ncells) {
            line = byLine.front; byLine.popFront();
            cells ~= new USGCell(line, i);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "boundaries: %d", &nboundaries);
        boundaries.length = 0;
        foreach (i; 0 .. nboundaries) {
            line = byLine.front; byLine.popFront();
            boundaries ~= new BoundaryFaceSet(line);
        }
        link_faces_back_to_cells();
    } // end read_from_gzip_file()

    override void read_from_raw_binary_file(string fileName, double scale=1.0)
    // This function, together with the constructors (from a binary file) for
    // the classes USGFace, USGCell and BoundaryFaceSet (above),
    // define the Eilmer4 native format.
    // For output, there is the matching write_to_raw_binary_file() below.
    // scale = unit length in metres
    {
        File fin = File(fileName, "rb");
        string expected_header = "unstructured_grid 1.0";
        char[] found_header = new char[expected_header.length];
        fin.rawRead(found_header);
        if (found_header != expected_header) {
            throw new Error("UnstructuredGrid.read_from_raw_binary_file(): " ~
                            "unexpected header: " ~ to!string(found_header));
        }
        int[1] buf1; fin.rawRead(buf1);
        int label_length = buf1[0];
        if (label_length > 0) {
            char[] found_label = new char[label_length];
            fin.rawRead(found_label);
            label = to!string(found_label);
        }
        int[5] buf5; fin.rawRead(buf5);
        dimensions = buf5[0];
        nvertices = buf5[1];
        nfaces = buf5[2];
        ncells = buf5[3];
        nboundaries = buf5[4];
        double[3] xyz;
        double x, y, z;
        vertices.length = 0;
        foreach (i; 0 .. nvertices) {
            fin.rawRead(xyz);
            vertices ~= Vector3(scale*xyz[0], scale*xyz[1], scale*xyz[2]);
        }
        faces.length = 0;
        foreach (i; 0 .. nfaces) {
            auto myFace = new USGFace(fin);
            faceIndices[makeFaceTag(myFace.vtx_id_list)] = faces.length;
            faces ~= myFace;
        }
        cells.length = 0;
        foreach (i; 0 .. ncells) {
            cells ~= new USGCell(fin, i);
        }
        boundaries.length = 0;
        foreach (i; 0 .. nboundaries) {
            boundaries ~= new BoundaryFaceSet(fin);
        }
        link_faces_back_to_cells();
    } // end read_from_raw_binary_file()

    void read_from_su2_file(string fileName,
                            bool true_centroids,
                            double scale=1.0,
                            bool expect_gmsh_order_for_wedges=false)
    // Information on the su2 file format from
    // https://su2code.github.io/docs_v7/Mesh-File/
    // scale = unit length in metres
    {
        writef("Importing: %s ", fileName);
        auto f = File(fileName, "r");
        string getHeaderContent(string target)
        // Helper function to proceed through file, line-by-line,
        // looking for a particular header line.
        // Returns the content from the header line and leaves the file
        // at the next line to be read, presumably with expected data.
        {
            while (!f.eof) {
                auto line = f.readln().strip();
                if (canFind(line, target)) {
                    auto tokens = line.split("=");
                    return tokens[1].strip();
                }
            } // end while
            return ""; // didn't find the target
        }
        dimensions = to!int(getHeaderContent("NDIME"));
        //
        // Since we are not sure of the order of the blocks (elements or points),
        // let's be safe and look from the beginning.
        f.rewind();
        ncells = to!size_t(getHeaderContent("NELEM"));
        cells.length = ncells;
        foreach(i; 0 .. ncells) {
            // For very large files, emit a character periodically,
            // to show that work is being done.
            if ((i % 100_000) == 0) { stdout.write("."); stdout.flush();}
            auto lineContent = f.readln().strip();
            // We are expecting the following items on the line:
            // a. vtk-element-type
            // b. vertex indices for the element
            // c. cell index
            // If the cell index is not available, use the current value of i.
            auto tokens = lineContent.split();
            int vtk_element_type = to!int(tokens[0]);
            size_t indx = i; // default index if not provided on line
            int nv_expected = n_vertices_for_vtk_element(vtk_element_type);
            size_t[] vtx_id_list;
            if (tokens.length == 2+nv_expected) {
                // All items are present, including the cell index at the end.
                indx = to!size_t(tokens[$-1]);
                foreach(j; 1 .. tokens.length-1) { vtx_id_list ~= to!size_t(tokens[j]); }
            } else if (tokens.length == 1+nv_expected) {
                // Presume that the cell index is not available at the end of the line.
                foreach(j; 1 .. tokens.length) { vtx_id_list ~= to!size_t(tokens[j]); }
            } else {
                string msg = "Incorrect number of items on line for element.";
                msg ~= format("\nCell i: %d line:\"%s\"", i, lineContent);
                throw new GeometryException(msg);
            }
            USGCell_type cell_type = convert_cell_type(vtk_element_type);
            if (cell_type == USGCell_type.none) {
                string msg = "Unknown element type.";
                msg ~= format("\nCell i: %d line:\"%s\"", i, lineContent);
                throw new GeometryException(msg);
            }
            if (cell_type == USGCell_type.wedge && expect_gmsh_order_for_wedges) {
                // We assume that we have picked up the vertex indices for a wedge
                // in GMSH order but we need to store them in VTK order.
                vtx_id_list = [vtx_id_list[1], vtx_id_list[0], vtx_id_list[2],
                               vtx_id_list[4], vtx_id_list[3], vtx_id_list[5]];
            }
            size_t[] face_id_list; // empty, so far
            int[] outsign_list; // empty, so far
            // Once we have created the full list of vertices, we will be able to make
            // the set of faces and then come back to filling in the face_id_list and
            // outsign_list for each cell.
            cells[indx] = new USGCell(cell_type, indx, vtx_id_list, face_id_list, outsign_list);
        } // end foreach i .. ncells
        //
        size_t uninitialized_cell_count = 0;
        foreach(i; 0 .. cells.length) {
            if (!cells[i]) { ++uninitialized_cell_count; }
        }
        if (uninitialized_cell_count > 0) {
            writefln("Warning: There are %d initialized cells and %d uninitialized cells.",
                     ncells, uninitialized_cell_count);
            // [TODO] If we have any uninitialized cells, we should
            // compress the array to eliminate empty elements, however,
            // the following code is tolerant of uninitialized cells.
        }
        //
        // Now, look for the vertices.
        // It may be that the block of vertices is located in the SU2 file
        // before the block that defines the elements, so let's look from
        // the start of the file.
        f.rewind();
        nvertices = to!size_t(getHeaderContent("NPOIN"));
        vertices.length = nvertices;
        foreach(i; 0 .. nvertices) {
            // For very large files, emit a character periodically,
            // to show that work is being done.
            if ((i % 100_000) == 0) { stdout.write("."); stdout.flush();}
            auto tokens = f.readln().strip().split();
            auto ntokens = tokens.length;
            double x=0.0; double y=0.0; double z = 0.0; size_t indx = 0;
            // We expect to see a coordinate value for each dimension plus
            // an integer index for the vertex.
            // The index may or may notbe present.
            if (dimensions == 2) {
                x = scale * to!double(tokens[0]);
                y = scale * to!double(tokens[1]);
                if (ntokens == 3) {
                    indx = to!size_t(tokens[2]);
                } else {
                    indx = i;
                }
            } else {
                assert(dimensions == 3, "invalid dimensions");
                x = scale * to!double(tokens[0]);
                y = scale * to!double(tokens[1]);
                z = scale * to!double(tokens[2]);
                if (ntokens == 4) {
                    indx = to!size_t(tokens[3]);
                } else {
                    indx = i;
                }
            }
            vertices[indx].set(x, y, z);
        } // end foreach i .. nvertices
        //
        if (dimensions == 2) {
            // In 2D, the flow solver code assumes that the cycle of vertices
            // defining the cell is counterclockwise, looking down the negative z-axis,
            // toward the x,y-plane.
            // If any of the cells that have been read from the file presently have
            // their normal vector pointing in the negative z-direction,
            // we should reverse the order of their vertex list.
            foreach(c; cells) {
                switch (c.cell_type) {
                case USGCell_type.triangle:
                    Vector3 a = Vector3(vertices[c.vtx_id_list[1]]);
                    a -= vertices[c.vtx_id_list[0]];
                    Vector3 b = Vector3(vertices[c.vtx_id_list[2]]);
                    b -= vertices[c.vtx_id_list[0]];
                    Vector3 n;
                    cross(n,a,b);
                    if (n.z < 0.0) { reverse(c.vtx_id_list); }
                    break;
                case USGCell_type.quad:
                    Vector3 a = Vector3(vertices[c.vtx_id_list[1]]);
                    a -= vertices[c.vtx_id_list[0]];
                    a += vertices[c.vtx_id_list[2]];
                    a -= vertices[c.vtx_id_list[3]];
                    Vector3 b = Vector3(vertices[c.vtx_id_list[3]]);
                    b -= vertices[c.vtx_id_list[0]];
                    b += vertices[c.vtx_id_list[2]];
                    b -= vertices[c.vtx_id_list[1]];
                    Vector3 n;
                    cross(n,a,b);
                    if (n.z < 0.0) { reverse(c.vtx_id_list); }
                    break;
                default:
                    throw new GeometryException("invalid cell type for 2D grid");
                } // end switch
            } // end foreach c
        } // end if (dimensions == 2)
        //
        // Now that we have the full list of cells and vertices assigned to each cell,
        // we can construct the faces between cells and along the boundaries.
        //
        void add_face_to_cell(ref USGCell cell, size_t[] corners, ref Vector3 cell_centroid)
        {
            // If the face is new, we add it to the list, else we use the face
            // already stored within the list, so that it may no longer be
            // outward pointing for this cell.
            //
            size_t[] my_vtx_id_list;
            foreach(i; corners) { my_vtx_id_list ~= cell.vtx_id_list[i]; }
            size_t face_indx = 0;
            string faceTag = makeFaceTag(my_vtx_id_list);
            if (faceTag in faceIndices) {
                // Use the face already present.
                face_indx = faceIndices[faceTag];
            } else {
                // Since we didn't find the face already, construct it.
                face_indx = faces.length;
                faces ~= new USGFace(my_vtx_id_list);
                faceIndices[faceTag] = face_indx;
            }
            cell.face_id_list ~= face_indx;
            //
            // Note that, from this point, we work with the face from the collection
            // which, if it was an existing face, will not have the same orientation
            // as indicated by the cycle of vertices passed in to this function call.
            //
            // Pointers to all vertices on face.
            Vector3*[] v;
            foreach (i; 0 .. corners.length) {
                v ~= &vertices[faces[face_indx].vtx_id_list[i]];
            }
            //
            // Determine if cell is on left- or right-side of face.
            //
            bool onLeft;
            switch (corners.length) {
            case 2:
                // Two-dimensional simulation:
                // The face is a line in the xy-plane.
                onLeft = on_left_of_xy_line(*(v[0]), *(v[1]), cell_centroid);
                break;
            case 3:
                // Three-dimensional simulation:
                // The face has been defined with the vertices being in
                // a counter-clockwise cycle when viewed from outside of the cell.
                // We expect negative volumes for the pyramid having the cell's
                // midpoint as its apex.
                onLeft = (tetrahedron_volume(*(v[0]), *(v[1]), *(v[2]), cell_centroid) < 0.0);
                break;
            case 4:
                // Three-dimensional simulation:
                // The face has been defined with the vertices being in
                // a counter-clockwise cycle when viewed from outside of the cell.
                // We expect negative volumes for the pyramid having the cell's
                // midpoint as its apex.
                Vector3 vmid = 0.25*( *(v[0]) + *(v[1]) + *(v[2]) + *(v[3]) );
                onLeft = (tetragonal_dipyramid_volume(*(v[0]), *(v[1]), *(v[2]), *(v[3]),
                                                      vmid, cell_centroid) < 0.0);
                break;
            default:
                throw new GeometryException("invalid number of corners on face.");
            }
            if (onLeft) {
                if (faces[face_indx].left_cell) {
                    throw new GeometryException("face already has a left cell");
                }
                faces[face_indx].left_cell = cell;
                cell.outsign_list ~=  1;
            } else {
                if (faces[face_indx].right_cell) {
                    throw new GeometryException("face already has a right cell");
                }
                faces[face_indx].right_cell = cell;
                cell.outsign_list ~= -1;
            }
            return;
        } // end add_face_to_cell()
        //
        foreach(i, cell; cells) {
            // For a very large grid, emit a character periodically,
            // to show that work is being done.
            if ((i % 100_000) == 0) { stdout.write("."); stdout.flush();}
            if (!cell) continue;
            // Attach the faces to each cell. In 2D, faces are defined as lines.
            // As we progress along the line the face normal is pointing to the right.
            // In 3D, a counter-clockwise cycles of points plus the right-hand rule
            // define the face normal. Whether a face points out of or into a cell
            // will be determined and remembered when we add the face to the cell.
            Vector3 centroid;
            number volume;
            compute_cell_properties(i, true_centroids, centroid, volume);
            if (dimensions == 2) {
                switch(cell.cell_type) {
                case USGCell_type.triangle:
                    add_face_to_cell(cell, [0,1], centroid);
                    add_face_to_cell(cell, [1,2], centroid);
                    add_face_to_cell(cell, [2,0], centroid);
                    break;
                case USGCell_type.quad:
                    add_face_to_cell(cell, [2,3], centroid); // north
                    add_face_to_cell(cell, [1,2], centroid); // east
                    add_face_to_cell(cell, [0,1], centroid); // south
                    add_face_to_cell(cell, [3,0], centroid); // west
                    break;
                default:
                    throw new GeometryException("invalid cell type in 2D");
                }
            } else {
                assert(dimensions == 3, "invalid dimensions");
                switch(cell.cell_type) {
                case USGCell_type.tetra:
                    add_face_to_cell(cell, [0,1,2], centroid);
                    add_face_to_cell(cell, [0,1,3], centroid);
                    add_face_to_cell(cell, [1,2,3], centroid);
                    add_face_to_cell(cell, [2,0,3], centroid);
                    break;
                case USGCell_type.hexahedron:
                    add_face_to_cell(cell, [2,3,7,6], centroid); // north
                    add_face_to_cell(cell, [1,2,6,5], centroid); // east
                    add_face_to_cell(cell, [1,0,4,5], centroid); // south
                    add_face_to_cell(cell, [0,3,7,4], centroid); // west
                    add_face_to_cell(cell, [4,5,6,7], centroid); // top
                    add_face_to_cell(cell, [0,1,2,3], centroid); // bottom
                    break;
                case USGCell_type.wedge:
                    add_face_to_cell(cell, [0,1,2], centroid);
                    add_face_to_cell(cell, [3,4,5], centroid);
                    add_face_to_cell(cell, [0,2,5,3], centroid);
                    add_face_to_cell(cell, [0,3,4,1], centroid);
                    add_face_to_cell(cell, [1,4,5,2], centroid);
                    break;
                case USGCell_type.pyramid:
                    add_face_to_cell(cell, [0,1,2,3], centroid);
                    add_face_to_cell(cell, [0,1,4], centroid);
                    add_face_to_cell(cell, [1,2,4], centroid);
                    add_face_to_cell(cell, [2,3,4], centroid);
                    add_face_to_cell(cell, [3,0,4], centroid);
                    break;
                default:
                    throw new GeometryException("invalid cell type in 3D");
                }
            }
            Vector3 varea = vectorAreaOfCell(cell);
            double relTol=1.0e-9; double absTol=1.0e-9;
            if (!approxEqualVectors(varea, Vector3(0,0,0), relTol, absTol)) {
                string errMsg = format("cell id=%d, non-zero vector area=%s\n   cell is %s",
                                       i, to!string(varea), cell.toIOString());
                foreach (j, fid; cell.face_id_list) {
                    errMsg ~= format("\n    j=%d outsign=%d vectorArea=%s", j, cell.outsign_list[j],
                                     vectorAreaOfFacet(faces[fid].vtx_id_list));
                }
                throw new GeometryException(errMsg);
            }
        } // end foreach cell
        nfaces = faces.length;
        //
        // Now that we have a full set of cells and faces,
        // make lists of the boundary faces.
        //
        foreach (myface; faces) { myface.is_on_boundary = false; } // will remark boundary faces below
        // Note that we do not rewind the file because
        // we expect the boundary sets last, always.
        nboundaries = to!size_t(getHeaderContent("NMARK"));
        foreach(i; 0 .. nboundaries) {
            string tag = getHeaderContent("MARKER_TAG");
            size_t[] face_id_list;
            int[] outsign_list;
            size_t nelem = to!size_t(getHeaderContent("MARKER_ELEMS"));
            foreach(j; 0 .. nelem) {
                auto tokens = f.readln().strip().split();
                int vtk_type = to!int(tokens[0]);
                size_t[] my_vtx_id_list;
                foreach(k; 1 .. tokens.length) { my_vtx_id_list ~= to!size_t(tokens[k]); }
                string faceTag = makeFaceTag(my_vtx_id_list);
                if (faceTag in faceIndices) {
                    size_t face_indx = faceIndices[faceTag];
                    USGFace my_face = faces[face_indx];
                    assert(my_face.left_cell || my_face.right_cell, "face is not properly connected");
                    if (my_face.left_cell && !(my_face.right_cell)) {
                        outsign_list ~= 1;
                    } else if ((!my_face.left_cell) && my_face.right_cell) {
                        outsign_list ~= -1;
                    } else {
                        throw new GeometryException("appears to be an interior face");
                    }
                    face_id_list ~= face_indx;
                    my_face.is_on_boundary = true;
                } else {
                    throw new GeometryException("cannot find face in collection");
                }
            } // end foreach j
            boundaries ~= new BoundaryFaceSet(tag, face_id_list, outsign_list);
        } // end foreach i
        //
        // At this point, all interior faces should have a cell attached to both
        // left- and right-sides while boundary faces should have a cell attached
        // to one-side only.
        // If there are any faces that have only one cell attached and are not in
        // a boundary-set already, assign them to the orphan-faces boundary set.
        size_t[] orphan_face_id_list;
        int[] orphan_face_outsign_list;
        foreach (i, myface; faces) {
            if (myface.is_on_boundary) continue;
            if (myface.left_cell && myface.right_cell) continue; // interior face with two cells
            orphan_face_id_list ~= i;
            if (myface.left_cell && !(myface.right_cell)) {
                orphan_face_outsign_list ~= 1;
            } else if (!(myface.left_cell) && myface.right_cell) {
                orphan_face_outsign_list ~= -1;
            }
        } // end foreach
        if (orphan_face_id_list.length > 0) {
            boundaries ~= new BoundaryFaceSet("ORPHAN_FACES", orphan_face_id_list,
                                              orphan_face_outsign_list);
            nboundaries = boundaries.length;
        }
        //
        Vector3 varea = Vector3(0,0,0);
        foreach (b; boundaries) {
            foreach (j, fid; b.face_id_list) {
                varea += b.outsign_list[j] * vectorAreaOfFacet(faces[fid].vtx_id_list);
            }
        }
        double relTol=1.0e-9; double absTol=1.0e-9;
        if (!approxEqualVectors(varea, Vector3(0,0,0), relTol, absTol)) {
            string errMsg = format("SU2 grid has non-zero bounding vector area=", varea);
            throw new GeometryException(errMsg);
        }
        //
        // If we arrive here, the import of the SU2 grid seems to have been successful
        writefln(" (done) scale=%s dim=%s ncells=%d nbfs= %d", scale, dimensions, ncells, boundaries.length);
        return;
    } // end read_from_su2_text_file()

    void read_from_vtk_text_file(string fileName, double scale=1.0)
    // scale = unit length in metres
    {
        string[] tokens;
        auto f = File(fileName, "r");
        read_VTK_header_line("vtk", f);
        label = f.readln().strip();
        read_VTK_header_line("ASCII", f);
        read_VTK_header_line("UNSTRUCTURED_GRID", f);
        tokens = f.readln().strip().split();
        throw new Error("read_from_vtk_text_file() is not finished, yet.");
        // For VTK files, we need to work out how to extract
        // the topological details for the incoming grid.
    } // end read_from_vtk_text_file()

    override void write_to_gzip_file(string fileName)
    // This function, together with the toIOstring methods for
    // the classes USGFace, USGCell and BoundaryFaceSet (way above),
    // defines the Eilmer4 native format.
    // For input, there is the matching read_from_gzip_file(), above.
    {
        auto fout = new GzipOut(fileName);
        fout.compress("unstructured_grid 1.0\n");
        fout.compress(format("label: %s\n", label));
        fout.compress(format("dimensions: %d\n", dimensions));
        fout.compress(format("vertices: %d\n", nvertices));
        foreach (v; vertices) {
            fout.compress(format("%.18e %.18e %.18e\n", v.x.re, v.y.re, v.z.re));
        }
        fout.compress(format("faces: %d\n", nfaces));
        foreach (f; faces) { fout.compress(f.toIOString ~ "\n"); }
        fout.compress(format("cells: %d\n", ncells));
        foreach (c; cells) { fout.compress(c.toIOString ~ "\n"); }
        fout.compress(format("boundaries: %d\n", boundaries.length));
        foreach (b; boundaries) { fout.compress(b.toIOString ~ "\n"); }
        fout.finish();
    } // end write_to_gzip_file()

    override void write_to_raw_binary_file(string fileName)
    // This function, together with the toIOstring methods for
    // the classes USGFace, USGCell and BoundaryFaceSet (way above),
    // defines the Eilmer4 native binary format.
    // For input, there is the matching read_from_raw_binary_file(), above.
    {
        File fout = File(fileName, "wb");
        fout.rawWrite(to!(char[])("unstructured_grid 1.0"));
        int[1] buf1; buf1[0] = to!int(label.length); fout.rawWrite(buf1);
        if (label.length > 0) { fout.rawWrite(to!(char[])(label)); }
        int[5] buf5; buf5[0] = to!int(dimensions);
        buf5[1] = to!int(nvertices);
        buf5[2] = to!int(nfaces);
        buf5[3] = to!int(ncells);
        buf5[4] = to!int(boundaries.length);
        fout.rawWrite(buf5);
        double[3] xyz;
        foreach (v; vertices) {
            xyz[0] = v.x.re; xyz[1] = v.y.re; xyz[2] = v.z.re;
            fout.rawWrite(xyz);
        }
        foreach (f; faces) { f.rawWrite(fout); }
        foreach (c; cells) { c.rawWrite(fout); }
        foreach (b; boundaries) { b.rawWrite(fout); }
        fout.close();
    } // end write_to_raw_binary_file()

    override void write_to_vtk_file(string fileName)
    {
        auto f = File(fileName, "w");
        f.writeln("# vtk DataFile Version 2.0");
        f.writeln(label);
        f.writeln("ASCII");
        f.writeln("");
        f.writeln("DATASET UNSTRUCTURED_GRID");
        f.writefln("POINTS %d float", nvertices);
        foreach (v; vertices) { f.writefln("%.18e %.18e %.18e", v.x.re, v.y.re, v.z.re); }
        f.writeln("");
        int n_ints = 0; // number of integers to describe connections
        foreach (c; cells) { n_ints += 1 + c.vtx_id_list.length; }
        f.writefln("CELLS %d %d", ncells, n_ints);
        foreach (c; cells) {
            f.write(c.vtx_id_list.length);
            foreach (vtx_id; c.vtx_id_list) { f.write(" ", vtx_id); }
            f.writeln();
        }
        f.writefln("CELL_TYPES %d", ncells);
        foreach (c; cells) {
            f.writeln(vtk_element_types[c.cell_type]);
        }
        f.close();
    } // end write_to_vtk_file()

    override void write_to_su2_file(string fileName, double scale=1.0,
                                    bool use_gmsh_order_for_wedges=true)
    {
        auto f = File(fileName, "w");
        f.writefln("NDIME= %d", dimensions);
        f.writeln("");
        f.writefln("NELEM= %d", ncells);
        foreach (i, c; cells) {
            f.write(vtk_element_types[c.cell_type]);
            if (c.cell_type == USGCell_type.wedge && use_gmsh_order_for_wedges) {
                size_t[] vtx_id_list = c.vtx_id_list.dup();
                // We have the vertex indices for a wedge stored in VTK order
                // but we are requested to output them in GMSH order.
                vtx_id_list = [vtx_id_list[1], vtx_id_list[0], vtx_id_list[2],
                               vtx_id_list[4], vtx_id_list[3], vtx_id_list[5]];
                foreach (vtx_id; vtx_id_list) { f.write(" ", vtx_id); }
            } else {
                foreach (vtx_id; c.vtx_id_list) { f.write(" ", vtx_id); }
            }
            f.writefln(" %d", i); // cell index
        }
        f.writeln("");
        f.writefln("NPOIN= %d", nvertices);
        foreach (i, v; vertices) {
            if (dimensions == 3) {
                f.writefln("%.18e %.18e %.18e %d", v.x.re/scale, v.y.re/scale, v.z.re/scale, i);
            } else {
                f.writefln("%.18e %.18e %d", v.x.re/scale, v.y.re/scale, i);
            }
        }
        f.writeln("");
        f.writefln("NMARK= %d", nboundaries);
        foreach (b; boundaries) {
            f.writefln("MARKER_TAG= %s", b.tag);
            f.writefln("MARKER_ELEMS= %d", b.face_id_list.length);
            foreach (faceIndx; b.face_id_list) {
                auto vtx_id_list = faces[faceIndx].vtx_id_list;
                if (dimensions == 3) {
                    switch (vtx_id_list.length) {
                    case 3: f.writef("%d", VTKElement.triangle); break;
                    case 4: f.writef("%d", VTKElement.quad); break;
                    default: f.write("oops, should not have this number of vertices");
                    }
                } else {
                    // In 2D, we have only lines as the boundary elements.
                    f.writef("%d", VTKElement.line);
                }
                foreach (vtx_id; vtx_id_list) { f.writef(" %d", vtx_id); }
                f.writeln();
            }
        }
        f.close();
    } // end write_to_su2_file()

    void write_openFoam_polyMesh(string topLevelDir)
    // The OpenFoam mesh format is defined primarily in terms of the face polygons.
    // There is a required order for the faces in which all internal faces appear first,
    // followed by the boundary sets of faces.  This order is required so that each
    // boundary set can be simply identified as a starting face and a number of faces.
    // There is also a constraint that each face's unit normal points out of its
    // "owning" cell.  For an internal cell, between two cells, the owning cell is
    // the one with the smaller index.
    //
    // Files that will be written:
    // <topLevelDir>/polyMesh/points
    //                       /faces
    //                       /owner
    //                       /neighbour
    //                       /boundary
    {
        // First, we will construct a description of the unstructured mesh
        // in terms of quantities that OpenFOAM uses.
        //
        // Construct dictionary of cell indices so that we can look up
        // the indices of the left and right cells.
        //
        int[USGCell] cellId; foreach (i, c; cells) { cellId[c] = to!int(i); }
        //
        // Set up the list of internal-to-block face ids and
        // lists of owner and neighbour cells.
        //
        size_t[] internal_face_id_list; // limited to internal faces only
        int[] owner_cell_id_list; // an entry for all faces
        int[] neighbour_cell_id_list; // an entry for all faces
        foreach (i, face; faces) {
            int owner_id = -1;
            int neighbour_id = -1;
            if (face.left_cell is null) {
                owner_id = cellId[face.right_cell];
            } else if (face.right_cell is null) {
                owner_id = cellId[face.left_cell];
            } else {
                internal_face_id_list ~= i;
                if (cellId[face.left_cell] < cellId[face.right_cell]) {
                    owner_id = cellId[face.left_cell];
                    neighbour_id = cellId[face.right_cell];
                } else {
                    owner_id = cellId[face.right_cell];
                    neighbour_id = cellId[face.left_cell];
                }
            }
            assert(owner_id >= 0, "face seems to be not owned by a cell"); // [TODO] more info
            owner_cell_id_list ~= owner_id;
            neighbour_cell_id_list ~= neighbour_id;
        }
        size_t face_tally = internal_face_id_list.length;
        foreach (b; boundaries) { face_tally += b.face_id_list.length; }
        assert(face_tally == faces.length, "mismatch in number of faces in boundary sets");
        //
        // We need to sift through the list of cells and check the orientations
        // of all of the attached faces.
        // OpenFoam expects each face to be pointing out of the owner cell
        // so, where we currently have the outsign negative for the owner cell,
        // we need to reverse the order when writing the vertex ids for the face.
        //
        bool[] reverse_vtx_order_for_face;
        reverse_vtx_order_for_face.length = faces.length;
        bool[] face_already_checked;
        foreach (fid; 0 .. faces.length) { face_already_checked ~= false; }
        foreach (i, c; cells) {
            foreach (j, fid; c.face_id_list) {
                USGFace f = faces[fid];
                if (face_already_checked[fid]) continue;
                if (f.left_cell is null) {
                    // block-boundary cell
                    assert(c.outsign_list[j] == -1, "expected inward-pointing face");
                    reverse_vtx_order_for_face[fid] = true;
                    face_already_checked[fid] = true;
                    continue;
                }
                if (f.right_cell is null) {
                    // block-boundary cell
                    assert(c.outsign_list[j] == 1, "expected outward-pointing face");
                    reverse_vtx_order_for_face[fid] = false;
                    face_already_checked[fid] = true;
                    continue;
                }
                // From this point, assume an interior face with two attached cells.
                size_t other_cell_id;
                if (c.outsign_list[j] == 1) {
                    other_cell_id = cellId[f.right_cell];
                } else {
                    other_cell_id = cellId[f.left_cell];
                }
                if (((i < other_cell_id) && (c.outsign_list[j] == 1)) ||
                    ((i > other_cell_id) && (c.outsign_list[j] == -1))) {
                    // Current face orientaion is consistent with openFoam arrangement.
                    reverse_vtx_order_for_face[fid] = false;
                } else {
                    // Current face orientation is opposite to openFoam arrangement.
                    reverse_vtx_order_for_face[fid] = true;
                } // end if
                face_already_checked[fid] = true;
            } // end foreach fid
        } // end foreach c
        //
        // Now, we are ready to write the files.
        //
        string polyMeshDir = topLevelDir ~ "/polyMesh";
        if (!exists(polyMeshDir)) { mkdirRecurse(polyMeshDir); }
        //
        auto f = File(polyMeshDir~"/points", "w");
        f.writeln("FoamFile\n{");
        f.writeln(" version   2.0;");
        f.writeln(" format    ascii;");
        f.writeln(" class     vectorField;");
        f.writeln(" location  \"constant/polyMesh\";");
        f.writeln(" object    points;");
        f.writeln("}");
        f.writefln("%d\n(", vertices.length);
        foreach (i,v; vertices) { f.writefln(" (%.18e %.18e %.18e)", v.x.re, v.y.re, v.z.re); }
        f.writeln(")");
        f.close();
        //
        f = File(polyMeshDir~"/faces", "w");
        f.writeln("FoamFile\n{");
        f.writeln(" version   2.0;");
        f.writeln(" format    ascii;");
        f.writeln(" class     faceList;");
        f.writeln(" location  \"constant/polyMesh\";");
        f.writeln(" object    faces;");
        f.writeln("}");
        string index_str(const ref USGFace f, bool reverse_vtx_order)
        {
            string str = " " ~ to!string(f.vtx_id_list.length) ~ "(";
            auto my_vtx_id_list = f.vtx_id_list.dup();
            if (reverse_vtx_order) { reverse(my_vtx_id_list); }
            foreach (j, vtx_id; my_vtx_id_list) {
                if (j > 0) { str ~= " "; }
                str ~= to!string(vtx_id);
            }
            str ~= ")";
            return str;
        }
        f.writefln("%d\n(", faces.length);
        foreach (i; internal_face_id_list) {
            f.writeln(index_str(faces[i], reverse_vtx_order_for_face[i]));
        }
        foreach (b; boundaries) {
            foreach (i; b.face_id_list) {
                f.writeln(index_str(faces[i], reverse_vtx_order_for_face[i]));
            }
        }
        f.writeln(")");
        f.close();
        //
        f = File(polyMeshDir~"/owner", "w");
        f.writeln("FoamFile\n{");
        f.writeln(" version   2.0;");
        f.writeln(" format    ascii;");
        f.writeln(" class     labelList;");
        f.writeln(" location  \"constant/polyMesh\";");
        f.writeln(" object    owner;");
        f.writeln("}");
        f.writefln("%d\n(", faces.length);
        foreach (i; internal_face_id_list) { f.writefln(" %d", owner_cell_id_list[i]); }
        foreach (b; boundaries) {
            foreach (i; b.face_id_list) { f.writefln(" %d", owner_cell_id_list[i]); }
        }
        f.writeln(")");
        f.close();
        //
        f = File(polyMeshDir~"/neighbour", "w");
        f.writeln("FoamFile\n{");
        f.writeln(" version   2.0;");
        f.writeln(" format    ascii;");
        f.writeln(" class     labelList;");
        f.writeln(" location  \"constant/polyMesh\";");
        f.writeln(" object    neighbour;");
        f.writeln("}");
        f.writefln("%d\n(", internal_face_id_list.length);
        foreach (i; internal_face_id_list) { f.writefln(" %d", neighbour_cell_id_list[i]); }
        f.writeln(")");
        f.close();
        //
        f = File(polyMeshDir~"/boundary", "w");
        f.writeln("FoamFile\n{");
        f.writeln(" version   2.0;");
        f.writeln(" format    ascii;");
        f.writeln(" class     polyBoundaryMesh;");
        f.writeln(" location  \"constant/polyMesh\";");
        f.writeln(" object    boundary;");
        f.writeln("}");
        size_t bndry_count = 0; // for non-empty boundary sets
        foreach (b; boundaries) { if (b.face_id_list.length > 0) { ++bndry_count; } }
        f.writefln("%d\n(", bndry_count);
        size_t startFace = internal_face_id_list.length;
        foreach (i, b; boundaries) {
            if (b.face_id_list.length == 0) continue; // skip past empty sets
            string label = b.tag; // Nominally, we would like to use the assigned tag
            // but, if it is not useful, make up something that is.
            if (label.length == 0) { label = format("boundary%04d", i); }
            f.writefln(" %s\n {", label);
            f.writefln("  type      %s;", "patch");
            f.writefln("  nFaces    %d;", b.face_id_list.length);
            f.writefln("  startFace %d;", startFace);
            f.writeln(" }");
            startFace += b.face_id_list.length; // for the next boundary
        }
        f.writeln(")");
        f.close();
        //
        // Finally, write the face sets using openFoam index order for the faces.
        // So far these faceSets  are named for the objects in the boundary file.
        string setsDir = topLevelDir ~ "/polyMesh/sets";
        if (!exists(setsDir)) { mkdirRecurse(setsDir); }
        startFace = internal_face_id_list.length;
        foreach (i, b; boundaries) {
            size_t nf = b.face_id_list.length;
            if (nf == 0) continue; // skip past empty sets
            string objName = b.tag; // Nominally, we would like to use the assigned tag
            // but, if it is not useful, make up something that is.
            if (objName.length == 0) { label = format("boundary%04d", i); }
            f = File(setsDir~"/"~objName, "w");
            f.writeln("FoamFile\n{");
            f.writeln(" version   2.0;");
            f.writeln(" format    ascii;");
            f.writeln(" class     faceSet;");
            f.writeln(" location  \"constant/polyMesh/sets\";");
            f.writefln(" object    %s;", objName);
            f.writeln("}");
            f.writefln("%d\n(", nf);
            foreach (j; 0 .. nf) { f.writefln(" %d", j+startFace); }
            f.writeln(")");
            startFace += nf; // for the next boundary
        } // end foreach boundary
    } // end write_openFoam_polyMesh()

    UnstructuredGrid joinGrid(const UnstructuredGrid other,
                              double relTol=1.0e-6, double absTol=1.0e-9,
                              bool checkBoundariesOnly=true,
                              bool openFoam=false, int dimensions=3)
    {
        // We consider *this as the "master" grid and will join
        // a copy of the unique bits of the other grid to it.
        //
        if (this == other) return this; // nothing to do
        // Some checking of the incoming data.
        assert(nvertices == vertices.length, "wrong number of vertices in master grid");
        assert(other.nvertices == other.vertices.length, "wrong number of vertices in other grid");
        assert(nfaces == faces.length, "wrong number of faces in master grid");
        assert(other.nfaces == other.faces.length, "wrong number of faces in other grid");
        assert(ncells == cells.length, "wrong number of cells in master grid");
        assert(other.ncells == other.cells.length, "wrong number of cells in other grid");
        //
        // Merge the unique points from the other grid and keep a record of where we put them.
        size_t[] new_vtx_ids;
        if (!checkBoundariesOnly) {
            // Search all vertices in grid.
            new_vtx_ids.length = other.vertices.length;
            foreach (i, vtx; other.vertices) {
                bool found = false;
                size_t jsave;
                foreach (j, v; vertices) {
                    if (approxEqualVectors(v, vtx, relTol, absTol)) {
                        found = true;
                        jsave = j; // Remember where we found it.
                        break;
                    }
                }
                if (!found) {
                    // We have a new vertex for the master grid.
                    vertices ~= vtx;
                    new_vtx_ids[i] = vertices.length - 1;
                } else {
                    // This vertex was already in master grid.
                    new_vtx_ids[i] = jsave;
                }
            }
            nvertices = vertices.length;
        } else {
            // Search only vertices on boundaries of the grids.
            //
            // First, collect vertices on boundaries that need to be compared.
            size_t[] vtx_ids_boundary;
            foreach (i,b; boundaries) {
                // For OpenFOAM 2-D meshes the top and bottom boundaries
                // (index 4 and 5 and mutiples) can be skipped because
                // we will not be joining grids across those boundaries.
                if (openFoam && dimensions == 2  && ( (i%6)==4  || (i%6)==5) ) { continue; }
                foreach (j,f; b.face_id_list) {
                    vtx_ids_boundary ~= faces[f].vtx_id_list;
                }
            }
            vtx_ids_boundary = array(vtx_ids_boundary.sort().uniq());
            //
            size_t[] vtx_ids_other_boundary;
            foreach (i,b; other.boundaries) {
                // For OpenFOAM 2-D meshes the top and bottom boundaries
                // (index 4 and 5) can be skipped.
                if (openFoam && dimensions == 2 && (i == 4 || i == 5)) { continue; }
                foreach (f; b.face_id_list) {
                    vtx_ids_other_boundary ~= other.faces[f].vtx_id_list;
                }
            }
            vtx_ids_other_boundary = array(vtx_ids_other_boundary.sort().uniq());
            //
            // Compare vertics on boundaries and add them to master grid if new.
            new_vtx_ids.length = other.vertices.length;
            foreach (i; vtx_ids_other_boundary) {
                bool found = false;
                size_t jsave;
                Vector3 vtx;
                vtx = other.vertices[i];
                foreach (j; vtx_ids_boundary) {
                    Vector3 v;
                    v = vertices[j];
                    if (approxEqualVectors(v, vtx, relTol, absTol)) {
                        found = true;
                        jsave = j; // Remember where we found it.
                        break;
                    }
                }
                if (!found) {
                    // We have a new vertex for the master grid.
                    vertices ~= vtx;
                    new_vtx_ids[i] = vertices.length - 1;
                } else {
                    // This vertex was already in master grid.
                    new_vtx_ids[i] = jsave;
                }
            }
            //
            // Add internal vertices from other to master grid
            foreach (i, vtx; other.vertices) {
                if ( !canFind(vtx_ids_other_boundary, i)) {
                    // We have a new vertex for the master grid.
                    vertices ~= vtx;
                    new_vtx_ids[i]= vertices.length - 1;
                }
            }
            nvertices = vertices.length;
        }
        //
        //
        // Make sure that we have a full dictionary of the faces in the master-grid,
        // in preparation for merging the other collection of faces.
        foreach (i, f; faces) {
            string faceTag = makeFaceTag(f.vtx_id_list);
            if (faceTag !in faceIndices) { faceIndices[faceTag] = i; }
        }
        //
        // Merge the faces from the other grid and keep a record of where, in the master list,
        // we store them.  Also decide if a face already in the master list needs to be flipped
        // to match the other description.
        //
        size_t[] new_face_ids; new_face_ids.length = other.faces.length;
        bool[] flip_new_face; flip_new_face.length = other.faces.length;
        foreach (i, f; other.faces) {
            size_t[] new_vtx_id_list;
            foreach (vid; f.vtx_id_list) { new_vtx_id_list ~= new_vtx_ids[vid]; }
            string faceTag = makeFaceTag(new_vtx_id_list);
            if (faceTag !in faceIndices) {
                faces ~= new USGFace(new_vtx_id_list);
                auto j = faces.length - 1;
                faceIndices[faceTag] = j;
                new_face_ids[i] = j;
                flip_new_face[i] = false;
            } else {
                // Face is already in collection.
                auto j = faceIndices[faceTag];
                new_face_ids[i] = j;
                flip_new_face[i] = !sameOrientation(faces[j].vtx_id_list, new_vtx_id_list);
            }
        }
        nfaces = faces.length;
        //
        // Merge all of the cells from the other grid, using vertices and faces
        // stored in the merged lists.
        //
        foreach (i, c; other.cells) {
            size_t[] new_vtx_id_list;
            foreach (vid; c.vtx_id_list) { new_vtx_id_list ~= new_vtx_ids[vid]; }
            size_t[] new_face_id_list;
            int[] new_outsign_list;
            foreach (j, fid; c.face_id_list) {
                new_face_id_list ~= new_face_ids[fid];
                if (flip_new_face[fid]) {
                    new_outsign_list ~= -c.outsign_list[j];
                } else {
                    new_outsign_list ~= c.outsign_list[j];
                }
            }
            cells ~= new USGCell(c.cell_type, i, new_vtx_id_list,
                                 new_face_id_list, new_outsign_list);
        }
        ncells = cells.length;
        //
        // At this point, we can work through the cells and update cell references
        // into the left_cell and right_cell variables for each interface.
        // We should be able to check that the outsign values are consistent.
        //
        foreach (f; faces) { f.left_cell = null; f.right_cell = null; }
        foreach (i, c; cells) {
            foreach (j, fid; c.face_id_list) {
                switch (c.outsign_list[j]) {
                case 1:
                    assert(faces[fid].left_cell is null, "Oops, seem to already have a cell on left.");
                    faces[fid].left_cell = c;
                    break;
                case -1:
                    assert(faces[fid].right_cell is null, "Oops, seem to already have a cell on right.");
                    faces[fid].right_cell = c;
                    break;
                default:
                    throw new GeometryException("Oops, we seem to have an invalid value for outsign.");
                }
            }
        }
        foreach (i, f; faces) {
            if (f.left_cell is null && f.right_cell is null) {
                writefln("Warning, face[%d] in the merged collection is unattached.", i);
            }
        }
        //
        // Sift through the original boundary sets and eliminate faces
        // that have become "internal" (i.e. have non-null left_cell and right_cell).
        // Note: Modified June/2022 by NNG to remove empty FaceSets
        //
        
        BoundaryFaceSet[] new_boundaries;
        foreach (i, b; boundaries) {
            size_t[] new_face_id_list;
            int[] new_outsign_list;
            foreach (j, fid; b.face_id_list) {
                auto f = faces[fid];
                if (f.left_cell !is null && f.right_cell !is null) {
                    // Omit face from boundary set.
                } else {
                    new_face_id_list ~= fid;
                    new_outsign_list ~= b.outsign_list[j];
                }
            }
            if (new_face_id_list.length!=0) {
                new_boundaries ~= new BoundaryFaceSet(b.tag, new_face_id_list, new_outsign_list);
            }
        }
        boundaries = new_boundaries;

        //
        // Merge the other boundary sets into the master collection.
        //
        foreach (i, b; other.boundaries) {
            size_t[] new_face_id_list;
            int[] new_outsign_list;
            foreach (j, fid; b.face_id_list) {
                // Remember that fid is the original id in the other grid and
                // that we need to merge into the master collection with the new id.
                auto f = faces[new_face_ids[fid]];
                if (f.left_cell !is null && f.right_cell !is null) {
                    // Omit face from boundary set.
                } else {
                    // Add the face to the set using the new face id.
                    new_face_id_list ~= new_face_ids[fid];
                    if (flip_new_face[fid]) {
                        new_outsign_list ~= -b.outsign_list[j];
                    } else {
                        new_outsign_list ~= b.outsign_list[j];
                    }
                }
            }
            if (new_face_id_list.length!=0){
                boundaries ~= new BoundaryFaceSet(b.tag, new_face_id_list, new_outsign_list);
            }
        }
        // Need to update nboundaries after adding new boundaries
        nboundaries = boundaries.length;
        return this; // allows us to chain joinGrid calls
    } // end joinGrid()

    void writeStats(double relTol=1.0e-9, double absTol=1.0e-9)
    {
        writefln("UnstructuredGrid: label=%s", label);
        writefln("  nvertices=%d nfaces=%d ncells=%d", vertices.length, faces.length, cells.length);
        foreach (i, b; boundaries) {
            writefln("  BoundarySet[%d] tag=%s nfaces=%d", i, b.tag, b.face_id_list.length);
        }
        writeln("Check that each of the cells is closed.");
        bool allOK = true;
        foreach (i, c; cells) {
            Vector3 varea = vectorAreaOfCell(c);
            if (!approxEqualVectors(varea, Vector3(0,0,0), relTol, absTol)) {
                writefln("  cell id=%d, non-zero vector area=%s", i, to!string(varea));
                allOK = false;
            }
        }
        if (allOK) { writeln("  OK"); }
        //
        writeln("Check that the region bounded by the grid is closed.");
        Vector3 varea = Vector3(0,0,0);
        foreach (b; boundaries) {
            foreach (j, fid; b.face_id_list) {
                varea += b.outsign_list[j] * vectorAreaOfFacet(faces[fid].vtx_id_list);
            }
        }
        if (!approxEqualVectors(varea, Vector3(0,0,0), relTol, absTol)) {
            writeln("  Region non-zero bounding vector area=", varea);
        } else {
            writeln("  OK");
        }
        //
        return;
    } // end writeStats()

} // end class UnstructuredGrid


//-----------------------------------------------------------------

version(usgrid_test) {
    import util.msg_service;
    int main() {
        auto p00 = Vector3(0.0, 0.1);
        auto p10 = Vector3(1.0, 0.1);
        auto p11 = Vector3(1.0, 1.1);
        auto p01 = Vector3(0.0, 1.1);
        auto my_patch = new CoonsPatch(p00, p10, p11, p01);
        auto cf = [new LinearFunction(), new LinearFunction(),
                   new LinearFunction(), new LinearFunction()];
        auto my_grid = new StructuredGrid(my_patch, 11, 21, cf);
        my_grid.set_tag(Face.north, "my-special-boundary");
        string repr = to!string(my_grid);
        // writeln("my_grid=", repr);
        assert(repr.canFind("my-special-boundary"), failedUnitTest());
        assert(approxEqualVectors(*my_grid[5,5], Vector3(0.5, 0.35)),
               failedUnitTest());
        // writeln("grid point 5 5 at x=", my_grid[5,5].x, " y=", my_grid[5,5].y);
        auto usg = new UnstructuredGrid(my_grid);
        string usg_repr = to!string(usg);
        // writeln("usg=", usg_repr);
        assert(usg_repr.canFind("my-special-boundary"), failedUnitTest());
        usg.sort_cells_into_bins(10, 10);
        Vector3 my_point = Vector3(0.55, 0.325);
        size_t cell_indx = 0; bool found = false;
        usg.find_enclosing_cell(my_point, cell_indx, found);
        // writeln("Search for cell enclosing my_point= ", my_point);
        // if (found) {
        //     writeln("    cell found, index= ", cell_indx);
        //     writeln("    cell barycentre= ", usg.cell_barycentre(cell_indx));
        // } else {
        //     writeln("    cell not found");
        // }
        assert(found == true, failedUnitTest());
        assert(cell_indx == 45, failedUnitTest());
        assert(approxEqualVectors(usg.cell_barycentre(cell_indx), Vector3(0.55, 0.325)),
               failedUnitTest());
        return 0;
    }
} // end sgrid_test
