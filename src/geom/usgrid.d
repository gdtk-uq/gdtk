/**
 * usgrid.d -- unstructured-grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-11-02 First code
 */

module usgrid;

import std.string;
import std.array;
import std.conv;
import std.stdio;
import std.algorithm;
import std.format;
import std.math;
import gzip;

import geom;
import gpath;
import surface;
import volume;
import univariatefunctions;
import grid;
import sgrid;

import paver;

//-----------------------------------------------------------------

enum USGCell_type {
    // We will use VTK names and vertex ordering, to align with su2 file format.
    none = 0,
    triangle = 1,
    quad = 2,
    polygon = 3,
    tetra = 4,
    wedge = 5,
    hexahedron = 6,
    pyramid = 7
}
string[] cell_names = ["none", "triangle", "quad", "polygon",
		       "tetra", "wedge", "hexahedron", "pyramid"];

int[] cell_type_for_vtk = [0, VTKElement.triangle, VTKElement.quad, VTKElement.polygon,
			   VTKElement.tetra, VTKElement.wedge, VTKElement.hexahedron,
			   VTKElement.pyramid];

USGCell_type cell_type_from_name(string name)
{
    switch (name) {
    case "none": return USGCell_type.none;
    case "triangle": return USGCell_type.triangle;
    case "quad": return USGCell_type.quad;
    case "polygon": return USGCell_type.polygon;
    case "tetra": return USGCell_type.tetra;
    case "wedge": return USGCell_type.wedge;
    case "hexahedron": return USGCell_type.hexahedron;
    case "pyramid": return USGCell_type.pyramid;
    default:
	throw new Error(text("Invalid USGCell type name: ", name));
    }
}
    
class USGFace {
public:
    size_t[] vtx_id_list;
    // Note that, when construction a grid from a set of vertices and lists
    // of vertex indices, we will be able to tell if a face in a boundary list
    // is facing in or out by whether it has a right or left cell.
    // This information is filled in when constructing a new mesh but is not copied
    // in the USGFace copy constructor.  Beware.
    USGCell left_cell;
    USGCell right_cell;

    this(size_t[] vtx_id_list) 
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

    string toIOString()
    // Defines the format for output to a text stream.
    {
	string str = to!string(vtx_id_list.length);
	foreach (vtx_id; vtx_id_list) { str ~= " " ~ to!string(vtx_id); }
	return str;
    }
} // end class USGFace

class USGCell {
public:
    USGCell_type cell_type;
    size_t[] vtx_id_list;
    size_t[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this(USGCell_type cell_type, size_t[] vtx_id_list,
	 size_t[] face_id_list, int[] outsign_list)
    {
	this.cell_type = cell_type;
	this.vtx_id_list = vtx_id_list.dup();
	this.face_id_list = face_id_list.dup();
	this.outsign_list = outsign_list.dup();
    }

    this(const USGCell other)
    {
	cell_type = other.cell_type;
	vtx_id_list = other.vtx_id_list.dup();
	face_id_list = other.face_id_list.dup();
	outsign_list = other.outsign_list.dup();
    }

    this(string str)
    // Defines the format for input from a text stream.
    {
	auto tokens = str.strip().split();
	int itok = 0;
	cell_type = cell_type_from_name(tokens[itok++]);
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
    }
    
    this(string str) 
    // Defines the format for input from a text stream.
    {
	auto tokens = str.strip().split();
	if (tokens.length == 1) {
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
	}
    }
    
    this(const BoundaryFaceSet other)
    {
	tag = other.tag;
	face_id_list = other.face_id_list.dup();
	outsign_list = other.outsign_list.dup();
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
} // end class BoundaryFaceSet


class UnstructuredGrid : Grid {
public:
    size_t nfaces, nboundaries;
    USGFace[] faces;
    USGCell[] cells;
    BoundaryFaceSet[] boundaries;

    //Paved Grid Constructor by Heather Muir, 2016.
    this(const Vector3[] boundary, BoundaryFaceSet[] in_boundaries, const string new_label="")
    {
    	double[][] boundary_points;
	foreach(p; boundary){
	    boundary_points ~= [p._p[0], p._p[1], p._p[2]];
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
	    this.faces ~= new USGFace(vtx_id_list);
	}
	foreach(c; CELL_LIST){
	    c.auto_cell_type();
	    if(c.cell_type != "quad"){
		throw new Error(text("paver generated a non-quad cell"));
	    } else {
		USGCell_type cell_type = cell_type_from_name(c.cell_type);
		this.cells ~= new USGCell(cell_type, c.point_IDs, c.face_IDs, c.outsigns);
	    }
	}	
    }//end paved grid constructor

    this(const UnstructuredGrid other, const string new_label="")
    {
	super(Grid_t.unstructured_grid, other.dimensions,
	      ((new_label == "")? other.label : new_label));
	nvertices = other.nvertices;
	ncells = other.ncells;
	nfaces = other.nfaces;
	nboundaries = other.nboundaries;
	foreach(v; other.vertices) { vertices ~= Vector3(v); }
	foreach(f; other.faces) { faces ~= new USGFace(f); }
	foreach(c; other.cells) { cells ~= new USGCell(c); }
	foreach(b; other.boundaries) { boundaries ~= new BoundaryFaceSet(b); }
    }

    this(const StructuredGrid sg, const string new_label="")
    {
	super(Grid_t.unstructured_grid, sg.dimensions,
	      ((new_label == "")? sg.label : new_label));
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
		boundaries ~= new BoundaryFaceSet(face_name[ib]);
		// 0=north, 1=east, 2=south, 3=west
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
	    size_t[][] iface_id;
	    iface_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		iface_id[i].length = sg.njv-1;
	    }
	    foreach (j; 0 .. sg.njv-1) {
		foreach (i; 0 .. sg.niv) {
		    faces ~= new USGFace([vtx_id[i][j], vtx_id[i][j+1]]);
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
		    faces ~= new USGFace([vtx_id[i+1][j], vtx_id[i][j]]);
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
		    auto cell_faces = [jface_id[i][j+1], // north
				       iface_id[i+1][j], // east
				       jface_id[i][j], // south
				       iface_id[i][j]]; // west
		    auto outsigns = [+1, +1, -1, -1];
		    auto my_cell = new USGCell(USGCell_type.quad, cell_vertices,
					       cell_faces, outsigns);
		    cells ~= my_cell;
		    // Now that we have the new cell, we can make the connections
		    // from the faces back to the cell.
		    faces[jface_id[i][j+1]].left_cell = my_cell;
		    faces[iface_id[i+1][j]].left_cell = my_cell;
		    faces[jface_id[i][j]].right_cell = my_cell;
		    faces[iface_id[i][j]].right_cell = my_cell;
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
		boundaries ~= new BoundaryFaceSet(face_name[ib]);
		// 0=north, 1=east, 2=south, 3=west, 4=top, 5=bottom
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
			faces ~= new USGFace([vtx_id[i][j][k], vtx_id[i][j+1][k],
					      vtx_id[i][j+1][k+1], vtx_id[i][j][k+1]]);
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
			faces ~= new USGFace([vtx_id[i][j][k], vtx_id[i+1][j][k],
					      vtx_id[i+1][j][k+1], vtx_id[i][j][k+1]]);
			jface_id[i][j][k] = faces.length - 1;
			if (j == 0) {
			    boundaries[Face.south].face_id_list ~= jface_id[i][j][k];
			    boundaries[Face.south].outsign_list ~= +1;
			}
			if (j == sg.njv-1) {
			    boundaries[Face.north].face_id_list ~= jface_id[i][j][k];
			    boundaries[Face.north].outsign_list ~= -1;
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
			faces ~= new USGFace([vtx_id[i][j][k], vtx_id[i+1][j][k],
					      vtx_id[i+1][j+1][k], vtx_id[i][j+1][k]]);
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
			auto cell_faces = [jface_id[i][j+1][k], // north
					   iface_id[i+1][j][k], // east
					   jface_id[i][j][k], // south
					   iface_id[i][j][k], // west
					   kface_id[i][j][k+1], // top
					   kface_id[i][j][k]]; // bottom
			auto outsigns = [-1, +1, +1, -1, +1, -1];
			auto my_cell = new USGCell(USGCell_type.hexahedron, cell_vertices,
						   cell_faces, outsigns);
			cells ~= my_cell;
			// Now that we have the new cell, we can make the connections
			// from the faces back to the cell.
			faces[jface_id[i][j+1][k]].left_cell = my_cell;
			faces[iface_id[i+1][j][k]].left_cell = my_cell;
			faces[jface_id[i][j][k]].right_cell = my_cell;
			faces[iface_id[i][j][k]].right_cell = my_cell;
			faces[kface_id[i][j][k+1]].left_cell = my_cell;
			faces[kface_id[i][j][k]].right_cell = my_cell;
		    }
		}
	    }
	} // end if dimensions
	assert(nvertices == vertices.length, "mismatch in number of vertices");
	assert(nfaces == faces.length, "mismatch in number of faces");
	assert(ncells == cells.length, "mismatch in number of cells");
    } // end constructor from StructuredGrid object

    // Imported grid.
    this(string fileName, string fmt, string new_label="")
    {
	super(Grid_t.unstructured_grid, 0, new_label);
	// dimensions will be reset on reading grid
	switch (fmt) {
	case "gziptext": read_from_gzip_file(fileName); break;
	case "su2text": read_from_su2_text_file(fileName); break;
	case "vtktext": read_from_vtk_text_file(fileName); break;
	case "vtkxml": throw new Error("Reading from VTK XML format not implemented.");
	default: throw new Error("Import an UnstructuredGrid, unknown format: " ~ fmt);
	}
	if (new_label != "") { label = new_label; }
    }

    UnstructuredGrid dup() const
    {
	return new UnstructuredGrid(this);
    }

    // -----------------------------
    // Indexing and location methods.
    // -----------------------------
    
    override Vector3* opIndex(size_t i, size_t j, size_t k=0)
    in {
	assert (i < nvertices, text("index i=", i, " is invalid, nvertices=", nvertices));
	assert (j == 0, text("index j=", j, " is invalid for unstructured grid"));
	assert (k == 0, text("index k=", k, " is invalid for unstructured grid"));
    }
    body {
	return &(vertices[i]);
    }

    override Vector3* opIndex(size_t indx)
    in {
	assert (indx < nvertices,
		text("index indx=", indx, " is invalid, nvertices=", nvertices));
    }
    body {
	return &(vertices[indx]);
    }

    override size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const
    in {
	assert (i < ncells, text("index i=", i, " is invalid, ncells=", ncells));
	assert (j == 0, text("index j=", j, " is invalid for unstructured grid"));
	assert (k == 0, text("index k=", k, " is invalid for unstructured grid"));
    }
    body {
	return cells[i].vtx_id_list.dup();
    }

    override size_t[] get_vtx_id_list_for_cell(size_t indx) const
    in {
	assert (indx < ncells,
		text("index indx=", indx, " is invalid, ncells=", ncells));
    }
    body {
	return cells[indx].vtx_id_list.dup();
    }

    bool findFaceIndex(const size_t[] id_list, ref size_t indx)
    // Given a list of vertex ids, determine is we already have a face
    // constructed from those vertices, and returns its index.
    //
    // This is a slow, brute-force search with lots of allocations.
    // If we have to use it a lot, if may be better to store hashes
    // of the list of vertex-id values at face construction time and
    // search on those.
    {
	bool found = false;
	size_t[] sorted_id_list = id_list.dup();
	foreach(i; 0 .. faces.length) {
	    if (faces[i].vtx_id_list.length != sorted_id_list.length) continue;
	    size_t[] sorted_list_for_face = faces[i].vtx_id_list.dup();
	    if (equal(sorted_id_list, sorted_list_for_face)) {
		found = true; indx = i;
	    }
	}
	return found;
    }
    
    // ------------------------
    // Import-from-file methods.
    // ------------------------
    
    override void read_from_gzip_file(string fileName)
    // This function, together with the constructors (from strings) for
    // the classes USGFace, USGCell and BoundaryFaceSet (above),
    // define the Eilmer4 native format.
    // For output, there is the matching write_to_gzip_file() below.
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
	    vertices ~= Vector3(x, y, z);
	}
	line = byLine.front; byLine.popFront();
	formattedRead(line, "faces: %d", &nfaces);
	faces.length = 0;
	foreach (i; 0 .. nfaces) {
	    line = byLine.front; byLine.popFront();
	    faces ~= new USGFace(line);
	}
	line = byLine.front; byLine.popFront();
	formattedRead(line, "cells: %d", &ncells);
	cells.length = 0;
	foreach (i; 0 .. ncells) {
	    line = byLine.front; byLine.popFront();
	    cells ~= new USGCell(line);
	}
	line = byLine.front; byLine.popFront();
	formattedRead(line, "boundaries: %d", &nboundaries);
	boundaries.length = 0;
	foreach (i; 0 .. nboundaries) {
	    line = byLine.front; byLine.popFront();
	    boundaries ~= new BoundaryFaceSet(line);
	}
    } // end read_from_gzip_file()

    void read_from_su2_text_file(string fileName)
    // Information on the su2 file format from
    // https://github.com/su2code/SU2/wiki/Mesh-File
    {
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
	writeln("dimensions=", dimensions);
	ncells = to!size_t(getHeaderContent("NELEM"));
	writeln("ncells=", ncells);
	cells.length = ncells;
	foreach(i; 0 .. ncells) {
	    auto lineContent = f.readln().strip();
	    auto tokens = lineContent.split();
	    int vtk_element_type = to!int(tokens[0]);
	    size_t indx = to!size_t(tokens[$-1]);
	    size_t[] vtx_id_list;
	    foreach(j; 1 .. tokens.length-1) { vtx_id_list ~= to!size_t(tokens[j]); }
	    USGCell_type cell_type = USGCell_type.none;
	    switch (vtk_element_type) {
	    case VTKElement.triangle: cell_type = USGCell_type.triangle; break;
	    case VTKElement.quad: cell_type = USGCell_type.quad; break;
	    case VTKElement.tetra: cell_type = USGCell_type.tetra; break;
	    case VTKElement.wedge: cell_type = USGCell_type.wedge; break;
	    case VTKElement.hexahedron: cell_type = USGCell_type.hexahedron; break;
	    case VTKElement.pyramid: cell_type = USGCell_type.pyramid; break;
	    default:
		throw new Exception("unknown element type for line: "~to!string(lineContent));
	    }
	    size_t[] face_id_list; // empty, so far
	    int[] outsign_list; // empty, so far
	    // Once we have created the full list of vertices, we will be able to make
	    // the set of faces and then come back to filling in the face_id_list and
	    // outsign_list for each cell.
	    cells[indx] = new USGCell(cell_type, vtx_id_list, face_id_list, outsign_list);
	} // end foreach i .. ncells
	foreach(i; 0 .. cells.length) {
	    if (!cells[i]) { writeln("Warning: uninitialized cell at index: ", i); }
	    // [TODO] if we have any uninitialized cells,
	    // we should compress the array to eliminate empty elements.
	}
	nvertices = to!size_t(getHeaderContent("NPOIN"));
	writeln("nvertices=", nvertices);
	vertices.length = nvertices;
	foreach(i; 0 .. nvertices) {
	    auto tokens = f.readln().strip().split();
	    double x=0.0; double y=0.0; double z = 0.0; size_t indx = 0;
	    if (dimensions == 2) {
		x = to!double(tokens[0]);
		y = to!double(tokens[1]);
		indx = to!size_t(tokens[2]);
	    } else {
		assert(dimensions == 3, "invalid dimensions");
		x = to!double(tokens[0]);
		y = to!double(tokens[1]);
		z = to!double(tokens[2]);
		indx = to!size_t(tokens[3]);
	    }
	    Vector3* vtx = &vertices[indx];
	    vtx.refx = x; vtx.refy = y; vtx.refz = z;
	    // writeln("indx=", indx, " vtx=", *vtx);
	} // end foreach i .. nvertices
	//
	// Now that we have the full list of cells and vertices assigned to each cell,
	// we can construct the faces between cells and along the boundaries.
	//
	void add_linear_face_for_2D(ref USGCell cell, size_t idx0, size_t idx1)
	{
	    size_t[] my_vtx_id_list = [cell.vtx_id_list[idx0], cell.vtx_id_list[idx1]];
	    size_t face_indx = 0;
	    if (!findFaceIndex(my_vtx_id_list, face_indx)) {
		// Since we didn't find the face already, construct it.
		face_indx = faces.length;
		faces ~= new USGFace(my_vtx_id_list);
	    }
	    cell.face_id_list ~= face_indx;
	    // Note that, from this point, we work with the face from the collection
	    // which, if it was an existing face, will not have the same orientation
	    // as indicated by the cycle of vertices passed in to this function call.
	    Vector3* v0 = &vertices[faces[face_indx].vtx_id_list[0]];
	    Vector3* v1 = &vertices[faces[face_indx].vtx_id_list[1]];
	    Vector3 pmid = Vector3(vertices[cell.vtx_id_list[0]]);
	    foreach(i; 1 .. cell.vtx_id_list.length) {
		pmid += vertices[cell.vtx_id_list[i]];
	    }
	    pmid /= cell.vtx_id_list.length;
	    if (on_left_of_xy_line(*v0, *v1, pmid)) {
		if (faces[face_indx].left_cell) {
		    throw new Exception("face already has a left cell");
		}
		faces[face_indx].left_cell = cell;
		cell.outsign_list ~=  1;
	    } else {
		if (faces[face_indx].right_cell) {
		    throw new Exception("face already has a right cell");
		}
		faces[face_indx].right_cell = cell;
		cell.outsign_list ~= -1;
	    }
	    return;
	} // end add_linear_face_for_2D()
	void add_quadrilateral_face_for_3D(ref USGCell cell, size_t idx0, size_t idx1,
					   size_t idx2, size_t idx3)
	{
	    size_t[] my_vtx_id_list = [cell.vtx_id_list[idx0], cell.vtx_id_list[idx1],
				       cell.vtx_id_list[idx2], cell.vtx_id_list[idx3]];
	    size_t face_indx = 0;
	    if (!findFaceIndex(my_vtx_id_list, face_indx)) {
		// Since we didn't find the face already, construct it.
		face_indx = faces.length;
		faces ~= new USGFace(my_vtx_id_list);
	    }
	    cell.face_id_list ~= face_indx;
	    // Note that, from this point, we work with the face from the collection
	    // which, if it was an existing face, will not have the same orientation
	    // as indicated by the cycle of vertices passed in to this function call.
	    Vector3* v0 = &vertices[faces[face_indx].vtx_id_list[0]];
	    Vector3* v1 = &vertices[faces[face_indx].vtx_id_list[1]];
	    Vector3* v2 = &vertices[faces[face_indx].vtx_id_list[2]];
	    Vector3* v3 = &vertices[faces[face_indx].vtx_id_list[3]];
	    Vector3 vmid = 0.25*(*v0 + *v1 + *v2 + *v3);
	    Vector3 cell_mid = Vector3(vertices[cell.vtx_id_list[0]]);
	    foreach(i; 1 .. cell.vtx_id_list.length) {
		cell_mid += vertices[cell.vtx_id_list[i]];
	    }
	    cell_mid /= cell.vtx_id_list.length;
	    if (tetragonal_dipyramid_volume(*v0, *v1, *v2, *v3, vmid, cell_mid) < 0.0) {
		if (faces[face_indx].left_cell) {
		    throw new Exception("face already has a left cell");
		}
		faces[face_indx].left_cell = cell;
		cell.outsign_list ~=  1;
	    } else {
		if (faces[face_indx].right_cell) {
		    throw new Exception("face already has a right cell");
		}
		faces[face_indx].right_cell = cell;
		cell.outsign_list ~= -1;
	    }
	    return;
	} // end add_quadrilateral_face_for_3D()
	foreach(cell; cells) {
	    if (!cell) continue;
	    if (dimensions == 2) {
		switch(cell.cell_type) {
		case USGCell_type.triangle:
		    add_linear_face_for_2D(cell, 0, 1);
		    add_linear_face_for_2D(cell, 1, 2);
		    add_linear_face_for_2D(cell, 2, 0);
		    break;
		case USGCell_type.quad:
		    add_linear_face_for_2D(cell, 2, 3); // north
		    add_linear_face_for_2D(cell, 1, 2); // east
		    add_linear_face_for_2D(cell, 0, 1); // south
		    add_linear_face_for_2D(cell, 3, 0); // west
		    break;
		default:
		    throw new Exception("invalid cell type in 2D");
		}
	    } else {
		assert(dimensions == 3, "invalid dimensions");
		switch(cell.cell_type) {
		case USGCell_type.tetra:
		    throw new Exception("not implemented yet");
		case USGCell_type.hexahedron:
		    add_quadrilateral_face_for_3D(cell, 2, 3, 7, 6); // north
		    add_quadrilateral_face_for_3D(cell, 1, 2, 6, 5); // east
		    add_quadrilateral_face_for_3D(cell, 1, 0, 4, 5); // south
		    add_quadrilateral_face_for_3D(cell, 0, 3, 7, 4); // west
		    add_quadrilateral_face_for_3D(cell, 4, 5, 6, 7); // top
		    add_quadrilateral_face_for_3D(cell, 0, 1, 2, 3); // bottom
		    break;
		case USGCell_type.wedge:
		    throw new Exception("not implemented yet");
		case USGCell_type.pyramid:
		    throw new Exception("not implemented yet");
		default:
		    throw new Exception("invalid cell type in 3D");
		}
	    }
	} // end foreach cell
	//
	// Now that we have a full set of cells and faces,
	// make lists of the boundary faces.
	//
	nboundaries = to!size_t(getHeaderContent("NMARK"));
	foreach(i; 0 .. nboundaries) {
	    string tag = getHeaderContent("MARKER_TAG");
	    writeln("boundary i=", i, " tag=", tag);
	    size_t[] face_id_list;
	    int[] outsign_list;
	    size_t nelem = to!size_t(getHeaderContent("MARKER_ELEMS"));
	    // writeln("nelem=", nelem);
	    foreach(j; 0 .. nelem) {
		auto tokens = f.readln().strip().split();
		int vtk_type = to!int(tokens[0]);
		size_t[] my_vtx_id_list;
		foreach(k; 1 .. tokens.length) { my_vtx_id_list ~= to!size_t(tokens[k]); }
		size_t face_indx = 0;
		if (findFaceIndex(my_vtx_id_list, face_indx)) {
		    USGFace my_face = faces[face_indx];
		    assert(my_face.left_cell || my_face.right_cell,
			   "face is not properly connected");
		    if (my_face.left_cell && !(my_face.right_cell)) {
			outsign_list ~= 1;
		    } else if ((!my_face.left_cell) && my_face.right_cell) {
			outsign_list ~= -1;
		    } else {
			throw new Exception("appears to be an interior face");
		    }
		    face_id_list ~= face_indx;
		} else {
		    throw new Exception("cannot find face in collection");
		}
	    } // end foreach j
	    boundaries ~= new BoundaryFaceSet(tag, face_id_list, outsign_list);
	} // end foreach i
    } // end read_from_su2_text_file()

    void read_from_vtk_text_file(string fileName)
    {
	string[] tokens;
	auto f = File(fileName, "r");
	read_VTK_header_line("vtk", f);
	label = f.readln().strip();
	read_VTK_header_line("ASCII", f);
	read_VTK_header_line("USTRUCTURED_GRID", f);
	tokens = f.readln().strip().split();
	throw new Error("read_from_vtk_text_file() is not finished, yet.");
	// For VTK files, we need to work out how to extract 
	// the topological details for the incoming grid.
    } // end read_from_vtk_text_file()

    override void write_to_gzip_file(string fileName)
    // This function, together with the toIOstring methods for
    // the classes USGFace, USGCell and BoundaryFaceSet (way above)
    // define the Eilmer4 native format.
    // For input, there is the matching read_from_gzip_file(), above.
    {
	auto fout = new GzipOut(fileName);
	fout.compress("unstructured_grid 1.0\n");
	fout.compress(format("label: %s\n", label));
	fout.compress(format("dimensions: %d\n", dimensions));
	fout.compress(format("vertices: %d\n", nvertices));
	foreach (v; vertices) {
	    fout.compress(format("%.18e %.18e %.18e\n", v.x, v.y, v.z));
	}
	fout.compress(format("faces: %d\n", nfaces));
	foreach (f; faces) { fout.compress(f.toIOString ~ "\n"); }
	fout.compress(format("cells: %d\n", ncells));
	foreach (c; cells) { fout.compress(c.toIOString ~ "\n"); }
	fout.compress(format("boundaries: %d\n", boundaries.length));
	foreach (b; boundaries) { fout.compress(b.toIOString ~ "\n"); }
	fout.finish();
    } // end write_to_gzip_file()

    override void write_to_vtk_file(string fileName)
    {
	auto f = File(fileName, "w");
	f.writeln("# vtk DataFile Version 2.0");
	f.writeln(label);
	f.writeln("ASCII");
	f.writeln("");
	f.writeln("DATASET UNSTRUCTURED_GRID");
	f.writefln("POINTS %d float", nvertices);
	foreach (v; vertices) { f.writefln("%.18e %.18e %.18e", v.x, v.y, v.z); }
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
	    f.writeln(cell_type_for_vtk[c.cell_type]);
	}
	f.close();
    } // end write_to_vtk_file()

} // end class UnstructuredGrid



