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
import std.format;
import std.math;
import gzip;

import geom;
import gpath;
import surface;
import volume;
import univariatefunctions;
import sgrid;

//-----------------------------------------------------------------

enum USGCell_type {
    none = 0,
    triangle = 1,
    quad = 2,
    polygon = 3,
    tetra = 4,
    prism = 5,
    hexahedron = 6
}
string[] cell_names = ["none", "triangle", "quad", "polygon",
		       "tetra", "prism", "hexahedron"];

int[] cell_type_for_vtk = [0, VTKElement.triangle, VTKElement.quad, VTKElement.polygon,
			   VTKElement.tetra, 0, VTKElement.hexahedron];

USGCell_type cell_type_from_name(string name)
{
    switch (name) {
    case "none": return USGCell_type.none;
    case "triangle": return USGCell_type.triangle;
    case "quad": return USGCell_type.quad;
    case "polygon": return USGCell_type.polygon;
    case "tetra": return USGCell_type.tetra;
    case "prism": return USGCell_type.prism;
    case "hexahedron": return USGCell_type.hexahedron;
    default:
	throw new Error(text("Invalid USGCell type name: ", name));
    }
}
    
class USGFace {
public:
    size_t[] vtx_id_list;

    this(size_t[] vtx_id_list) 
    {
	this.vtx_id_list = vtx_id_list.dup();
    }

    this(const USGFace other)
    {
	vtx_id_list = other.vtx_id_list.dup();
    }

    this(string str)
    {
	auto tokens = str.strip().split();
	vtx_id_list.length = to!int(tokens[0]);
	foreach(i; 0 .. vtx_id_list.length) { vtx_id_list[i] = to!int(tokens[i+1]); }
    }

    string toIOString()
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

    this(string str) 
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
		    faces ~= new USGFace([vtx_id[i][j], vtx_id[i+1][j]]);
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
		    cells ~= new USGCell(USGCell_type.quad, cell_vertices,
					 cell_faces, outsigns);
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
			auto outsigns = [+1, +1, -1, -1, +1, -1];
			cells ~= new USGCell(USGCell_type.hexahedron, cell_vertices,
					     cell_faces, outsigns);
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
	case "text":
	    read_from_text_file(fileName, false);
	    break;
	case "gziptext":
	    read_from_gzip_file(fileName);
	    break;
	case "vtk":
	    read_from_text_file(fileName, true);
	    break;
	case "vtkxml":
	    throw new Error("Reading from VTK XML format not implemented.");
	default:
	    throw new Error("Import an UnstructuredGrid, unknown format: " ~ fmt);
	}
	if (new_label != "") { label = new_label; }
    }

    UnstructuredGrid dup() const
    {
	return new UnstructuredGrid(this);
    }

    override ref Vector3 opIndex(size_t i, size_t j, size_t k=0)
    in {
	assert (i < nvertices, text("index i=", i, " is invalid, nvertices=", nvertices));
	assert (j == 0, text("index j=", j, " is invalid for unstructured grid"));
	assert (k == 0, text("index k=", k, " is invalid for unstructured grid"));
    }
    body {
	return vertices[i];
    }

    override ref Vector3 opIndex(size_t indx)
    in {
	assert (indx < nvertices,
		text("index indx=", indx, " is invalid, nvertices=", nvertices));
    }
    body {
	return vertices[indx];
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

    void read_from_text_file(string fileName, bool vtkHeader=true)
    {
	// [TODO] It may be convenient to import grids from other preprocessing
	// programs via their generic text format.
	// We've left this function here as a place-holder.
	// For VTK files, we need to work out how to extract 
	// the topological dimensions for the incoming grid.
	// Maybe all cell types being triangle or quad for 2-dimensional grids.
	string[] tokens;
	auto f = File(fileName, "r");
	if (vtkHeader) {
	    read_VTK_header_line("vtk", f);
	    label = f.readln().strip();
	    read_VTK_header_line("ASCII", f);
	    read_VTK_header_line("USTRUCTURED_GRID", f);
	} else {
	    tokens = f.readln().strip().split();
	}
	throw new Error("read_from_text_file() is not implemented, yet.");
    } // end read_from_text_file()

    override void read_from_gzip_file(string fileName)
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

    override void write_to_gzip_file(string fileName)
    // This function essentially defines the Eilmer4 native format.
    {
	auto fout = new GzipOut(fileName);
	fout.compress("unstructured_grid 1.0\n");
	fout.compress(format("label: %s\n", label));
	fout.compress(format("dimensions: %d\n", dimensions));
	fout.compress(format("vertices: %d\n", nvertices));
	foreach (v; vertices) {
	    fout.compress(format("%20.16e %20.16e %20.16e\n", v.x, v.y, v.z));
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
	foreach (v; vertices) { f.writefln("%20.16e %20.16e %20.16e", v.x, v.y, v.z); }
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

