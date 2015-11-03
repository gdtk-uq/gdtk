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
    int[] vtx_id_list;

    this(int[] vtx_id_list) 
    {
	this.vtx_id_list = vtx_id_list.dup();
    }

    this(const USGFace other)
    {
	vtx_id_list = other.vtx_id_list.dup();
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
    int[] vtx_id_list;
    int[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this(USGCell_type cell_type, int[] vtx_id_list,
	 int[] face_id_list, int[] outsign_list)
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
    int[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this(string tag) 
    {
	this.tag = tag;
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

class UnstructuredGrid {
public:
    int dimensions;
    int nvertices, ncells, nfaces;
    Vector3[] vertices;
    USGFace[] faces;
    USGCell[] cells;
    BoundaryFaceSet[] boundaries;
    string label;

    this(const UnstructuredGrid other)
    {
	dimensions = other.dimensions;
	nvertices = other.nvertices;
	ncells = other.ncells;
	nfaces = other.nfaces;
	foreach(v; other.vertices) { vertices ~= Vector3(v); }
	foreach(f; other.faces) { faces ~= new USGFace(f); }
	foreach(c; other.cells) { cells ~= new USGCell(c); }
	foreach(b; other.boundaries) { boundaries ~= new BoundaryFaceSet(b); }
	label = other.label;
    }

    this(const StructuredGrid sg, int dimensions)
    {
	this.dimensions = dimensions;
	label = sg.label;
	if (dimensions == 2) {
	    if (sg.nkv != 1) {
		throw new Error("UnstructuredGrid(): 2D structured grid has nkv!=1");
	    }
	    nvertices = sg.niv * sg.njv;
	    nfaces = (sg.niv)*(sg.njv-1) + (sg.niv-1)*(sg.njv);
	    ncells = (sg.niv-1)*(sg.njv-1);
	    foreach(ib; 0 .. 4) { boundaries ~= new BoundaryFaceSet(face_name[ib]); }
	    // 0=north, 1=east, 2=south, 3=west
	    // vertices
	    int[][] vtx_id;
	    vtx_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		vtx_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
		    vertices ~= Vector3(sg.grid[i][j][0]);
		    vtx_id[i][j] = vertices.length - 1;
		}
	    }
	    // i-faces
	    int[][] iface_id;
	    iface_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		iface_id[i].length = sg.njv-1;
		foreach (j; 0 .. sg.njv-1) {
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
	    int[][] jface_id;
	    jface_id.length = sg.niv - 1;
	    foreach (i; 0 .. sg.niv-1) {
		jface_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
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
	    foreach (i; 0 .. sg.niv-1) {
		foreach (j; 0 .. sg.njv-1) {
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
	    foreach(ib; 0 .. 6) { boundaries ~= new BoundaryFaceSet(face_name[ib]); }
	    // 0=north, 1=east, 2=south, 3=west, 4=top, 5=bottom
	    // vertices
	    int[][][] vtx_id;
	    vtx_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		vtx_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
		    vtx_id[i][j].length = sg.nkv;
		    foreach (k; 0 .. sg.nkv) {
			vertices ~= Vector3(sg.grid[i][j][k]);
			vtx_id[i][j][k] = vertices.length - 1;
		    }
		}
	    }
	    // i-faces
	    int[][][] iface_id;
	    iface_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		iface_id[i].length = sg.njv-1;
		foreach (j; 0 .. sg.njv-1) {
		    iface_id[i][j].length = sg.nkv-1;
		    foreach (k; 0 .. sg.nkv-1) {
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
	    int[][][] jface_id;
	    jface_id.length = sg.niv-1;
	    foreach (i; 0 .. sg.niv-1) {
		jface_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
		    jface_id[i][j].length = sg.nkv-1;
		    foreach (k; 0 .. sg.nkv-1) {
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
	    int[][][] kface_id;
	    kface_id.length = sg.niv-1;
	    foreach (i; 0 .. sg.niv-1) {
		kface_id[i].length = sg.njv-1;
		foreach (j; 0 .. sg.njv-1) {
		    kface_id[i][j].length = sg.nkv;
		    foreach (k; 0 .. sg.nkv) {
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
	    foreach (i; 0 .. sg.niv-1) {
		foreach (j; 0 .. sg.njv-1) {
		    foreach (k; 0 .. sg.nkv-1) {
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
    this(string fileName, int dimensions, GridFileFormat fmt, string label="empty-label")
    {
	this.dimensions = dimensions;
	this.label = label;
	final switch (fmt) {
	case GridFileFormat.text:
	    read_from_text_file(fileName);
	    break;
	case GridFileFormat.gziptext:
	    read_from_gzip_file(fileName);
	    break;
	case GridFileFormat.vtk:
	    read_from_text_file(fileName, true);
	    break;
	case GridFileFormat.vtkxml:
	    throw new Error("Reading from VTK XML format not implemented.");
	}
    }

    UnstructuredGrid dup() const
    {
	return new UnstructuredGrid(this);
    }

    void read_from_text_file(string fileName, bool vtkHeader=true)
    {
	string[] tokens;
	auto f = File(fileName, "r");
	if (vtkHeader) {
	    read_VTK_header_line("vtk", f);
	    label = f.readln().strip();
	    read_VTK_header_line("ASCII", f);
	    read_VTK_header_line("STRUCTURED_GRID", f);
	    tokens = read_VTK_header_line("DIMENSIONS", f);
	} else {
	    tokens = f.readln().strip().split();
	}
	// niv = to!int(tokens[0]);
	// njv = to!int(tokens[1]);
	// nkv = to!int(tokens[2]);
	// resize_array();
	// foreach (k; 0 .. nkv) {
	//     foreach (j; 0 .. njv) {
	// 	foreach (i; 0 .. niv) {
	// 	    tokens = f.readln().strip().split();
	// 	    try {
	// 		grid[i][j][k].refx = to!double(tokens[0]);
	// 		grid[i][j][k].refy = to!double(tokens[1]);
	// 		grid[i][j][k].refz = to!double(tokens[2]);
	// 	    } catch (Exception e) {
	// 		throw new Error(text("Failed to read grid file at "
	// 				     "i=", i, " j=", j, " k=", k,
	// 				     "tokens=", tokens, "exception=", e));
	// 	    }
	// 	} // foreach i
	//     } // foreach j
	// } // foreach k
    } // end read_from_text_file()

    void read_from_gzip_file(string fileName)
    {
	auto byLine = new GzipByLine(fileName);
	auto line = byLine.front; byLine.popFront();
	// formattedRead(line, "%d %d %d", &niv, &njv, &nkv);
	// resize_array();
	// double x, y, z;
	// foreach (k; 0 .. nkv) {
	//     foreach (j; 0 .. njv) {
	// 	foreach (i; 0 .. niv) {
	// 	    line = byLine.front; byLine.popFront();
	// 	    // Note that the line starts with whitespace.
	// 	    formattedRead(line, " %g %g %g", &x, &y, &z);
	// 	    grid[i][j][k].refx = x;
	// 	    grid[i][j][k].refy = y;
	// 	    grid[i][j][k].refz = z;
	// 	} // foreach i
	//     } // foreach j
	// } // foreach k
    } // end read_from_gzip_file()

    void write_to_vtk_file(string fileName)
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

    void write_to_gzip_file(string fileName)
    // This function essentially defines the Eilmer4 native format.
    {
	auto fout = new GzipOut(fileName);
	fout.compress(format("label %s\n", label));
	fout.compress(format("vertices %d\n", nvertices));
	foreach (v; vertices) {
	    fout.compress(format("%20.16e %20.16e %20.16e\n", v.x, v.y, v.z));
	}
	fout.compress(format("faces %d\n", nfaces));
	foreach (f; faces) { fout.compress(f.toIOString ~ "\n"); }
	fout.compress(format("cells %d\n", ncells));
	foreach (c; cells) { fout.compress(c.toIOString ~ "\n"); }
	fout.compress(format("boundaries %d\n", boundaries.length));
	foreach (b; boundaries) { fout.compress(b.toIOString ~ "\n"); }
	fout.finish();
    } // end write_to_gzip_file()

} // end class UnstructuredGrid

