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

class USGFace {
public:
    int[] vtx_id_list;

    this(int[] vtx_id_list) 
    {
	this.vtx_id_list = vtx_id_list.dup();
    }
} // end class USGFace

class USGCell {
public:
    int[] vtx_id_list;
    int[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this(int[] vtx_id_list, int[] face_id_list, int[] outsign_list)
    {
	this.vtx_id_list = vtx_id_list.dup();
	this.face_id_list = face_id_list.dup();
	this.outsign_list = outsign_list.dup();
    }
} // end class USGCell

class BoundaryFaceSet {
public:
    int[] face_id_list;
    int[] outsign_list; // +1 face normal is outward; -1 face normal is inward

    this() {}
} // end class BoundaryFaceSet

class UnstructuredGrid {
public:
    int dimensions;
    int nvtx, ncell, nface;
    Vector3[] vtx;
    USGFace[] face;
    USGCell[] cell;
    BoundaryFaceSet[] bndry;
    string label;

    this(const UnstructuredGrid other)
    {
	// [TODO] implement
    }

    this(const StructuredGrid sg, int dimensions)
    {
	this.dimensions = dimensions;
	label = sg.label;
	if (dimensions == 2) {
	    if (sg.nkv != 1) {
		throw new Error("UnstructuredGrid(): 2D structured grid has nkv!=1");
	    }
	    nvtx = sg.niv * sg.njv;
	    nface = (sg.niv)*(sg.njv-1) + (sg.niv-1)*(sg.njv);
	    ncell = (sg.niv-1)*(sg.njv-1);
	    bndry.length = 4; // 0=north, 1=east, 2=south, 3=west
	    // vertices
	    int[][] vtx_id;
	    vtx_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		vtx_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
		    vtx ~= Vector3(sg.grid[i][j][0]);
		    vtx_id[i][j] = vtx.length - 1;
		}
	    }
	    // i-faces
	    int[][] iface_id;
	    iface_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		iface_id[i].length = sg.njv-1;
		foreach (j; 0 .. sg.njv-1) {
		    face ~= new USGFace([vtx_id[i][j], vtx_id[i][j+1]]);
		    iface_id[i][j] = face.length - 1;
		    if (i == 0) { // west
			bndry[3].face_id_list ~= iface_id[i][j];
			bndry[3].outsign_list ~= -1;
		    }
		    if (i == sg.niv-1) { // east
			bndry[1].face_id_list ~= iface_id[i][j];
			bndry[1].outsign_list ~= +1;
		    }
		}
	    }
	    // j-faces
	    int[][] jface_id;
	    jface_id.length = sg.niv - 1;
	    foreach (i; 0 .. sg.niv-1) {
		jface_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
		    face ~= new USGFace([vtx_id[i][j], vtx_id[i+1][j]]);
		    jface_id[i][j] = face.length - 1;
		    if (j == 0) { // south
			bndry[2].face_id_list ~= jface_id[i][j];
			bndry[2].outsign_list ~= -1;
		    }
		    if (j == sg.njv-1) { // north
			bndry[0].face_id_list ~= jface_id[i][j];
			bndry[0].outsign_list ~= +1;
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
		    cell ~= new USGCell(cell_vertices, cell_faces, outsigns);
		}
	    }
	} else {
	    // Assume dimensions == 3
	    nvtx = sg.niv * sg.njv * sg.nkv;
	    nface = (sg.niv)*(sg.njv-1)*(sg.nkv-1) +
		(sg.niv-1)*(sg.njv)*(sg.nkv-1) +
		(sg.niv-1)*(sg.njv-1)*(sg.nkv);
	    ncell = (sg.niv-1)*(sg.njv-1)*(sg.nkv-1);
	    bndry.length = 6;
	    // vertices
	    int[][][] vtx_id;
	    vtx_id.length = sg.niv;
	    foreach (i; 0 .. sg.niv) {
		vtx_id[i].length = sg.njv;
		foreach (j; 0 .. sg.njv) {
		    vtx_id[i][j].length = sg.nkv;
		    foreach (k; 0 .. sg.nkv) {
			vtx ~= Vector3(sg.grid[i][j][k]);
			vtx_id[i][j][k] = vtx.length - 1;
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
			face ~= new USGFace([vtx_id[i][j][k], vtx_id[i][j+1][k],
					     vtx_id[i][j+1][k+1], vtx_id[i][j][k+1]]);
			iface_id[i][j][k] = face.length - 1;
			if (i == 0) { // west
			    bndry[3].face_id_list ~= iface_id[i][j][k];
			    bndry[3].outsign_list ~= -1;
			}
			if (i == sg.niv-1) { // east
			    bndry[1].face_id_list ~= iface_id[i][j][k];
			    bndry[1].outsign_list ~= +1;
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
			face ~= new USGFace([vtx_id[i][j][k], vtx_id[i+1][j][k],
					     vtx_id[i+1][j][k+1], vtx_id[i][j][k+1]]);
			jface_id[i][j][k] = face.length - 1;
			if (j == 0) { // south
			    bndry[2].face_id_list ~= jface_id[i][j][k];
			    bndry[2].outsign_list ~= -1;
			}
			if (j == sg.njv-1) { // north
			    bndry[0].face_id_list ~= jface_id[i][j][k];
			    bndry[0].outsign_list ~= +1;
			}
		    }
		}
	    }
	    // k-faces
	    int[][][] kface_id;
	    kface_id.length = sg.niv-1;
	    foreach (i; 0 .. sg.niv-1) {
		kface_id[i].length = sg.njv-1;
		foreach (j; 0 .. sg.njv) {
		    kface_id[i][j].length = sg.nkv;
		    foreach (k; 0 .. sg.nkv) {
			face ~= new USGFace([vtx_id[i][j][k], vtx_id[i+1][j][k],
					     vtx_id[i+1][j+1][k], vtx_id[i][j+1][k]]);
			kface_id[i][j][k] = face.length - 1;
			if (k == 0) { // bottom
			    bndry[5].face_id_list ~= kface_id[i][j][k];
			    bndry[5].outsign_list ~= -1;
			}
			if (k == sg.nkv-1) { // top
			    bndry[4].face_id_list ~= kface_id[i][j][k];
			    bndry[4].outsign_list ~= +1;
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
			cell ~= new USGCell(cell_vertices, cell_faces, outsigns);
		    }
		}
	    }
	} // end if dimensions
	// [TODO] assert numbers of faces and vertices are correct
    } // end constructor from StructuredGrid object

    // Imported grid.
    this(string fileName, int dimensions, GridFileFormat fmt, string label="empty-label")
    {
	this.dimensions = dimensions;
	this.label = label;
	final switch (fmt) {
	case GridFileFormat.text:
	    read_grid_from_text_file(fileName);
	    break;
	case GridFileFormat.gziptext:
	    read_grid_from_gzip_file(fileName);
	    break;
	case GridFileFormat.vtk:
	    read_grid_from_text_file(fileName, true);
	    break;
	case GridFileFormat.vtkxml:
	    throw new Error("Reading from VTK XML format not implemented.");
	}
    }

    UnstructuredGrid dup() const
    {
	return new UnstructuredGrid(this);
    }

    void read_grid_from_text_file(string fileName, bool vtkHeader=true)
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
    } // end read_grid_from_text_file()

    void read_grid_from_gzip_file(string fileName)
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
    } // end read_grid_from_gzip_file()

    void write_to_text_file(string fileName, bool vtkHeader=true)
    {
	auto f = File(fileName, "w");
	// if (vtkHeader) {
	//     f.writeln("# vtk DataFile Version 2.0");
	//     f.writeln(label);
	//     f.writeln("ASCII");
	//     f.writeln("");
	//     f.writeln("DATASET STRUCTURED_GRID");
	//     f.writefln("DIMENSIONS %d %d %d", niv, njv, nkv);
	//     f.writefln("POINTS %d float", (niv * njv * nkv));
	// } else {
	//     f.writefln("%d %d %d  # ni nj nk", niv, njv, nkv);
	// }
	// foreach (k; 0 .. nkv) {
	//     foreach (j; 0 .. njv) {
	// 	foreach (i; 0 .. niv) {
	// 	    f.writefln("%20.16e %20.16e %20.16e", 
	// 		       grid[i][j][k].x, grid[i][j][k].y, grid[i][j][k].z);
	// 	}
	//     }
	// }
    } // end write_to_text_file()

    void write_to_gzip_file(string fileName)
    // This function essentially defines the Eilmer4 native format.
    {
	auto f = new GzipOut(fileName);
	auto writer = appender!string();
	// formattedWrite(writer, "%d %d %d  # ni nj nk\n", niv, njv, nkv);
	// f.compress(writer.data);
	// foreach (k; 0 .. nkv) {
	//     foreach (j; 0 .. njv) {
	// 	foreach (i; 0 .. niv) {
	// 	    writer = appender!string();
	// 	    formattedWrite(writer, "%20.16e %20.16e %20.16e\n", 
	// 			   grid[i][j][k].x, grid[i][j][k].y, grid[i][j][k].z);
	// 	    f.compress(writer.data);
	// 	}
	//     }
	// }
	f.finish();
    } // end write_grid_to_gzip_file()

} // end class UnstructuredGrid

