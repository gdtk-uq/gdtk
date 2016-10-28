/**
 * sgrid.d -- structured-grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-02-22 First code ported from e3_grid.py
 */

module sgrid;

import std.algorithm;
import std.string;
import std.array;
import std.conv;
import std.stdio;
import std.format;
import std.math;
import gzip;

import grid;
import geom;
import gpath;
import surface;
import volume;
import univariatefunctions;

//-----------------------------------------------------------------

class StructuredGrid : Grid {
public:
    size_t[] vtx_id; // used to hold the single-index id for each vertex.
    // So that the structured-grid can look like an unstructured grid.

    // Blank grid, ready for import of data.
    this(size_t niv, size_t njv, size_t nkv=1, string label="")
    in {
	assert(niv > 0 && njv > 0 && nkv > 0);
    }
    body {
	// infer dimensions from numbers of vertices in each direction
	int dim = 3;
	if (nkv == 1) {
	    if (njv == 1) { dim = 1; } else { dim = 2; }
	} 
	super(Grid_t.structured_grid, dim, label);
	this.niv = niv; this.njv = njv; this.nkv = nkv;
	switch (dim) {
	case 1: ncells = niv-1; break;
	case 2: ncells = (niv-1)*(njv-1); break;
	case 3: ncells = (niv-1)*(njv-1)*(nkv-1); break;
	default: assert(0);
	}
	nvertices = niv*njv*nkv;
	vertices.length = nvertices;
	vtx_id.length = nvertices;
	// Standard order of vertices.
	size_t ivtx = 0;
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    vtx_id[single_index(i,j,k)] = ivtx;
		    ivtx++;
		}
	    }
	}
    } // end this()

    // 1D grid, built on Path.
    this(const Path my_path, size_t niv,
	 const(UnivariateFunction) clusterf = new LinearFunction(0.0, 1.0),
	 string label="")
    {
	this(niv, 1, 1, label);
	make_grid_from_path(my_path, clusterf);
    }

    // 2D grid, built on Parametric surface.
    this(const ParametricSurface surf, size_t niv, size_t njv,
	 const(UnivariateFunction)[] clusterf, string label="")
    {
	this(niv, njv, 1, label);
	// Any unspecified clustering functions default to the linear identity.
	while (clusterf.length < 4) clusterf ~= new LinearFunction(0.0, 1.0);
	make_grid_from_surface(surf, clusterf);
    }

    // 3D grid, built on ParametricVolume.
    this(const ParametricVolume pvolume, size_t niv, size_t njv, size_t nkv,
	 const(UnivariateFunction)[] clusterf, string label="")
    {
	this(niv, njv, nkv, label);
	// Any unspecified clustering functions default to the linear identity.
	while (clusterf.length < 12) clusterf ~= new LinearFunction(0.0, 1.0);
	make_grid_from_volume(pvolume, clusterf);
    }

    // Imported grid.
    this(string fileName, string fmt, string label="")
    {
	this(1, 1, 1, label); // these settings will be reset on actually reading the data file
	switch (fmt) {
	case "text":
	    read_from_text_file(fileName);
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
	    throw new Error("Reading StructuredGrid, unknown format: " ~ fmt);
	}
    }

    this(const StructuredGrid other)
    {
	this(other.niv, other.njv, nkv = other.nkv, other.label);
	foreach (i; 0 .. vertices.length) {
	    vertices[i].refx = other.vertices[i].x;
	    vertices[i].refy = other.vertices[i].y;
	    vertices[i].refz = other.vertices[i].z;
	}
    }

    StructuredGrid dup() const
    {
	return new StructuredGrid(this);
    }
	
    override Vector3* opIndex(size_t i, size_t j, size_t k=0)
    in {
	assert (i < niv, text("index i=", i, " is invalid, niv=", niv));
	assert (j < njv, text("index j=", j, " is invalid, njv=", njv));
	assert (k < nkv, text("index k=", k, " is invalid, nkv=", nkv));
    }
    body {
	return &(vertices[single_index(i,j,k)]);
    }

    override Vector3* opIndex(size_t indx)
    in {
	assert (indx < niv*njv*nkv,
		text("index indx=", indx, " is invalid, niv*njv*nkv=", niv*njv*nkv));
    }
    body {
	return &(vertices[indx]);
    }

    override size_t number_of_vertices_for_cell(size_t i)
    {
	switch (dimensions) {
	case 1: return 2;
	case 2: return 4;
	case 3: return 8;
	default: return 0; // oops
	}
    }

    override int vtk_element_type_for_cell(size_t i)
    {
	switch (dimensions) {
	case 1: return VTKElement.line;
	case 2: return VTKElement.quad;
	case 3: return VTKElement.hexahedron;
	default: return 0; // oops
	}
    }
    
    override size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const
    in {
	size_t nic, njc, nkc;
	switch (dimensions) {
	case 1:
	    nic = niv - 1;
	    assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
	    assert (j == 0, text("index j=", j, " is invalid for 1D grid"));
	    assert (k == 0, text("index k=", k, " is invalid for 1D grid"));
	    break;
	case 2:
	    nic = niv - 1; njc = njv - 1;
	    assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
	    assert (j < njc, text("index j=", j, " is invalid, njc=", njc));
	    assert (k == 0, text("index k=", k, " is invalid for 2D grid"));
	    break;
	case 3:
	    nic = niv - 1; njc = njv - 1; nkc = nkv - 1;
	    assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
	    assert (j < njc, text("index j=", j, " is invalid, njc=", njc));
	    assert (k < nkc, text("index k=", k, " is invalid, nkc=", nkc));
	    break;
	default: assert(0);
	}
    }
    body {
	switch (dimensions) {
	case 1: return [vtx_id[single_index(i,j,k)],
			vtx_id[single_index(i+1,j,k)]];
	case 2: return [vtx_id[single_index(i,j,k)],
			vtx_id[single_index(i+1,j,k)],
			vtx_id[single_index(i+1,j+1,k)],
			vtx_id[single_index(i,j+1,k)]];
	case 3: return [vtx_id[single_index(i,j,k)],
			vtx_id[single_index(i+1,j,k)], 
			vtx_id[single_index(i+1,j+1,k)],
			vtx_id[single_index(i,j+1,k)],
			vtx_id[single_index(i,j,k+1)],
			vtx_id[single_index(i+1,j,k+1)], 
			vtx_id[single_index(i+1,j+1,k+1)],
			vtx_id[single_index(i,j+1,k+1)]];
	default: assert(0);
	}
    }

    override size_t[] get_vtx_id_list_for_cell(size_t indx) const
    {
	size_t nic, njc, nkc;
	switch (dimensions) {
	case 1: nic = niv-1; njc = 1; nkc = 1; break;
	case 2: nic = niv-1; njc = njv-1; nkc = 1; break;
	case 3: nic = niv-1; njc = njv-1; nkc = njv-1; break;
	default: assert(0);
	}
	size_t k = indx / (nic*njc);
	indx -= k * (nic * njc);
	size_t j = indx / nic;
	indx -= j * nic;
	return get_vtx_id_list_for_cell(indx, j, k);
    }

    StructuredGrid subgrid(size_t i0, size_t ni,
			   size_t j0=0, size_t nj=1,
			   size_t k0=0, size_t nk=1) const
    // Partition on vertex indices.
    {
	if (i0+ni > niv)
	    throw new Error(text("Sgrid.subgrid overrun i0=",i0,", ni=",ni,", niv=",niv));
	if (j0+nj > njv)
	    throw new Error(text("Sgrid.subgrid overrun j0=",j0,", nj=",nj,", njv=",njv));
	if (k0+nk > nkv)
	    throw new Error(text("Sgrid.subgrid overrun k0=",k0,", nk=",nk,", nkv=",nkv));
	auto new_grd = new StructuredGrid(ni, nj, nk);
	foreach (i; 0 .. ni) {
	    foreach (j; 0 .. nj) {
		foreach (k; 0 .. nk) {
		    new_grd[i,j,k].set(vertices[single_index(i0+i,j0+j,k0+k)]);
		}
	    }
	}
	return new_grd;
    } // end subgrid()

    override Grid get_boundary_grid(size_t boundary_indx)
    // Returns the grid defining a particular boundary of the original grid.
    // For an 3D block, a 2D surface grid will be returned, with index directions
    // as defined on the debugging cube.
    {
	size_t new_niv, new_njv, new_nkv;
	string new_label = label;
	if (dimensions == 3) {
	    // 1. Use orientation to decide size of new 2D grid.
	    // When looking at the 2D grid alone, the new i-direcion starts west and
	    // progresses to east, while the new j-direction starts south and progresses north.
	    final switch (boundary_indx) {
	    case Face.north:
		new_niv = niv; new_njv = nkv; new_nkv = 1; new_label ~= "-north"; break;
	    case Face.south:
		new_niv = niv; new_njv = nkv; new_nkv = 1; new_label ~= "-south"; break;
	    case Face.west:
		new_niv = njv; new_njv = nkv; new_nkv = 1; new_label ~= "-west"; break;
	    case Face.east:
		new_niv = njv; new_njv = nkv; new_nkv = 1; new_label ~= "-east"; break;
	    case Face.top:
		new_niv = niv; new_njv = njv; new_nkv = 1; new_label ~= "-top"; break;
	    case Face.bottom:
		new_niv = niv; new_njv = njv; new_nkv = 1; new_label ~= "-bottom"; break;
	    } // end switch boundary_indx
	} else {
	    throw new Exception("Extraction from 2D grid not implemented yet.");
	}
	// 2. prepare the empty grid.
	auto bgrid = new StructuredGrid(new_niv, new_njv, new_nkv, label);
	// 3. and fill it.
	if (dimensions == 3) {
	    final switch (boundary_indx) {
	    case Face.north:
		foreach (i; 0 .. new_niv) {
		    foreach (j; 0 .. new_njv) { bgrid.vertices[i+new_niv*j].set(vertices[single_index(i,njv-1,j)]); }
		}
		break;
	    case Face.south:
		foreach (i; 0 .. new_niv) {
		    foreach (j; 0 .. new_njv) { bgrid.vertices[i+new_niv*j].set(vertices[single_index(i,0,j)]); }
		}
		break;
	    case Face.west:
		foreach (i; 0 .. new_niv) {
		    foreach (j; 0 .. new_njv) { bgrid.vertices[i+new_niv*j].set(vertices[single_index(0,i,j)]); }
		}
		break;
	    case Face.east:
		foreach (i; 0 .. new_niv) {
		    foreach (j; 0 .. new_njv) { bgrid.vertices[i+new_niv*j].set(vertices[single_index(niv-1,i,j)]); }
		}
		break;
	    case Face.top:
		foreach (i; 0 .. new_niv) {
		    foreach (j; 0 .. new_njv) { bgrid.vertices[i+new_niv*j].set(vertices[single_index(i,j,nkv-1)]); }
		}
		break;
	    case Face.bottom:
		foreach (i; 0 .. niv) {
		    foreach (j; 0 .. njv) { bgrid.vertices[i+new_niv*j].set(vertices[single_index(i,j,0)]); }
		}
		break;
	    } // end switch boundary_indx
	} else {
	    throw new Exception("Extraction from 2D grid not implemented yet.");
	}
	return bgrid;
    } // end get_boundary_grid()

    override size_t[] get_list_of_boundary_cells(size_t boundary_indx)
    // Prepares list of cells indicies that match the boundary grid selected by
    // the function get_boundary_grid().  See above.
    {
	size_t new_nic, new_njc, new_nkc;
	size_t[] cellList;
	if (dimensions == 3) {
	    // 1. Use orientation to decide size of new 2D grid.
	    final switch (boundary_indx) {
	    case Face.north:
		new_nic = niv-1; new_njc = nkv-1; new_nkc = 1; break;
	    case Face.south:
		new_nic = niv-1; new_njc = nkv-1; new_nkc = 1; break;
	    case Face.west:
		new_nic = njv-1; new_njc = nkv-1; new_nkc = 1; break;
	    case Face.east:
		new_nic = njv-1; new_njc = nkv-1; new_nkc = 1; break;
	    case Face.top:
		new_nic = niv-1; new_njc = njv-1; new_nkc = 1; break;
	    case Face.bottom:
		new_nic = niv-1; new_njc = njv-1; new_nkc = 1; break;
	    } // end switch boundary_indx
	} else {
	    throw new Exception("Extraction from 2D grid not implemented yet.");
	}
	cellList.length = new_nic * new_njc * new_nkc;
	if (dimensions == 3) {
	    size_t single_cell_index(size_t i, size_t j, size_t k) { return i + (niv-1)*(j + (njv-1)*k); }
	    final switch (boundary_indx) {
	    case Face.north:
		foreach (i; 0 .. new_nic) {
		    foreach (j; 0 .. new_njc) { cellList[i+new_nic*j] = single_cell_index(i,njv-2,j); }
		}
		break;
	    case Face.south:
		foreach (i; 0 .. new_nic) {
		    foreach (j; 0 .. new_njc) { cellList[i+new_nic*j] = single_cell_index(i,0,j); }
		}
		break;
	    case Face.west:
		foreach (i; 0 .. new_nic) {
		    foreach (j; 0 .. new_njc) { cellList[i+new_nic*j] = single_cell_index(0,i,j); }
		}
		break;
	    case Face.east:
		foreach (i; 0 .. new_nic) {
		    foreach (j; 0 .. new_njc) { cellList[i+new_nic*j] = single_cell_index(niv-2,i,j); }
		}
		break;
	    case Face.top:
		foreach (i; 0 .. new_nic) {
		    foreach (j; 0 .. new_njc) { cellList[i+new_nic*j] = single_cell_index(i,j,nkv-2); }
		}
		break;
	    case Face.bottom:
		foreach (i; 0 .. new_nic) {
		    foreach (j; 0 .. new_njc) { cellList[i+new_nic*j] = single_cell_index(i,j,0); }
		}
		break;
	    } // end switch boundary_indx
	} else {
	    throw new Exception("Extraction from 2D grid not implemented yet.");
	}
	return cellList;
    } // end get_list_of_boundary_cells()
    
    void make_grid_from_path(const Path pth,
			     const(UnivariateFunction) clusterf)
    {
	// First, set up clustered parameter values.
        double[] r = clusterf.distribute_parameter_values(niv);
	// Now, work through the mesh, one point at a time,
        // and create the actual vertex coordinates in Cartesian space.
        size_t k = 0; size_t j = 0;
	foreach (i; 0 .. niv) {
	    Vector3 p = pth(r[i]);
	    this[i,j,k].set(p);
	}
    } // end make_grid_from_path()

    void make_grid_from_surface(const ParametricSurface surf,
				const(UnivariateFunction)[] clusterf)
    {
	// First, set up clustered parameter values along each edge.
        double[] rNorth = clusterf[0].distribute_parameter_values(niv);
        double[] sEast = clusterf[1].distribute_parameter_values(njv);
        double[] rSouth = clusterf[2].distribute_parameter_values(niv);
        double[] sWest = clusterf[3].distribute_parameter_values(njv);
	// Now, work through the mesh, one point at a time,
        // blending the stretched parameter values
        // and creating the actual vertex coordinates in Cartesian space.
        size_t k = 0;
        foreach (j; 0 .. njv) {
            double s = to!double(j) / (njv - 1);
	    foreach (i; 0 .. niv) {
                double r = to!double(i) / (niv - 1);
                double sdash = (1.0-r) * sWest[j] + r * sEast[j]; 
                double rdash = (1.0-s) * rSouth[i] + s * rNorth[i];
                Vector3 p = surf(rdash, sdash);
                this[i,j,k].set(p);
	    }
	}
    } // end make_grid_from_surface()

    void make_grid_from_volume(const ParametricVolume pvolume,
			       const(UnivariateFunction)[] clusterf)
    // Given a parametric volume, create the grid via TFI.
    //
    // The clustering information always comes from the 12 edges.
    {
	// First, set up clustered parameter values along each edge.
        double[] r01 = clusterf[0].distribute_parameter_values(niv);
        double[] s12 = clusterf[1].distribute_parameter_values(njv);
        double[] r32 = clusterf[2].distribute_parameter_values(niv);
        double[] s03 = clusterf[3].distribute_parameter_values(njv);
	//
	double[] r45 = clusterf[4].distribute_parameter_values(niv);
	double[] s56 = clusterf[5].distribute_parameter_values(njv);
	double[] r76 = clusterf[6].distribute_parameter_values(niv);
	double[] s47 = clusterf[7].distribute_parameter_values(njv);
	//
	double[] t04 = clusterf[8].distribute_parameter_values(nkv);
	double[] t15 = clusterf[9].distribute_parameter_values(nkv);
	double[] t26 = clusterf[10].distribute_parameter_values(nkv);
	double[] t37 = clusterf[11].distribute_parameter_values(nkv);
	//
	// Now, work through the mesh, one point at a time,
        // blending the stretched parameter values
        // and creating the actual vertex coordinates in Cartesian space.
        foreach (k; 0 .. nkv) {
	    double t = to!double(k) / (nkv - 1);
	    foreach (j; 0 .. njv) {
		double s = to!double(j) / (njv - 1);
		foreach (i; 0 .. niv) {
		    double r = to!double(i) / (niv - 1);
                    double tdash = (1.0-r)*(1.0-s)*t04[k] + r*s*t26[k] + 
			(1.0-s)*r*t15[k] + s*(1.0-r)*t37[k];
                    double sdash = (1.0-t)*(1.0-r)*s03[j] + t*r*s56[j] + 
			(1.0-t)*r*s12[j] + t*(1-r)*s47[j];
                    double rdash = (1.0-s)*(1.0-t)*r01[i] + s*t*r76[i] + 
			(1.0-s)*t*r45[i] + s*(1.0-t)*r32[i];
		    Vector3 p = pvolume(rdash, sdash, tdash);
		    this[i,j,k].set(p);
		} // i
	    } // j
	} // k
    } // end make_grid_from_volume()

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
	niv = to!int(tokens[0]);
	njv = to!int(tokens[1]);
	nkv = to!int(tokens[2]);
	if (nkv == 1) {
	    if (njv == 1) {
		ncells = niv-1;
	    } else {
		ncells = (niv-1)*(njv-1);
	    }
	} else {
	    ncells = (niv-1)*(njv-1)*(nkv-1);
	}
	nvertices = niv*njv*nkv;
	vertices.length = nvertices;
	vtx_id.length = nvertices;
	// Standard order of vertices.
	size_t ivtx = 0;
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    tokens = f.readln().strip().split();
		    try {
			this[i,j,k].set(to!double(tokens[0]),
					to!double(tokens[1]),
					to!double(tokens[2]));
		    } catch (Exception e) {
			throw new Error(text("Failed to read grid file at " ~
					     "i=", i, " j=", j, " k=", k,
					     "tokens=", tokens, "exception=", e));
		    }
		    vtx_id[single_index(i,j,k)] = ivtx;
		    ivtx++;
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end read_grid_from_text_file()

    override void read_from_gzip_file(string fileName, double scale=1.0)
    // scale = unit length in metres
    {
	auto byLine = new GzipByLine(fileName);
	auto line = byLine.front; byLine.popFront();
	string format_version;
	formattedRead(line, "structured_grid %s", &format_version);
	if (format_version != "1.0") {
	    throw new Error("StructuredGrid.read_from_gzip_file(): " ~
			    "format version found: " ~ format_version); 
	}
	line = byLine.front; byLine.popFront();
	formattedRead(line, "label: %s", &label);
	line = byLine.front; byLine.popFront();
	formattedRead(line, "dimensions: %d", &dimensions);
	line = byLine.front; byLine.popFront();
	formattedRead(line, "niv: %d", &niv);
	line = byLine.front; byLine.popFront();
	formattedRead(line, "njv: %d", &njv);
	line = byLine.front; byLine.popFront();
	formattedRead(line, "nkv: %d", &nkv);
	if (nkv == 1) {
	    if (njv == 1) {
		ncells = niv-1;
	    } else {
		ncells = (niv-1)*(njv-1);
	    }
	} else {
	    ncells = (niv-1)*(njv-1)*(nkv-1);
	}
	nvertices = niv*njv*nkv;
	vertices.length = nvertices;
	vtx_id.length = nvertices;
	// Standard order of vertices.
	size_t ivtx = 0;
	double x, y, z;
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    line = byLine.front; byLine.popFront();
		    // Note that the line starts with whitespace.
		    formattedRead(line, " %g %g %g", &x, &y, &z);
		    this[i,j,k].set(scale*x, scale*y, scale*z);
		    vtx_id[single_index(i,j,k)] = ivtx;
		    ivtx++;
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end read_grid_from_gzip_file()

    override void write_to_gzip_file(string fileName)
    // This function essentially defines the Eilmer4 native format.
    {
	auto f = new GzipOut(fileName);
	auto writer = appender!string();
	formattedWrite(writer, "structured_grid 1.0\n");
	formattedWrite(writer, "label: %s\n", label);
	formattedWrite(writer, "dimensions: %d\n", dimensions);
	formattedWrite(writer, "niv: %d\n", niv);
	formattedWrite(writer, "njv: %d\n", njv);
	formattedWrite(writer, "nkv: %d\n", nkv);
	f.compress(writer.data);
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    writer = appender!string();
		    formattedWrite(writer, "%.18e %.18e %.18e\n", 
				   this[i,j,k].x, this[i,j,k].y, this[i,j,k].z);
		    f.compress(writer.data);
		}
	    }
	}
	f.finish();
    } // end write_grid_to_gzip_file()

    override void write_to_vtk_file(string fileName)
    {
	auto f = File(fileName, "w");
	f.writeln("# vtk DataFile Version 2.0");
	f.writeln(label);
	f.writeln("ASCII");
	f.writeln("");
	f.writeln("DATASET STRUCTURED_GRID");
	f.writefln("DIMENSIONS %d %d %d", niv, njv, nkv);
	f.writefln("POINTS %d float", (niv * njv * nkv));
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    f.writefln("%.18e %.18e %.18e", 
			       this[i,j,k].x, this[i,j,k].y, this[i,j,k].z);
		}
	    }
	}
    } // end write_to_vtk_file()

    /**
     * joinGrid is used to join a supplied grid with
     * the current grid.
     * 
     * Parameters:
     *   gridToJoin   : the supplied grid to be joined with "this" (parent grid)
     *   joinLocation :  is w.r.t "this" grid.
     *                   eg. joinLocation == "north" then that means
     *                   that the gridToJoin is added at the north
     *                   boundary of "this" grid.
     *
     * Notes:
     * + Not all join combinations are possible.
     *   I've only implemented those that are of
     *   immediate use to me. (RJG, 2016-01-22)
     *
     * + Some joins can be achieved by switching
     *   which grid we treat as the parent.
     *   For example, if we want to join the west
     *   of block A to the east of block B, we might
     *   want: A.joinGrid(B, "west") but that option isn't
     *   available. Instead, treat block B as the parent:
     *      B.joinGrid(A, "east").
     *   Some creative thinking like this cuts down on a
     *   lot of implementation code.
     */
    void joinGrid(StructuredGrid gridToJoin, string joinLocation)
    {
	string[] allowedJoins = ["east", "imax",
				 "north", "jmax"];
	if ( find(allowedJoins, joinLocation).empty ) {
	    string errMsg = "Error in StructuredGrid.joinGrid.\n";
	    errMsg ~= "The specified joinLocation = " ~ joinLocation ~ " is not supported.\n";
	    throw new Error(errMsg);
	}
	    
	// Begin by testing if this join is possible
	if (dimensions != 2) {
	    throw new Error("StructuredGrid.joinGrid only implemented for 2D grids.");
	}
	if ( (joinLocation == "east") || (joinLocation == "imax") ) {
	    if ( njv != gridToJoin.njv ) {
		string errMsg = "Error in StructureGrid.joinGrid.\n";
		errMsg ~= "The number of vertices in the j-direction do not match when attempting an east-west join.\n";
		errMsg ~= format("The parent grid has njv= %d\n", njv);
		errMsg ~= format("The grid to be joined has njv= %d\n", gridToJoin.njv);
		throw new Error(errMsg);
	    }
	    if ( nkv != gridToJoin.nkv ) {
		string errMsg = "Error in StructureGrid.joinGrid.\n";
		errMsg ~= "The number of vertices in the k-direction do not match when attempting an east-west join.\n";
		errMsg ~= format("The parent grid has nkv= %d\n", nkv);
		errMsg ~= format("The grid to be joined has nkv= %d\n", gridToJoin.nkv);
		throw new Error(errMsg);
	    }
	}
	if ( (joinLocation == "north") || (joinLocation == "jmax") ) {
	    if ( niv != gridToJoin.niv ) {
		string errMsg = "Error in StructureGrid.joinGrid.\n";
		errMsg ~= "The number of vertices in the i-direction do not match when attempting a north-south join.\n";
		errMsg ~= format("The parent grid has niv= %d\n", niv);
		errMsg ~= format("The grid to be joined has niv= %d\n", gridToJoin.niv);
		throw new Error(errMsg);
	    }
	    if ( nkv != gridToJoin.nkv ) {
		string errMsg = "Error in StructureGrid.joinGrid.\n";
		errMsg ~= "The number of vertices in the k-direction do not match when attempting a north-south join.\n";
		errMsg ~= format("The parent grid has nkv= %d\n", nkv);
		errMsg ~= format("The grid to be joined has nkv= %d\n", gridToJoin.nkv);
		throw new Error(errMsg);
	    }
	}
	// Next we test that the vertices of the joined grids physically coincide (to within some tolerance)
	if ( (joinLocation == "east") || (joinLocation == "imax") ) {
	    foreach ( j; 0 .. njv ) {
		if ( !approxEqualVectors(*(this[niv-1,j]), *(gridToJoin[0,j])) ) {
		    string errMsg = "Error in StructuredGrid.joinGrid.\n";
		    errMsg ~= "At least one of vertices in the join do not coincide.";
		    errMsg ~= "Parent grid vertex: " ~ (*this[niv-1,j]).toString() ~ "\n";
		    errMsg ~= "Join grid vertex: " ~ (*gridToJoin[0,j]).toString() ~ "\n";
		    throw new Error(errMsg);
		}
	    }
	}
	if ( (joinLocation == "north") || (joinLocation == "jmax") ) {
	    foreach ( i; 0 .. niv ) {
		if ( !approxEqualVectors(*(this[i,njv-1]), *(gridToJoin[i,0])) ) {
		    string errMsg = "Error in StructuredGrid.joinGrid.\n";
		    errMsg ~= "At least one of vertices in the join do not coincide.";
		    errMsg ~= "Parent grid vertex: " ~ (*this[i,njv-1]).toString() ~ "\n";
		    errMsg ~= "Join grid vertex: " ~ (*gridToJoin[i,0]).toString() ~ "\n";
		    throw new Error(errMsg);
		}
	    }
	}
	// If the join appears valid, then we can resize the storage
	auto orig_niv = niv;
	auto orig_njv = njv;
	if ( (joinLocation == "east") || (joinLocation == "imax") ) {
	    // -1 because we don't duplicate the coincident vertices at the join
	    niv += gridToJoin.niv - 1;
	}
	if ( (joinLocation == "north") || (joinLocation == "jmax") ) {
	    // -1 because we don't duplicate the coincident vertices at the join
	    njv += gridToJoin.njv - 1; 
	}

	switch (dimensions) {
	case 1: ncells = niv-1; break;
	case 2: ncells = (niv-1)*(njv-1); break;
	case 3: ncells = (niv-1)*(njv-1)*(nkv-1); break;
	default: assert(0);
	}
	nvertices = niv*njv*nkv;
	vertices.length = nvertices;
	vtx_id.length = nvertices;
	// Now we need to add the new vertices.
	if ( (joinLocation == "east") || (joinLocation == "imax") ) {
	    foreach ( j; 0 .. gridToJoin.njv ) {
		foreach ( i; 1 .. gridToJoin.niv ) {
		    *(this[(i-1)+orig_niv,j]) = *(gridToJoin[i,j]);
		}
	    }
	}
	if ( (joinLocation == "north") || (joinLocation == "jmax") ) {
	    foreach ( j; 1 .. gridToJoin.njv ) {
		foreach ( i; 0 .. gridToJoin.niv ) {
		    *(this[i,(j-1)+orig_njv]) = *(gridToJoin[i,j]);
		}
	    }
	}
    } // end joinGrid

} // end class StructuredGrid

//-----------------------------------------------------------------
// Helper functions

// Set very small quantities to zero, exactly.
//
// This is intended primarily to avoid the bad behaviour of VTK
// when it is reading Float32 values that are *too* small.
// We have also come across unreasonably-small float values 
// in the context of reading GridPro files.
double uflowz(double q, double tiny=1.0e-30)
{
    return (fabs(q) > tiny) ? q: 0.0;
}

// Locate the line containing target and return the tokens on that line.
// Legacy-format VTK lines use spaces as delimiters.
string[] read_VTK_header_line(string target, File f)
{
    bool found = false; 
    string[] tokens;
    while (!found) {
        auto line = f.readln();
        if (line.length == 0) break; // presume end of file
	line = strip(line);
        if (indexOf(line, target, CaseSensitive.no) > -1) {
            tokens = split(line);
            found = true; break;
	}
    }
    if (!found) { 
        throw new Error(text("Did not find ", target, " while reading VTK grid file"));
    }
    return tokens;
} // end locate_VTK_header_line()

//-----------------------------------------------------------------

unittest {
    auto p00 = Vector3(0.0, 0.1);
    auto p10 = Vector3(1.0, 0.1);
    auto p11 = Vector3(1.0, 1.1);
    auto p01 = Vector3(0.0, 1.1);
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    auto cf = [new LinearFunction(), new LinearFunction(), 
	       new LinearFunction(), new LinearFunction()];
    auto my_grid = new StructuredGrid(my_patch, 11, 21, cf);
    assert(approxEqualVectors(*my_grid[5,5], Vector3(0.5, 0.35, 0.0)),
			      "StructuredGrid sample point");
    auto my_subgrid = my_grid.subgrid(4, 3, 4, 5);
    assert(approxEqualVectors(*my_subgrid[1,1], Vector3(0.5, 0.35, 0.0)),
			      "subgrid sample point");
}

//-----------------------------------------------------------------

StructuredGrid[] import_gridpro_grid(string fileName, double scale=1.0)
/+
Reads a complete Gridpro grid file, returns a list of StructuredGrids.

A complete Gridpro grid file contains multiple blocks. This function
will read through all blocks and store them as StructuredGrid objects.
These are returned by the function. Gridpro builds grids in the same
dimensions as the supplied geometry. Care should be taken with 
Gridpro grids built from CAD geometries which may typically be
in millimetres. In this case, the required 'scale' would be 0.001
to convert to metres for use in Eilmer.
    
:param fileName: name of Gridpro grid file
:param scale: a scale to convert supplied coordinates to metres
:returns: list of StructuredGrid object(s)

.. Author: Rowan J. Gollan
.. Date: 16-Aug-2012
.. Date: 2014-02-24 ported to D, PJ
+/
{
    auto f = File(fileName, "r");
    StructuredGrid[] grids;
    while (true) {
        auto line = f.readln().strip();
        if (line.length == 0) break;
        if (line[0] == '#') continue;
        auto tks = line.split();
        if (tks.length == 0) break; // Presumably reached end of file
        auto niv = to!int(tks[0]);
	auto njv = to!int(tks[1]);
	auto nkv = to!int(tks[2]);
	writeln("niv=", niv, " njv=", njv, " nkv=", nkv);
        auto mygrid = new StructuredGrid(niv, njv, nkv);
	foreach (i; 0 .. niv) {
            foreach (j; 0 .. njv) {
                foreach (k; 0 .. nkv) {
                    tks = f.readln().strip().split();
                    mygrid[i,j,k].set(uflowz(scale*to!double(tks[0])),
				      uflowz(scale*to!double(tks[1])),
				      uflowz(scale*to!double(tks[2])));
		}
	    }
	}
	grids ~= mygrid;
    } // end while
    return grids;
} // end import_gridpro_grid()
