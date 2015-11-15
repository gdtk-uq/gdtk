/**
 * sgrid.d -- structured-grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-02-22 First code ported from e3_grid.py
 */

module sgrid;

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

//-----------------------------------------------------------------

enum Grid_t {structured_grid, unstructured_grid}

string gridTypeName(Grid_t gt)
{
    final switch (gt) {
    case Grid_t.structured_grid: return "structured_grid";
    case Grid_t.unstructured_grid: return "unstructured_grid";
    }
}

Grid_t gridTypeFromName(string name)
{
    switch (name) {
    case "structured_grid": return Grid_t.structured_grid;
    case "unstructured_grid": return Grid_t.unstructured_grid;
    default: throw new Error("Unknown type of grid: " ~ name);
    }
}

class Grid {
    Grid_t grid_type;
    int dimensions; // 2 or 3
    string label;
    size_t ncells;
    size_t nvertices;
    Vector3[] vertices;
    size_t[] vtx_id;
    
    this(Grid_t grid_type, int dimensions, string label="")
    {
	this.grid_type = grid_type;
	this.dimensions = dimensions;
	this.label = label;
    }

    abstract ref Vector3 opIndex(size_t i, size_t j, size_t k=0);
    abstract ref Vector3 opIndex(size_t indx);
    abstract size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const; 
    abstract size_t[] get_vtx_id_list_for_cell(size_t indx) const;
    abstract void read_from_gzip_file(string fileName);
    abstract void write_to_gzip_file(string fileName);
    abstract void write_to_vtk_file(string fileName);
}

//-----------------------------------------------------------------

class StructuredGrid : Grid {
public:
    int niv, njv, nkv;

    size_t single_index(size_t i, size_t j, size_t k=0) const
    in {
	assert (i < niv, text("index i=", i, " is invalid, niv=", niv));
	assert (j < njv, text("index j=", j, " is invalid, njv=", njv));
	assert (k < nkv, text("index k=", k, " is invalid, nkv=", nkv));
    }
    body {
	return i + niv*(j + njv*k);
    }

    // Blank grid, ready for import of data.
    this(int niv, int njv, int nkv=1, string label="")
    {
	int dim = (nkv == 1) ? 2 : 3; // infer dimensions
	super(Grid_t.structured_grid, dim, label);
	this.niv = niv; this.njv = njv; this.nkv = nkv;
	if (dim == 2) {
	    ncells = (niv-1)*(njv-1);
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
		    vtx_id[single_index(i,j,k)] = ivtx;
		    ivtx++;
		}
	    }
	}
    }

    // 2D grid, built on Parametric surface.
    this(const ParametricSurface surf, int niv, int njv,
	 const(UnivariateFunction)[] clusterf, string label="")
    {
	this(niv, njv, 1, label);
	// Any unspecified clustering functions default to the linear identity.
	while (clusterf.length < 4) clusterf ~= new LinearFunction(0.0, 1.0);
	make_grid_from_surface(surf, clusterf);
    }

    // 3D grid, built on ParametricVolume.
    this(const ParametricVolume pvolume, int niv, int njv, int nkv,
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
	this(0, 0, 0, label); // these settings will be reset on actually reading the data file
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
	
    override ref Vector3 opIndex(size_t i, size_t j, size_t k=0)
    in {
	assert (i < niv, text("index i=", i, " is invalid, niv=", niv));
	assert (j < njv, text("index j=", j, " is invalid, njv=", njv));
	assert (k < nkv, text("index k=", k, " is invalid, nkv=", nkv));
    }
    body {
	return vertices[single_index(i,j,k)];
    }

    override ref Vector3 opIndex(size_t indx)
    in {
	assert (indx < niv*njv*nkv,
		text("index indx=", indx, " is invalid, niv*njv*nkv=", niv*njv*nkv));
    }
    body {
	return vertices[indx];
    }

    override size_t[] get_vtx_id_list_for_cell(size_t i, size_t j, size_t k=0) const
    in {
	size_t nic = niv - 1;
	assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
	size_t njc = njv - 1;
	assert (j < njc, text("index j=", j, " is invalid, njc=", njc));
	if (dimensions == 2) {
	    assert (k == 0, text("index k=", k, " is invalid for 2D grid"));
	} else {
	    size_t nkc = nkv - 1;
	    assert (k < nkc, text("index k=", k, " is invalid, nkc=", nkc));
	}
    }
    body {
	if (dimensions == 2) {
	    return [vtx_id[single_index(i,j,k)],
		    vtx_id[single_index(i+1,j,k)],
		    vtx_id[single_index(i+1,j+1,k)],
		    vtx_id[single_index(i,j+1,k)]];
	} else {
	    return [vtx_id[single_index(i,j,k)],
		    vtx_id[single_index(i+1,j,k)], 
		    vtx_id[single_index(i+1,j+1,k)],
		    vtx_id[single_index(i,j+1,k)],
		    vtx_id[single_index(i,j,k+1)],
		    vtx_id[single_index(i+1,j,k+1)], 
		    vtx_id[single_index(i+1,j+1,k+1)],
		    vtx_id[single_index(i,j+1,k+1)]];
	}
    }

    override size_t[] get_vtx_id_list_for_cell(size_t indx) const
    {
	size_t nic = niv - 1;
	size_t njc = njv - 1;
	size_t nkc = 1;
	if (dimensions == 3) {
	    nkc = nkv - 1;
	}
	size_t k = indx / (nic*njc);
	indx -= k * (nic * njc);
	size_t j = indx / nic;
	indx -= j * nic;
	return get_vtx_id_list_for_cell(indx, j, k);
    }

    StructuredGrid subgrid(size_t i0, size_t ni,
			   size_t j0, size_t nj,
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
		    new_grd[i,j,k].refx = vertices[single_index(i0+i,j0+j,k0+k)].x;
		    new_grd[i,j,k].refy = vertices[single_index(i0+i,j0+j,k0+k)].y;
		    new_grd[i,j,k].refz = vertices[single_index(i0+i,j0+j,k0+k)].z;
		}
	    }
	}
	return new_grd;
    } // end subgrid()

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
        int k = 0;
        foreach (j; 0 .. njv) {
            double s = to!double(j) / (njv - 1);
	    foreach (i; 0 .. niv) {
                double r = to!double(i) / (niv - 1);
                double sdash = (1.0-r) * sWest[j] + r * sEast[j]; 
                double rdash = (1.0-s) * rSouth[i] + s * rNorth[i];
                Vector3 p = surf(rdash, sdash);
                this[i,j,k].refx = p.x;
                this[i,j,k].refy = p.y;
                this[i,j,k].refz = p.z;
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
		    this[i,j,k].refx = p.x;
		    this[i,j,k].refy = p.y;
		    this[i,j,k].refz = p.z;
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
	    ncells = (niv-1)*(njv-1);
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
			this[i,j,k].refx = to!double(tokens[0]);
			this[i,j,k].refy = to!double(tokens[1]);
			this[i,j,k].refz = to!double(tokens[2]);
		    } catch (Exception e) {
			throw new Error(text("Failed to read grid file at "
					     "i=", i, " j=", j, " k=", k,
					     "tokens=", tokens, "exception=", e));
		    }
		    vtx_id[single_index(i,j,k)] = ivtx;
		    ivtx++;
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end read_grid_from_text_file()

    override void read_from_gzip_file(string fileName)
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
	    ncells = (niv-1)*(njv-1);
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
		    this[i,j,k].refx = x;
		    this[i,j,k].refy = y;
		    this[i,j,k].refz = z;
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
		    formattedWrite(writer, "%20.16e %20.16e %20.16e\n", 
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
		    f.writefln("%20.16e %20.16e %20.16e", 
			       this[i,j,k].x, this[i,j,k].y, this[i,j,k].z);
		}
	    }
	}
    } // end write_to_vtk_file()

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
    assert(approxEqualVectors(my_grid[5,5], Vector3(0.5, 0.35, 0.0)),
			      "StructuredGrid sample point");
    auto my_subgrid = my_grid.subgrid(4, 3, 4, 5);
    assert(approxEqualVectors(my_subgrid[1,1], Vector3(0.5, 0.35, 0.0)),
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
                    mygrid[i,j,k].refx = uflowz(scale*to!double(tks[0]));
                    mygrid[i,j,k].refy = uflowz(scale*to!double(tks[1]));
                    mygrid[i,j,k].refz = uflowz(scale*to!double(tks[2]));
		}
	    }
	}
	grids ~= mygrid;
    } // end while
    return grids;
} // end import_gridpro_grid()
