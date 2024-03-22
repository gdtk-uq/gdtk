/**
 * sgrid.d -- structured-grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-02-22 First code ported from e3_grid.py
 */

module geom.grid.sgrid;

import std.algorithm;
import std.string;
import std.array;
import std.conv;
import std.stdio;
import std.format;
import std.math;
import std.range : iota;
import gzip;
import ntypes.complex;
import nm.number;
import nm.nelmin;

import geom;

//-----------------------------------------------------------------

class StructuredGrid : Grid {
public:
    size_t[] vtx_id; // used to hold the single-index id for each vertex.
    // So that the structured-grid can look like an unstructured grid.
    string[] bcTags; // A user-specified name for each boundary.

    // Blank grid, ready for import of data.
    this(size_t niv, size_t njv, size_t nkv=1, string label="")
    {
        // Infer dimensions from numbers of vertices in each direction.
        // The boundary tags are initially empty strings and may be set later.
        int dim = 3; if (nkv == 1) { dim = (njv == 1) ? 1 : 2; }
        super(Grid_t.structured_grid, dim, label);
        this.niv = niv; this.njv = njv; this.nkv = nkv;
        switch (dim) {
        case 1:
            ncells = niv-1;
            bcTags = ["", ""];
            break;
        case 2:
            ncells = (niv-1)*(njv-1);
            bcTags = ["", "", "", ""];
            break;
        case 3:
            ncells = (niv-1)*(njv-1)*(nkv-1);
            bcTags = ["", "", "", "", "", ""];
            break;
        default:
            throw new GeometryException("invalid number of dimensions");
        }
        nvertices = niv*njv*nkv;
        vertices.length = nvertices;
        vtx_id.length = nvertices;
        // Standard order of vertices is defined in grid.d.
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    size_t ivtx = single_index(i,j,k);
                    vtx_id[ivtx] = ivtx;
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
         const(UnivariateFunction)[] clusterf,
         string interpolation="simple",
         string label="")
    {
        this(niv, njv, 1, label);
        // Any unspecified clustering functions default to the linear identity.
        while (clusterf.length < 4) clusterf ~= new LinearFunction(0.0, 1.0);
        switch (interpolation) {
        case "simple":
            make_grid_from_surface(surf, clusterf);
            break;
        case "linear":
        case "bilinear":
            make_grid_from_surface_bilinear(surf, clusterf);
            break;
        default:
            throw new Error(text("Unknown parameter interpolation: ", interpolation));
        }
    }

    // 2D grid, with redistribution of parameters, built on Parametric surface.
    this(const ParametricSurface surf, size_t niv, size_t njv,
         const(UnivariateFunction)[] clusterf,
         const(double[4][4]) r_grid, const(double[4][4]) s_grid,
         string label="")
    {
        this(niv, njv, 1, label);
        // Any unspecified clustering functions default to the linear identity.
        while (clusterf.length < 4) { clusterf ~= new LinearFunction(0.0, 1.0); }
        make_grid_from_surface(surf, clusterf, r_grid, s_grid);
    }

    // 3D grid, built on ParametricVolume.
    this(const ParametricVolume pvolume, size_t niv, size_t njv, size_t nkv,
         const(UnivariateFunction)[] clusterf,
         string interpolation="simple",
         string label="")
    {
        this(niv, njv, nkv, label);
        // Any unspecified clustering functions default to the linear identity.
        while (clusterf.length < 12) { clusterf ~= new LinearFunction(0.0, 1.0); }
        switch (interpolation) {
        case "simple":
        case "tfi": // keep original name
            make_grid_from_volume(pvolume, clusterf);
            break;
        case "linear":
        case "trilinear":
            make_grid_from_volume_trilinear(pvolume, clusterf);
            break;
        default:
            throw new Error(text("Unknown parameter interpolation: ", interpolation));
        }
    }

    // Imported grid.
    this(string fileName, string fmt, string label="")
    {
        this(1, 1, 1, label); // these settings will be reset on actually reading the data file
        switch (fmt) {
        case "text":
            read_from_text_file(fileName, false);
            break;
        case "gziptext":
            read_from_gzip_file(fileName);
            break;
        case "rawbinary":
            read_from_raw_binary_file(fileName);
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
        foreach (i; 0 .. vertices.length) { vertices[i].set(other.vertices[i]); }
        bcTags = other.bcTags.dup();
    }

    StructuredGrid dup() const
    {
        return new StructuredGrid(this);
    }

    string set_tag(size_t boundary_indx, string new_tag)
    // Sets an individual boundary tag.
    {
        if (boundary_indx < bcTags.length) {
            bcTags[boundary_indx] = new_tag;
            return new_tag;
        } else {
            return "";
        }
    }

    int set_tags(string[string] tags)
    // Sets the boundary tags given the dictionary of tags.
    // Returns the number of tags successfully set.
    {
        string* p;
        int count = 0;
        p = "west" in tags;
        if (p !is null) { bcTags[Face.west] = *p; count++; }
        p = "east" in tags;
        if (p !is null) { bcTags[Face.east] = *p; count++; }
        p = "south" in tags;
        if (p !is null && dimensions >= 2) { bcTags[Face.south] = *p; count++; }
        p = "north" in tags;
        if (p !is null && dimensions >= 2) { bcTags[Face.north] = *p; count++; }
        p = "bottom" in tags;
        if (p !is null && dimensions == 3) { bcTags[Face.bottom] = *p; count++; }
        p = "top" in tags;
        if (p !is null && dimensions == 3) { bcTags[Face.top] = *p; count++; }
        return count;
    }

    override Vector3* opIndex(size_t i, size_t j, size_t k=0)
    in {
        assert (i < niv, text("index i=", i, " is invalid, niv=", niv));
        assert (j < njv, text("index j=", j, " is invalid, njv=", njv));
        assert (k < nkv, text("index k=", k, " is invalid, nkv=", nkv));
    }
    do {
        return &(vertices[single_index(i,j,k)]);
    }

    @nogc
    override Vector3* opIndex(size_t indx)
    in { assert (indx < niv*njv*nkv, "index indx is invalid"); }
    do {
        return &(vertices[indx]);
    }

    override string toString()
    {
        string repr = "StructuredGrid(";
        repr ~= format("dimensions=%d, niv=%d, njv=%d, nkv=%d, ncells=%d",
                       dimensions, niv, njv, nkv, ncells);
        repr ~= format(", vertices=[%s ... %s]", vertices[0], vertices[$-1]);
        repr ~= format(", bcTags=%s", bcTags);
        repr ~= ")";
        return repr;
    }

    @nogc
    override size_t number_of_vertices_for_cell(size_t i)
    {
        switch (dimensions) {
        case 1: return 2;
        case 2: return 4;
        case 3: return 8;
        default: return 0; // oops
        }
    }

    @nogc
    override int vtk_element_type_for_cell(size_t i)
    {
        switch (dimensions) {
        case 1: return VTKElement.line;
        case 2: return VTKElement.quad;
        case 3: return VTKElement.hexahedron;
        default: return 0; // oops
        }
    }

    @nogc
    override int get_cell_type(size_t i)
    {
        return -1; // for structured grid.
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
        default: throw new GeometryException("invalid number of dimensions");
        }
    }
    do {
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
        default: throw new GeometryException("invalid number of dimensions");
        }
    } // end get_vtx_id_list_for_cell()

    @nogc
    override void copy_vtx_id_list_for_cell(ref size_t[8] vtx_list_copy, ref size_t nvtx,
                                            size_t i, size_t j, size_t k=0) const
    {
        switch (dimensions) {
        case 1:
            vtx_list_copy[0] = vtx_id[single_index(i,j,k)];
            vtx_list_copy[1] = vtx_id[single_index(i+1,j,k)];
            nvtx = 2;
            break;
        case 2:
            vtx_list_copy[0] = vtx_id[single_index(i,j,k)];
            vtx_list_copy[1] = vtx_id[single_index(i+1,j,k)];
            vtx_list_copy[2] = vtx_id[single_index(i+1,j+1,k)];
            vtx_list_copy[3] = vtx_id[single_index(i,j+1,k)];
            nvtx = 4;
            break;
        case 3:
            vtx_list_copy[0] = vtx_id[single_index(i,j,k)];
            vtx_list_copy[1] = vtx_id[single_index(i+1,j,k)];
            vtx_list_copy[2] = vtx_id[single_index(i+1,j+1,k)];
            vtx_list_copy[3] = vtx_id[single_index(i,j+1,k)];
            vtx_list_copy[4] = vtx_id[single_index(i,j,k+1)];
            vtx_list_copy[5] = vtx_id[single_index(i+1,j,k+1)];
            vtx_list_copy[6] = vtx_id[single_index(i+1,j+1,k+1)];
            vtx_list_copy[7] = vtx_id[single_index(i,j+1,k+1)];
            nvtx = 8;
            break;
        default: throw new GeometryException("invalid number of dimensions");
        }
    } // end copy_vtx_id_list_for_cell()

    override size_t[] get_vtx_id_list_for_cell(size_t indx) const
    {
        size_t nic, njc, nkc;
        switch (dimensions) {
        case 1: nic = niv-1; njc = 1; nkc = 1; break;
        case 2: nic = niv-1; njc = njv-1; nkc = 1; break;
        case 3: nic = niv-1; njc = njv-1; nkc = njv-1; break;
        default: throw new GeometryException("invalid number of dimensions");
        }
        size_t k = indx / (nic*njc);
        indx -= k * (nic * njc);
        size_t j = indx / nic;
        indx -= j * nic;
        return get_vtx_id_list_for_cell(indx, j, k);
    } // end get_vtx_id_list_for_cell()

    @nogc
    override void copy_vtx_id_list_for_cell(ref size_t[8] vtx_list_copy, ref size_t nvtx,
                                            size_t indx) const
    {
        size_t nic, njc, nkc;
        switch (dimensions) {
        case 1: nic = niv-1; njc = 1; nkc = 1; break;
        case 2: nic = niv-1; njc = njv-1; nkc = 1; break;
        case 3: nic = niv-1; njc = njv-1; nkc = nkv-1; break;
        default: throw new GeometryException("invalid number of dimensions");
        }
        if (nic == 0 || njc == 0 || nkc == 0) {
            throw new GeometryException("Zero cells in at least one direction.");
        }
        size_t k = indx / (nic*njc);
        indx -= k * (nic * njc);
        size_t j = indx / nic;
        indx -= j * nic;
        copy_vtx_id_list_for_cell(vtx_list_copy, nvtx, indx, j, k);
    } // end copy_vtx_id_list_for_cell()

    StructuredGrid subgrid(size_t i0, size_t ni,
                           size_t j0=0, size_t nj=1,
                           size_t k0=0, size_t nk=1) const
    // Partition on vertex indices.
    {
        switch (dimensions) {
        case 1:
            if (ni < 2 || nj != 1 || nk != 1) {
                throw new Error(text("Sgrid.subgrid invalid ni=", ni, " nj=", nj, " nk=", nk,
                                     " for dimensions==1"));
            }
            break;
        case 2:
            if (ni < 2 || nj < 2 || nk != 1) {
                throw new Error(text("Sgrid.subgrid invalid ni=", ni, " nj=", nj, " nk=", nk,
                                     " for dimensions==2"));
            }
            break;
        case 3:
            if (ni < 2 || nj < 2 || nk < 2) {
                throw new Error(text("Sgrid.subgrid invalid ni=", ni, " nj=", nj, " nk=", nk,
                                     " for dimensions==3"));
            }
            break;
        default: throw new GeometryException("invalid number of dimensions");
        }
        if (i0+ni > niv) {
            throw new Error(text("Sgrid.subgrid overrun i0=",i0,", ni=",ni,", niv=",niv));
        }
        if (j0+nj > njv) {
            throw new Error(text("Sgrid.subgrid overrun j0=",j0,", nj=",nj,", njv=",njv));
        }
        if (k0+nk > nkv) {
            throw new Error(text("Sgrid.subgrid overrun k0=",k0,", nk=",nk,", nkv=",nkv));
        }
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
            throw new GeometryException("Extraction from 2D grid not implemented yet.");
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
            throw new GeometryException("Extraction from 2D grid not implemented yet.");
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
            final switch (boundary_indx) {
            case Face.north:
                new_nic = niv-1; new_njc = 1; new_nkc = 1; break;
            case Face.south:
                new_nic = niv-1; new_njc = 1; new_nkc = 1; break;
            case Face.west:
                new_nic = njv-1; new_njc = 1; new_nkc = 1; break;
            case Face.east:
                new_nic = njv-1; new_njc = 1; new_nkc = 1; break;
            }
        }
        cellList.length = new_nic * new_njc * new_nkc;
        size_t single_cell_index(size_t i, size_t j, size_t k) { return i + (niv-1)*(j + (njv-1)*k); }
        if (dimensions == 3) {
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
            final switch (boundary_indx) {
            case Face.north:
                foreach (i; 0 .. new_nic) cellList[i] = single_cell_index(i,njv-2,0);
                break;
            case Face.south:
                foreach (i; 0 .. new_nic) cellList[i] = single_cell_index(i,0,0);
                break;
            case Face.west:
                foreach (i; 0 .. new_nic) cellList[i] = single_cell_index(0,i,0);
                break;
            case Face.east:
                foreach (i; 0 .. new_nic) cellList[i] = single_cell_index(niv-2,i,0);
                break;
            }
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

    // Grid generation with simple interpolation of r,s parameters.
    void make_grid_from_surface(const ParametricSurface surf,
                                const(UnivariateFunction)[] clusterf)
    {
        // First, set up clustered parameter values along each edge.
        double[] sWest = clusterf[Face.west].distribute_parameter_values(njv);
        double[] sEast = clusterf[Face.east].distribute_parameter_values(njv);
        double[] rSouth = clusterf[Face.south].distribute_parameter_values(niv);
        double[] rNorth = clusterf[Face.north].distribute_parameter_values(niv);
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

    // Grid generation with bilinear interpolation of r,s parameters.
    void make_grid_from_surface_bilinear(const ParametricSurface surf,
                                         const(UnivariateFunction)[] clusterf)
    {
        // First, set up clustered parameter values along each edge.
        double[] sWest = clusterf[Face.west].distribute_parameter_values(njv);
        double[] sEast = clusterf[Face.east].distribute_parameter_values(njv);
        double[] rSouth = clusterf[Face.south].distribute_parameter_values(niv);
        double[] rNorth = clusterf[Face.north].distribute_parameter_values(niv);
        // Now, work through the mesh, one point at a time,
        // linearly interpolating the stretched parameter values from edge values
        // and creating the actual vertex coordinates in Cartesian space.
        // See PJ's workbook page 71,72 2021-11-11.
        size_t k = 0;
        foreach (j; 0 .. njv) {
            double s03 = sWest[j];
            double s12 = sEast[j];
            foreach (i; 0 .. niv) {
                double r01 = rSouth[i];
                double r32 = rNorth[i];
                double eta = (r01 + (r32-r01)*s03)/(1.0+(r32-r01)*(s03-s12));
                double sdash = (1.0-eta)*s03 + eta*s12;
                Vector3 p = surf(eta, sdash);
                this[i,j,k].set(p);
            }
        }
    } // end make_grid_from_surface_bilinear()

    // This flavour of the function has internal redistribution of the r,s parameters.
    void make_grid_from_surface(const ParametricSurface surf,
                                const(UnivariateFunction)[] clusterf,
                                const(double[4][4]) r_grid,
                                const(double[4][4]) s_grid)
    {
        // First, set up clustered parameter values along each edge.
        double[] sWest = clusterf[Face.west].distribute_parameter_values(njv);
        double[] sEast = clusterf[Face.east].distribute_parameter_values(njv);
        double[] rSouth = clusterf[Face.south].distribute_parameter_values(niv);
        double[] rNorth = clusterf[Face.north].distribute_parameter_values(niv);
        // Second, set up redistribution function to r,s.
        void remap(double r, double s, out double r_star, out double s_star)
        {
            // Tensor-product, cubic-Bezier interpolation of the gridded parameter values.
            double[4] rr; double[4] sr;
            foreach (j; 0 .. 4) {
                rr[j] = (1.0-r)^^3 * r_grid[0][j] + 3.0*r*(1.0-r)^^2 * r_grid[1][j] +
                    3.0*(1.0-r)*r*r * r_grid[2][j] + r^^3 * r_grid[3][j];
                sr[j] = (1.0-r)^^3 * s_grid[0][j] + 3.0*r*(1.0-r)^^2 * s_grid[1][j] +
                    3.0*(1.0-r)*r*r * s_grid[2][j] + r^^3 * s_grid[3][j];
            }
            r_star = (1.0-s)^^3 * rr[0] + 3.0*s*(1.0-s)^^2 * rr[1] + 3.0*(1.0-s)*s*s * rr[2] + s^^3 * rr[3];
            s_star = (1.0-s)^^3 * sr[0] + 3.0*s*(1.0-s)^^2 * sr[1] + 3.0*(1.0-s)*s*s * sr[2] + s^^3 * sr[3];
        }
        // Now, work through the mesh, one point at a time,
        // blending the stretched parameter values then remapping to rstar,sstar
        // and creating the actual vertex coordinates in Cartesian space.
        size_t k = 0;
        foreach (j; 0 .. njv) {
            double s = to!double(j) / (njv - 1);
            foreach (i; 0 .. niv) {
                double r = to!double(i) / (niv - 1);
                double sdash = (1.0-r)*sWest[j] + r*sEast[j];
                double rdash = (1.0-s)*rSouth[i] + s*rNorth[i];
                double rstar; double sstar;
                remap(rdash, sdash, rstar, sstar);
                Vector3 p = surf(rstar, sstar);
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
        // blending the stretched parameter values with a simple area-weighted interpolation
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

    void make_grid_from_volume_trilinear(const ParametricVolume pvolume,
                                         const(UnivariateFunction)[] clusterf)
    // Given a parametric volume, create the grid via trilinear interpolation.
    // Jens Kunze 2021-June.
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
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    double rb = ((r32[i]-r01[i])*s03[j] + r01[i]) /
                        (1.0 - (r32[i]-r01[i])*(s12[j]-s03[j]));
                    double rt =  ((r76[i]-r45[i])*s47[j] + r45[i]) /
                        (1.0 - (r76[i]-r45[i])*(s56[j]-s47[j]));
                    double tw = ((t37[k]-t04[k])*s03[j] + t04[k]) /
                        (1.0 - (t37[k]-t04[k])*(s47[j]-s03[j]));
                    double te =  ((t26[k]-t15[k])*s12[j] + t15[k]) /
                        (1.0 - (t26[k]-t15[k])*(s56[j]-s12[j]));
                    double rdash = ((rt-rb)*tw + rb) / (1.0 - (rt-rb)*(te-tw));
                    double tdash = (te-tw)*rdash + tw;
                    double sb = (s12[j]-s03[j])*rb + s03[j];
                    double st = (s56[j]-s47[j])*rt + s47[j];
                    double sdash = (st-sb)*tdash + sb;
                    Vector3 p = pvolume(rdash, sdash, tdash);
                    this[i,j,k].set(p);
                } // i
            } // j
        } // k
    } // end make_grid_from_volume_trilinear()

    void read_from_text_file(string fileName, bool vtkHeader=true)
    // Used to import the grid data from a third-party program.
    // We assume that the format is VTK-like.
    {
        string[] tokens;
        auto f = File(fileName, "r");
        if (vtkHeader) {
            read_VTK_header_line("vtk", f);
            label = f.readln().strip();
            read_VTK_header_line("ASCII", f);
            read_VTK_header_line("STRUCTURED_GRID", f);
            tokens = read_VTK_header_line("DIMENSIONS", f);
            niv = to!int(tokens[1]);
            njv = to!int(tokens[2]);
            nkv = to!int(tokens[3]);
            read_VTK_header_line("POINTS", f);
        } else {
            tokens = f.readln().strip().split();
            niv = to!int(tokens[0]);
            njv = to!int(tokens[1]);
            nkv = to!int(tokens[2]);
        }
        dimensions = 0; // start with none
        if (niv > 1 && njv > 1 && nkv > 1) {
            dimensions = 3;
        } else {
            if (niv > 1 && njv > 1) {
                dimensions = 2;
            } else {
                if (niv > 1) { dimensions = 1; }
            }
        }
        if (dimensions == 0) {
            throw new GeometryException(format("Invalid number of vertices" ~
                                               " niv=%d, njv=%d, nkv=%d.",
                                               niv, njv, nkv));
        }
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
    // This is the "native" grid file for Eilmer -- formatted text.
    // scale = unit length in metres
    {
        auto byLine = new GzipByLine(fileName);
        auto line = byLine.front; byLine.popFront();
        string format_version;
        formattedRead(line, "structured_grid %s", &format_version);
        if (!canFind(["1.0", "1.1"], format_version)) {
            throw new Error("StructuredGrid.read_from_gzip_file(): " ~
                            "invalid format version found: " ~ format_version);
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
        //
        // Boundary tags are in version 1.1 plus.
        switch (format_version) {
        case "1.0":
            break;
        case "1.1":
            line = byLine.front; byLine.popFront();
            size_t ntags; formattedRead(line, "ntags: %d", &ntags);
            bcTags.length = ntags;
            foreach (i; 0 .. ntags) {
                line = byLine.front; byLine.popFront();
                formattedRead(line, format("tag[%d]: %%s", i), &(bcTags[i]));
            }
            break;
        default:
            throw new Error("StructuredGrid.read_from_gzip_file(): " ~
                            "invalid format version: " ~ format_version);
        }
        byLine.range.f.close();
    } // end read_grid_from_gzip_file()

    override void read_from_raw_binary_file(string fileName, double scale=1.0)
    // This is the "native" grid file for Eilmer -- binary/raw data.
    // scale = unit length in metres
    {
        File f = File(fileName, "rb");
        string expected_header = "structured_grid 1.x";
        char[] found_header = new char[expected_header.length];
        f.rawRead(found_header);
        string format_version;
        formattedRead(found_header, "structured_grid %s", &format_version);
        if (!canFind(["1.0", "1.1"], format_version)) {
            throw new Error("StructuredGrid.read_from_raw_binary_file(): " ~
                            "invalid format version found: " ~ format_version);
        }
        int[1] buf1; f.rawRead(buf1);
        int label_length = buf1[0];
        if (label_length > 0) {
            char[] found_label = new char[label_length];
            f.rawRead(found_label);
            label = to!string(found_label);
        }
        int[4] buf4; f.rawRead(buf4);
        dimensions = buf4[0];
        niv = buf4[1]; njv = buf4[2]; nkv = buf4[3];
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
        double[3] xyz;
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    f.rawRead(xyz);
                    this[i,j,k].set(scale*xyz[0], scale*xyz[1], scale*xyz[2]);
                    vtx_id[single_index(i,j,k)] = ivtx;
                    ivtx++;
                } // foreach i
            } // foreach j
        } // foreach k
        //
        // Boundary tags are in version 1.1 plus.
        switch (format_version) {
        case "1.0":
            break;
        case "1.1":
            f.rawRead(buf1);
            size_t ntags = buf1[0];
            bcTags.length = ntags;
            foreach (i; 0 .. ntags) {
                f.rawRead(buf1);
                int tag_length = buf1[0];
                if (tag_length > 0) {
                    char[] found_tag = new char[tag_length];
                    f.rawRead(found_tag);
                    bcTags[i] = to!string(found_tag);
                }
            }
            break;
        default:
            throw new Error("StructuredGrid.read_from_raw_binary_file(): " ~
                            "invalid format version: " ~ format_version);
        }
        f.close();
    } // end read_grid_from_raw_binary_file()

    override void write_to_gzip_file(string fileName)
    // This function essentially defines the Eilmer4 "native" format for text output.
    {
        string format_version = "1.1";
        auto f = new GzipOut(fileName);
        auto writer = appender!string();
        formattedWrite(writer, "structured_grid %s\n", format_version);
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
                                   this[i,j,k].x.re, this[i,j,k].y.re, this[i,j,k].z.re);
                    f.compress(writer.data);
                }
            }
        }
        //
        // Boundary tags are in version 1.1 plus.
        switch (format_version) {
        case "1.0":
            break;
        case "1.1":
            writer = appender!string();
            formattedWrite(writer, "ntags: %d\n", bcTags.length);
            foreach (i; 0 .. bcTags.length) {
                formattedWrite(writer, "tag[%d]: %s\n", i, bcTags[i]);
            }
            f.compress(writer.data);
            break;
        default:
            throw new Error("StructuredGrid.write_to_gzip_file(): " ~
                            "invalid format version: " ~ format_version);
        }
        f.finish();
    } // end write_grid_to_gzip_file()

    override void write_to_raw_binary_file(string fileName)
    // This function essentially defines the Eilmer4 "native" raw-binary format.
    {
        File f = File(fileName, "wb");
        string format_version = "1.1";
        string fstr = format("structured_grid %s", format_version);
        f.rawWrite(to!(char[])(fstr));
        int[1] buf1; buf1[0] = to!int(label.length); f.rawWrite(buf1);
        if (label.length > 0) { f.rawWrite(to!(char[])(label)); }
        int[4] buf4; buf4[0] = to!int(dimensions);
        buf4[1] = to!int(niv); buf4[2] = to!int(njv); buf4[3] = to!int(nkv);
        f.rawWrite(buf4);
        double[3] xyz;
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    xyz[0] = this[i,j,k].x.re;
                    xyz[1] = this[i,j,k].y.re;
                    xyz[2] = this[i,j,k].z.re;
                    f.rawWrite(xyz);
                }
            }
        }
        //
        // Boundary tags are in version 1.1 plus.
        switch (format_version) {
        case "1.0":
            break;
        case "1.1":
            buf1[0] = to!int(bcTags.length); f.rawWrite(buf1);
            foreach (i; 0 .. bcTags.length) {
                buf1[0] = to!int(bcTags[i].length); f.rawWrite(buf1);
                f.rawWrite(to!(char[])(bcTags[i]));
            }
            break;
        default:
            throw new Error("StructuredGrid.write_to_raw_binary_file(): " ~
                            "invalid format version: " ~ format_version);
        }
        f.close();
    } // end write_grid_to_raw_binary_file()

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
                               this[i,j,k].x.re, this[i,j,k].y.re, this[i,j,k].z.re);
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
        default: throw new GeometryException("invalid number of dimensions");
        }
        nvertices = niv*njv*nkv;
        vertices.length = nvertices;
        vtx_id.length = nvertices;
        // Now we need to add the new vertices.
        size_t oldSingleIndex(size_t i, size_t j)
        {
            return i + orig_niv*j;
        }
        if ( (joinLocation == "east") || (joinLocation == "imax") ) {
            // Reshuffle points in existing grid
            for (auto j = njv-1; j > 0; --j) {
                for (int i = to!int(orig_niv)-1; i >= 0; --i) {
                   *(this[i, j]) = vertices[oldSingleIndex(i, j)];
                }
            }
            // Then add new points
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
        // Make the vtx_id list up-to-date
        size_t ivtx = 0;
        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    vtx_id[single_index(i,j,k)] = ivtx;
                    ivtx++;
                }
            }
        }
    } // end joinGrid

    StructuredGrid makeSlabGrid(Vector3 dz, bool symmetric=true, string label="")
    {
        assert(nkv == 1, "makeSlabGrid expected only 2D grid");
        string newlabel = (label.length > 0) ? label : this.label;
        StructuredGrid newg = new StructuredGrid(niv, njv, 2, newlabel);
        foreach (j; 0 ..njv) {
            foreach (i; 0 .. niv) {
                if (symmetric) {
                    *(newg[i,j,0]) = *(this[i,j,0]) - 0.5*dz;
                    *(newg[i,j,1]) = *(this[i,j,0]) + 0.5*dz;
                } else {
                    *(newg[i,j,0]) = *(this[i,j,0]);
                    *(newg[i,j,1]) = *(this[i,j,0]) + dz;
                }
            }
        }
        return newg;
    } // end makeSlabGrid()

    StructuredGrid makeSlabGrid(double dz, bool symmetric=true, string label="")
    {
        return makeSlabGrid(Vector3(0,0,dz), symmetric, label);
    }

    StructuredGrid makeWedgeGrid(double dtheta, bool symmetric=true, string label="", int nkv_new = 2)
    {
        assert(nkv == 1, "makeWedgeGrid expected only 2D grid");
        string newlabel = (label.length > 0) ? label : this.label;
        StructuredGrid newg = new StructuredGrid(niv, njv, nkv_new, newlabel);
        foreach (j; 0 ..njv) {
            foreach (i; 0 .. niv) {
                Vector3* p0 = this[i,j,0];
                // We want to rotate the point about the x-axis, according to the right-hand rule.
                // Angles are measured from the y-axis, positive as we swing around toward the z-axis.
                // Refer to PJ's workbook page 36, 2017-07-01
                double r = sqrt((p0.y)^^2 + (p0.z)^^2).re;
                double theta0 = atan2(p0.z, p0.y).re;
                double theta_slice = dtheta / (nkv_new - 1);
                if (symmetric) {
                    double theta1 = theta0-0.5*dtheta;
                    foreach (k; 0 .. nkv_new) {
                        newg[i, j, k].set(p0.x.re, r*cos(theta1 + k * theta_slice), r*sin(theta1 + k * theta_slice));
                    }
                } else {
                    *(newg[i,j,0]) = *p0;
                    foreach (k; 1 .. nkv_new) {
                        newg[i, j, k].set(p0.x.re, r*cos(theta0 + k * theta_slice), r*sin(theta0 + k * theta_slice));
                    }
                }
            }
        }
        return newg;
    } // end makeWedgeGrid()

    override void write_to_su2_file(string fileName, double scale=1.0,
                                    bool use_gmsh_order_for_wedges=true)
    {
        // We can write a unstructured grid file but I have no idea why you would want to.
        auto usg = new UnstructuredGrid(this);
        usg.write_to_su2_file(fileName, scale, use_gmsh_order_for_wedges);
    }

    double measure_of_badness()
    // Compute the cell badness measures for a given grid.
    // A grid that has a value of 1.0 is (likely) an ideal, orthogonal mesh.
    {
        double ratio_of_diagonals_avg = 0.0; double ratio_of_diagonals_max = 0.0;
        double nonorthogonality_avg = 0.0; double nonorthogonality_max = 0.0;
        double slenderness_avg = 0.0; double slenderness_max = 0.0;
        double area_min = double.max; double area_max = 0.0;
        double negative_area_penalty = 0.0;
        if (nkv == 1 && njv > 1) {
            // 2D surface mesh.
            // Compute measure on nonorthogonality at every vertex.
            size_t k = 0;
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    Vector3 a = *this[i,j,k];
                    Vector3 b; if (i < niv-1) { b.set(*this[i+1,j,k]); } else { b.set(*this[i-1,j,k]); }
                    Vector3 c; if (j < njv-1) { c.set(*this[i,j+1,k]); } else { b.set(*this[i,j-1,k]); }
                    Vector3 ab = b; ab -= a; ab.normalize();
                    Vector3 ac = c; ac -= a; ac.normalize();
                    double dotprod = fabs(ab.dot(ac).re);
                    nonorthogonality_avg += dotprod;
                    nonorthogonality_max = max(nonorthogonality_max, dotprod);
                }
            }
            size_t nvtx = njv*niv;
            nonorthogonality_avg /= nvtx; // average value per vertex
            foreach (j; 0 .. njv-1) {
                foreach (i; 0 .. niv-1) {
                    // Corners of the cell[i,j]
                    Vector3 p00 = *this[i,j,k]; Vector3 p10 = *this[i+1,j,k];
                    Vector3 p01 = *this[i,j+1,k]; Vector3 p11 = *this[i+1,j+1,k];
                    // Diagonals
                    Vector3 d0 = p00; d0 -= p11;
                    Vector3 d1 = p10; d1 -= p01;
                    double L20 = d0.dot(d0).re + 1.0e-16;
                    double L21 = d1.dot(d1).re + 1.0e-16;
                    double ratio_of_diagonals = (L20 > L21) ? L20/L21 : L21/L20;
                    ratio_of_diagonals_max = max(ratio_of_diagonals_max, ratio_of_diagonals);
                    ratio_of_diagonals_avg += ratio_of_diagonals;
                    // Slenderness
                    Vector3 centroid;
                    number xyplane_area, iLen, jLen, minLen;
                    xyplane_quad_cell_properties(p00, p10, p11, p01, centroid,
                                                 xyplane_area, iLen, jLen, minLen);
                    number perimeter = distance_between(p00, p10) + distance_between(p10, p11) +
                        distance_between(p11, p01) + distance_between(p01, p00);
                    number slenderness = (perimeter^^2)/(16*xyplane_area);
                    slenderness_max = max(slenderness_max, slenderness.re);
                    slenderness_avg += slenderness.re;
                    area_min = min(area_min, xyplane_area.re);
                    area_max = max(area_max, xyplane_area.re);
                    // Negative areas are a bad thing.
                    negative_area_penalty += (xyplane_area < 0.0) ? 10000.0 : 0.0;
                }
            }
            size_t ncells = (njv-1)*(niv-1);
            ratio_of_diagonals_avg /= ncells; // want average value per cell
            slenderness_avg /= ncells;
        } else {
            // 3D volume mesh.
            throw new Error("Not implemented.");
        }
        /+
        writefln("nonorthogonality max=%e avg=%g diagonals max=%e avg=%e"~
                 " slenderness max=%e avg=%e area_max/min=%e negative=%e",
                 nonorthogonality_max, nonorthogonality_avg,
                 ratio_of_diagonals_max, ratio_of_diagonals_avg,
                 slenderness_max, slenderness_avg,
                 (area_max/area_min), negative_area_penalty);
        +/
        return 7.0*nonorthogonality_max + 3.0*ratio_of_diagonals_max + slenderness_max +
            7.0*nonorthogonality_avg + 3.0*ratio_of_diagonals_avg + slenderness_avg +
            (area_max/area_min) + negative_area_penalty;
    } // end measure_of_badness()

    bool determine_rs_grids(const ParametricSurface surf, const(UnivariateFunction)[] clusterf,
                            ref double[4][4] r_grid, ref double[4][4] s_grid)
    {
        double objective(double[] data)
        {
            list_to_rs_grids(data, r_grid, s_grid);
            make_grid_from_surface(surf, clusterf, r_grid, s_grid);
            double obj = measure_of_badness();
            return obj;
        }
        double[] data; data.length = 16;
        rs_grids_to_list(r_grid, s_grid, data);
        double[] ddata; foreach (i; 0 .. data.length) { ddata ~= 0.05; }
        double fdata;
        int nfe, nres;
        bool conv_flag = nm.nelmin.minimize!(objective,double)(data, fdata, nfe, nres, ddata,
                                                               1.0e-3, 3, 800);
        /+
        debug {
            writeln("data = ", data, " fdata = ", fdata);
            writeln("convergence-flag = ", conv_flag);
            writeln("number-of-fn-evaluations = ", nfe);
            writeln("number-of-restarts = ", nres);
        }
        +/
        list_to_rs_grids(data, r_grid, s_grid);
        return conv_flag;
    } // end determine_rs_grids()

    void rotateGrid(double q0, double q1, double q2, double q3)
    {
        // Flynn Hack 2024-03-23
        /+ Rotate grid by applying the quaternion rotation of the form
                 q0 + q1 ihat + q2 jhat + q3 khat
            to all points. Some good info on quaternions and their relation
            to rotation matrices can be found in:
            https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation

        +/

        number norm = q0*q0 + q1*q1 + q2*q2 + q3*q3;

        if ((abs(norm) - 1.0) > 1e-8) {
            throw new Exception("Invalid quaternion provided.
                                Ensure norm is equal to 1.0.");
        }

        number r0, r1, r2, r3, r4, r5, r6, r7, r8;

        r0 = q0*q0 + q1*q1 - q2*q2 - q3*q3;
        r1 = 2*q1*q2 - 2*q0*q3;
        r2 = 2*q1*q3 + 2*q0*q2;
        r3 = 2*q1*q2 + 2*q0*q3;
        r4 = q0*q0 - q1*q1 + q2*q2 - q3*q3;
        r5 = 2*q2*q3 - 2*q0*q1;
        r6 = 2*q1*q3 - 2*q0*q2;
        r7 = 2*q2*q3 + 2*q0*q1;
        r8 = q0*q0 - q1*q1 - q2*q2 + q3*q3; // Implementation checked against
                                            // scipy rotation package.

        const(number[]) Rmat = [r0, r1, r2, r3, r4, r5, r6, r7, r8];

        foreach (k; 0 .. nkv) {
            foreach (j; 0 .. njv) {
                foreach (i; 0 .. niv) {
                    this[i,j,k].apply_matrix_transform(Rmat);
                }
            }
        }
    } // end rotateGrid

} // end class StructuredGrid

//-----------------------------------------------------------------
// Helper functions

void rs_grids_to_list(ref const(double[4][4]) r_grid, ref const(double[4][4]) s_grid,
                      ref double[] data)
// Mapping of redistribution r,s values from individual grids to parameter list.
// See PJ workbook page 63 2020-07-25.
// Should match function below.
{
    assert(data.length >= 16, "data array too short");
    data[0] = s_grid[0][1];
    data[1] = s_grid[0][2];
    data[2] = r_grid[1][0];
    data[3] = r_grid[1][1];
    data[4] = s_grid[1][1];
    data[5] = r_grid[1][2];
    data[6] = s_grid[1][2];
    data[7] = r_grid[1][3];
    data[8] = r_grid[2][0];
    data[9] = r_grid[2][1];
    data[10] = s_grid[2][1];
    data[11] = r_grid[2][2];
    data[12] = s_grid[2][2];
    data[13] = r_grid[2][3];
    data[14] = s_grid[3][1];
    data[15] = s_grid[3][2];
    return;
}

void list_to_rs_grids(ref const(double[]) data,
                      ref double[4][4] r_grid, ref double[4][4] s_grid)
// Mapping of redistribution r,s values from parameter list to individual grids.
// See PJ workbook page 63 2020-07-25.
// Should match function above.
{
    assert(data.length >= 16, "data array too short");
    r_grid[0][0] = 0.0;      s_grid[0][0] = 0.0;
    r_grid[0][1] = 0.0;      s_grid[0][1] = data[0];
    r_grid[0][2] = 0.0;      s_grid[0][2] = data[1];
    r_grid[0][3] = 0.0;      s_grid[0][3] = 1.0;
    r_grid[1][0] = data[2];  s_grid[1][0] = 0.0;
    r_grid[1][1] = data[3];  s_grid[1][1] = data[4];
    r_grid[1][2] = data[5];  s_grid[1][2] = data[6];
    r_grid[1][3] = data[7];  s_grid[1][3] = 1.0;
    r_grid[2][0] = data[8];  s_grid[2][0] = 0.0;
    r_grid[2][1] = data[9];  s_grid[2][1] = data[10];
    r_grid[2][2] = data[11]; s_grid[2][2] = data[12];
    r_grid[2][3] = data[13]; s_grid[2][3] = 1.0;
    r_grid[3][0] = 1.0;      s_grid[3][0] = 0.0;
    r_grid[3][1] = 1.0;      s_grid[3][1] = data[14];
    r_grid[3][2] = 1.0;      s_grid[3][2] = data[15];
    r_grid[3][3] = 1.0;      s_grid[3][3] = 1.0;
    return;
}

// Set very small quantities to zero, exactly.
//
// This is intended primarily to avoid the bad behaviour of VTK
// when it is reading Float32 values that are *too* small.
// We have also come across unreasonably-small float values
// in the context of reading GridPro files.
// Edit: (NNG 2020), additionally we now cap values that are too
// large as well. This is avoids a problem with paraviews ascii
// grid reader that causes memory errors when given values that
// invalid float32 values.
@nogc double uflowz(double q, double tiny=1.0e-30, double huge=1.0e30)
{
    double qsafe = q;
    if (fabs(q) < tiny) qsafe = 0.0;
    if (fabs(q) > huge) qsafe = copysign(huge, q);
    return qsafe;
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
} // end read_VTK_header_line()

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

StructuredGrid[] importPlot3DGrid(string fileName, int dim, double scale=1.0)
{
    if (dim != 2 && dim != 3) {
        string errMsg = "ERROR in importPlot3DGrids: 'dim' must be 2 or 3";
        throw new Error(errMsg);
    }

    auto f = File(fileName, "r");
    // Determine if we have multiblock file, or not.
    // A multiblock file will have a single integer on line
    // to give the number of blocks.
    // A single-block file will start immediately with numbers
    // that specify size of the grid in each dimension.
    auto line = f.readln().strip();
    auto tokens = line.split();
    int nBlocks = (tokens.length == 1) ? to!int(tokens[0]) : 1;
    // Next gather up the ni, nj and (possibly) nk dimensions.
    StructuredGrid[] grids;
    if ((nBlocks == 1) && (tokens.length > 1)) {
        // Already have read the line with grid size in each dimesion.
        auto niv = to!int(tokens[0]);
        auto njv = to!int(tokens[1]);
        auto nkv = (dim == 3) ? to!int(tokens[2]) : 1;
        grids ~= new StructuredGrid(niv, njv, nkv);
    }
    else {
        foreach (iBlk; 0 .. nBlocks) {
            line = f.readln().strip();
            tokens = line.split();
            auto niv = to!int(tokens[0]);
            auto njv = to!int(tokens[1]);
            auto nkv = (dim == 3) ? to!int(tokens[2]) : 1;
            grids ~= new StructuredGrid(niv, njv, nkv);
        }
    }
    // Now loop for all blocks, gathering x then y then z coordinates.
    // We only gather z for dim = 3. Plot3D considers a 2D grid as
    // one in the x-y plane. It is physically 2D. Not "logically" 2D
    // as a surface in 3D might be.
    //
    // Different programs seem to write the coordinate values with variable layout,
    // sometimes one per line, sometimes several per line.
    // Let's gather them all into one big list and then
    // assign them to the individual blocks.
    number[] coords;
    line = f.readln().strip();
    while (line.length > 0) {
        tokens = line.split();
        foreach (item; tokens) { coords ~= to!number(item); }
        line = f.readln().strip();
    }
    size_t ic =0;
    foreach (g; grids) {
        foreach (k; 0 .. g.nkv) {
            foreach (j; 0 .. g.njv) {
                foreach (i; 0 .. g.niv) {
                    g[i,j,k].x = coords[ic]; ic++;
                }
            }
        }
        foreach (k; 0 .. g.nkv) {
            foreach (j; 0 .. g.njv) {
                foreach (i; 0 .. g.niv) {
                    g[i,j,k].y = coords[ic]; ic++;
                }
            }
        }
        if (dim == 3) {
            foreach (k; 0 .. g.nkv) {
                foreach (j; 0 .. g.njv) {
                    foreach (i; 0 .. g.niv) {
                        g[i,j,k].z = coords[ic]; ic++;
                    }
                }
            }
        } // end if dim == 3
        else {
            foreach (k; 0 .. g.nkv) {
                foreach (j; 0 .. g.njv) {
                    foreach (i; 0 .. g.niv) {
                        g[i,j,k].z = 0.0;
                    }
                }
            }
        }
    } // end foreach g
    f.close();
    return grids;
} // end importPlot3DGrid()

void writeGridsAsPlot3D(string fname, StructuredGrid[] grids, int dim)
{
    if ( dim != 2 && dim != 3 ) {
        string errMsg = "ERROR in writeGridsAsPlot3D: 'dim' must be 2 or 3";
        throw new Error(errMsg);
    }
    auto f = File(fname, "w");
    f.writefln(" %d", grids.length); // Always write multi-block format
                                     // even if we only have a single block.
    foreach (g; grids) {
        if ( dim == 2 ) {
            f.writefln(" %d %d", g.niv, g.njv);
        }
        else {
            f.writefln(" %d %d %d", g.niv, g.njv, g.nkv);
        }
    }
    foreach (g; grids) {
        foreach (k; 0 .. g.nkv) {
            foreach (j; 0 .. g.njv) {
                foreach (i; 0 .. g.niv) {
                    f.writefln(" %.18e", g[i,j,k].x);
                }
            }
        }
        foreach (k; 0 .. g.nkv) {
            foreach (j; 0 .. g.njv) {
                foreach (i; 0 .. g.niv) {
                    f.writefln(" %.18e", g[i,j,k].y);
                }
            }
        }
        if ( dim == 3 ) {
            foreach (k; 0 .. g.nkv) {
                foreach (j; 0 .. g.njv) {
                    foreach (i; 0 .. g.niv) {
                        f.writefln(" %.18e", g[i,j,k].z);
                    }
                }
            }
        }
    }
    f.close();
} // end writeGridsAsPlot3D()

StructuredGrid rotate_gridpro_blocks(StructuredGrid grid,
                                     string rotateSouthToThis="South",
                                     string rotateWestToThis="West")
{
    StructuredGrid new_grid;

    if (rotateSouthToThis == "East" && rotateWestToThis == "South") {
    new_grid = new StructuredGrid(grid.njv, grid.niv);
        foreach (j; 0 .. grid.niv) {
            foreach (i; 0 .. grid.njv) {
                new_grid[i,j,0].set(grid[(grid.niv-1)-j,i].x,
                                    grid[(grid.niv-1)-j,i].y,
                                    grid[(grid.niv-1)-j,i].z);
            }
        }
    return new_grid;
    }

    else if (rotateSouthToThis == "North" && rotateWestToThis == "East") {
    new_grid = new StructuredGrid(grid.niv, grid.njv);
        foreach (j; 0 .. grid.njv) {
            foreach (i; 0 .. grid.niv) {
                new_grid[i,j,0].set(grid[(grid.niv-1)-i,(grid.njv-1)-j].x,
                                    grid[(grid.niv-1)-i,(grid.njv-1)-j].y,
                                    grid[(grid.niv-1)-i,(grid.njv-1)-j].z);
            }
        }
    return new_grid;
    }

    else if (rotateSouthToThis == "West" && rotateWestToThis == "North") {
    new_grid = new StructuredGrid(grid.njv, grid.niv);
        foreach (j; 0 .. grid.niv) {
            foreach (i; 0 .. grid.njv) {
                new_grid[i,j,0].set(grid[j,(grid.njv-1)-i].x,
                                    grid[j,(grid.njv-1)-i].y,
                                    grid[j,(grid.njv-1)-i].z);
            }
        }
    return new_grid;
    }
    return grid;
} // end rotate_gridpro_blocks()

StructuredGrid grid_faceswap(StructuredGrid grid,
                             bool swapNorthToSouth=false,
                             bool swapEastToWest=false)
{
    StructuredGrid new_grid;

    if (swapNorthToSouth == true) {
    new_grid = new StructuredGrid(grid.niv, grid.njv);
        foreach (j; 0 .. grid.njv) {
            foreach (i; 0 .. grid.niv) {
                new_grid[i,j,0].set(grid[i,(grid.njv-1)-j].x,
                                    grid[i,(grid.njv-1)-j].y,
                                    grid[i,(grid.njv-1)-j].z);
            }
        }
    }

    if (swapEastToWest == true) {
    new_grid = new StructuredGrid(grid.niv, grid.njv);
        foreach (j; 0 .. grid.njv) {
            foreach (i; 0 .. grid.niv) {
                new_grid[i,j,0].set(grid[(grid.niv-1)-i,j].x,
                                    grid[(grid.niv-1)-i,j].y,
                                    grid[(grid.niv-1)-i,j].z);
            }
        }
    }
    return new_grid;
} // end grid_faceswap()


//-----------------------------------------------------------------

version(sgrid_test) {
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
        my_grid.set_tag(Face.north, "north-face");
        my_grid.set_tags(["south":"south-face-here-sir", "west":"west-face-here-sir!"]);
        string repr = to!string(my_grid);
        // writeln("my_grid=", repr);
        assert(repr.canFind("north-face"), failedUnitTest());
        assert(repr.canFind("south-face-here-sir"), failedUnitTest());
        assert(repr.canFind("west-face-here-sir!"), failedUnitTest());
        assert(approxEqualVectors(*my_grid[5,5], Vector3(0.5, 0.35, 0.0)),
               failedUnitTest());
        auto my_subgrid = my_grid.subgrid(4, 3, 4, 5);
        assert(approxEqualVectors(*my_subgrid[1,1], Vector3(0.5, 0.35, 0.0)),
               failedUnitTest());

        assert(uflowz(1.0)==1.0, failedUnitTest());
        assert(uflowz(1.0e-32)==0.0, failedUnitTest());
        assert(uflowz(1.0e-29)==1.0e-29, failedUnitTest());
        assert(uflowz(1.0e+29)==1.0e+29, failedUnitTest());
        assert(uflowz(1.0e+32)==1.0e+30, failedUnitTest());

        assert(uflowz(-1.0)==-1.0, failedUnitTest());
        assert(uflowz(-1.0e-32)==0.0, failedUnitTest());
        assert(uflowz(-1.0e-29)==-1.0e-29, failedUnitTest());
        assert(uflowz(-1.0e+29)==-1.0e+29, failedUnitTest());
        assert(uflowz(-1.0e+32)==-1.0e+30, failedUnitTest());
        return 0;
    }
} // end sgrid_test
