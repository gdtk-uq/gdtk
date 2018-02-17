// mapped_cell_copy.d

module bc.ghost_cell_effect.mapped_cell_copy;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.algorithm;

import geom;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
import fvcore;
import fvinterface;
import fvcell;
import fluidblock;
import sfluidblock;
import gas;
import bc;

struct BlockAndCellId {
    size_t blkId;
    size_t cellId;

    this(size_t bid, size_t cid)
    {
	blkId = bid;
	cellId = cid;
    }
}

class GhostCellMappedCellCopy : GhostCellEffect {
public:
    // For each ghost cell associated with the current boundary,
    FVCell[] ghost_cells;
    // we will have a corresponding "mapped cell", also known as "source cell"
    // from which we will copy the flow conditions.
    // In the shared-memory flavour of the code, it is easy to get a direct
    // reference to each such mapped cell and store that for easy access.
    FVCell[] mapped_cells;
    // In the MPI-parallel code, we do not have such direct access and so
    // we store the integral ids of the source cell and block and send messages
    // to the source blocks to get the relevant geometry and flow data.
    // The mapped_file has one line per mapped cell and that is identified
    // by the faceTag of the boundary interface.
    BlockAndCellId[string][] mapped_cells_list;
    // size_t[] mapped_cell_ids;
    // size_t[] mapped_cell_blkids;
    //
    // Parameters for the calculation of the mapped-cell location.
    bool cell_mapping_from_file;
    string mapped_cells_filename;
    bool transform_position;
    Vector3 c0 = Vector3(0.0, 0.0, 0.0); // default origin
    Vector3 n = Vector3(0.0, 0.0, 1.0); // z-axis
    double alpha = 0.0; // rotation angle (radians) about specified axis vector
    Vector3 delta = Vector3(0.0, 0.0, 0.0); // default zero translation
    bool list_mapped_cells;
    // Parameters for the optional rotation of copied vector data.
    bool reorient_vector_quantities;
    double[] Rmatrix;

    this(int id, int boundary,
         bool cell_mapping_from_file,
         string mapped_cells_filename,
         bool transform_pos,
         ref const(Vector3) c0, ref const(Vector3) n, double alpha,
         ref const(Vector3) delta,
         bool list_mapped_cells,
         bool reorient_vector_quantities,
         ref const(double[]) Rmatrix)
    {
        super(id, boundary, "MappedCellCopy");
        this.cell_mapping_from_file = cell_mapping_from_file;
        this.mapped_cells_filename = mapped_cells_filename;
        this.transform_position = transform_pos;
        this.c0 = c0;
        this.n = n; this.n.normalize();
        this.alpha = alpha;
        this.delta = delta;
        this.list_mapped_cells = list_mapped_cells;
        this.reorient_vector_quantities = reorient_vector_quantities;
        this.Rmatrix = Rmatrix.dup();
    }

    override string toString() const
    { 
        string str = "MappedCellCopy(" ~
            "cell_mapping_from_file=" ~ to!string(cell_mapping_from_file) ~
            ", mapped_cells_filename=" ~ to!string(mapped_cells_filename) ~
            ", transform_position=" ~ to!string(transform_position) ~
            ", c0=" ~ to!string(c0) ~ 
            ", n=" ~ to!string(n) ~ 
            ", alpha=" ~ to!string(alpha) ~
            ", delta=" ~ to!string(delta) ~
            ", list_mapped_cells=" ~ to!string(list_mapped_cells) ~
            ", reorient_vector_quantities=" ~ to!string(reorient_vector_quantities) ~
            ", Rmatrix=[";
        foreach(i, v; Rmatrix) {
            str ~= to!string(v);
            str ~= (i < Rmatrix.length-1) ? ", " : "]";
        }
        str ~= ")";
        return str;
    }

    void set_up_cell_mapping()
    {
        if (cell_mapping_from_file) {
            set_up_cell_mapping_from_file();
        } else {
            set_up_cell_mapping_via_search();
	}
        foreach (i, mygc; ghost_cells) {
            mygc.copy_values_from(mapped_cells[i], CopyDataOption.grid);
	    if (list_mapped_cells) {
		writeln("    ghost-cell-pos=", to!string(mygc.pos[0]), 
			" mapped-cell-pos=", to!string(mapped_cells[i].pos[0]));
	    }
	}
    } // end set_up_cell_mapping()

    void set_up_cell_mapping_from_file()
    {
        // First, read the entire mapped_cells file.
	// The single mapped_cell file contains the indices mapped cells
	// for all boundary faces in all blocks.
	//
	// They are in sections labelled by the block id.
	// Each boundary face is identified by its "faceTag"
	// which is a string composed of the vertex indices, in ascending order.
	//
	// For the shared memory code, we only need the section for the block
	// associated with the current boundary.
	// For the MPI-parallel code, we need the mappings for all blocks,
	// so that we know what requests for data to expect from other blocks.
        //
	mapped_cells_list.length = GlobalConfig.nFluidBlocks;
        //
        if (!exists(mapped_cells_filename)) {
	    string msg = format("mapped_cells file %s does not exist.", mapped_cells_filename);
            throw new FlowSolverException(msg);
        }
        auto f = File(mapped_cells_filename, "r");
        string getHeaderContent(string target)
        {
	    // Helper function to proceed through file, line-by-line,
	    // looking for a particular header line.
	    // Returns the content from the header line and leaves the file
	    // at the next line to be read, presumably with expected data.
            while (!f.eof) {
                auto line = f.readln().strip();
                if (canFind(line, target)) {
                    auto tokens = line.split("=");
                    return tokens[1].strip();
                }
            } // end while
            return ""; // didn't find the target
        }
	foreach (bid; 0 .. GlobalConfig.nFluidBlocks) {
	    string txt = getHeaderContent(format("NMappedCells in BLOCK[%d]", bid));
	    if (!txt.length) {
		string msg = format("Did not find mapped cells section for block id=%d.", bid);
		throw new FlowSolverException(msg);
	    }
	    size_t nfaces  = to!size_t(txt);
	    foreach(i; 0 .. nfaces) {
		auto lineContent = f.readln().strip();
		auto tokens = lineContent.split();
		string faceTag = tokens[0];
		size_t src_blk_id = to!size_t(tokens[1]);
		size_t src_cell_id = to!size_t(tokens[2]);
		mapped_cells_list[bid][faceTag] = BlockAndCellId(src_blk_id, src_cell_id);
	    }
	} // end foreach bid
        //
        // For the shared-memory code, get references to the mapped-cells for the current block,
        // but only if the source block for each mapped cell is in localFluidBlocks.
        //
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid: 
            BoundaryCondition bc = blk.bc[which_boundary];
            foreach (i, face; bc.faces) {
                ghost_cells ~= (bc.outsigns[i] == 1) ? face.right_cell : face.left_cell;
                size_t[] my_vtx_list; foreach(vtx; face.vtx) { my_vtx_list ~= vtx.id; }
                string faceTag =  makeFaceTag(my_vtx_list);
                auto src_blk_id = mapped_cells_list[blk.id][faceTag].blkId;
                auto src_cell_id = mapped_cells_list[blk.id][faceTag].cellId;
                if (!find(GlobalConfig.localBlockIds, src_blk_id).empty) {
                    mapped_cells ~= globalFluidBlocks[src_blk_id].cells[src_cell_id];
                } else {
		    auto msg = format("block id %d is not in localFluidBlocks", src_blk_id);
                    throw new FlowSolverException(msg);
                }
            } // end foreach face
            break;
        case Grid_t.structured_grid:
            throw new Error("mapped cells from file not implemented for structured grids");
        }
    } // end set_up_cell_mapping_from_file()
    
    void set_up_cell_mapping_via_search()
    {
        // Stage-2 construction for this boundary condition,
        // for the situation when we haven't been told where to find our mapped cells.
        //
        // Needs to be called after the cell geometries have been computed,
        // because the search sifts through the cells in blocks
        // that happen to be in the local process.
        //
        // The search does not extend to cells in blocks in other MPI tasks.
        // If a search for the enclosing cell fails in the MPI context,
        // we will throw an exception rather than continuing the search
        // for the nearest cell.
        //
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid: 
            BoundaryCondition bc = blk.bc[which_boundary];
            foreach (i, face; bc.faces) {
                ghost_cells ~= (bc.outsigns[i] == 1) ? face.right_cell : face.left_cell;
            }
            break;
        case Grid_t.structured_grid:
            size_t i, j, k;
            auto blk = cast(SFluidBlock) this.blk;
            assert(blk !is null, "Oops, this should be an SFluidBlock object.");
            final switch (which_boundary) {
            case Face.north:
                j = blk.jmax;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        ghost_cells ~= blk.get_cell(i,j+1,k);
                        ghost_cells ~= blk.get_cell(i,j+2,k);
                    } // end i loop
                } // for k
                break;
            case Face.east:
                i = blk.imax;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i+1,j,k);
                        ghost_cells ~= blk.get_cell(i+2,j,k);
                    } // end j loop
                } // for k
                break;
            case Face.south:
                j = blk.jmin;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        ghost_cells ~= blk.get_cell(i,j-1,k);
                        ghost_cells ~= blk.get_cell(i,j-2,k);
                    } // end i loop
                } // for k
                break;
            case Face.west:
                i = blk.imin;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i-1,j,k);
                        ghost_cells ~= blk.get_cell(i-2,j,k);
                    } // end j loop
                } // for k
                break;
            case Face.top:
                k = blk.kmax;
                for (i = blk.imin; i <= blk.imax; ++i) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i,j,k+1);
                        ghost_cells ~= blk.get_cell(i,j,k+2);
                    } // end j loop
                } // for i
                break;
            case Face.bottom:
                k = blk.kmin;
                for (i = blk.imin; i <= blk.imax; ++i) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i,j,k-1);
                        ghost_cells ~= blk.get_cell(i,j,k-2);
                    } // end j loop
                } // for i
                break;
            } // end switch
        } // end switch blk.grid_type
        // Now that we have a collection of the local ghost cells,
        // locate the corresponding active cell so that we can later
        // copy that cell's flow state.
        if (list_mapped_cells) {
            writefln("Mapped cells for block[%d] boundary[%d]:", blk.id, which_boundary);
        }
        foreach (mygc; ghost_cells) {
            Vector3 ghostpos = mygc.pos[0];
            Vector3 mypos = ghostpos;
            if (transform_position) {
                Vector3 c1 = c0 + dot(n, (ghostpos - c0)) * n;
                Vector3 t1 = (ghostpos - c1);
                t1.normalize();
                Vector3 t2 = cross(n, t1);
                mypos = c1 + cos(alpha) * t1 + sin(alpha) * t2;
                mypos += delta;
            }
            // Because we need to access all of the gas blocks in the following search,
            // we have to run this set_up_cell_mapping function from a serial loop.
            // In parallel code, threads other than the main thread get uninitialized
            // versions of the localFluidBlocks array.
            //
            // First, attempt to find the enclosing cell at the specified position.
            bool found = false;
            foreach (ib, blk; localFluidBlocks) {
                found = false;
                size_t indx = 0;
                blk.find_enclosing_cell(mypos, indx, found);
                if (found) {
                    mapped_cells ~= blk.cells[indx];
                    break;
                }
            }
            version (mpi_parallel) {
                if (!found && GlobalConfig.in_mpi_context) {
                    string msg = "MappedCellCopy: search for mapped cell did not find an enclosing cell\n";
                    msg ~= "  at position " ~ to!string(mypos) ~ "\n";
                    msg ~= "  This may be because the appropriate cell is not in localFluidBlocks array.\n";
                    throw new FlowSolverException(msg);
                }
            }
            if (!found) {
                // Fall back to nearest cell search.
                FVCell closest_cell = localFluidBlocks[0].cells[0];
                Vector3 cellpos = closest_cell.pos[0];
                Vector3 dp = cellpos - mypos;
                double min_distance = abs(dp);
                foreach (blk; localFluidBlocks) {
                    foreach (cell; blk.cells) {
                        dp = cell.pos[0] - mypos;
                        double distance = abs(dp);
                        if (distance < min_distance) {
                            closest_cell = cell;
                            min_distance = distance;
                        }
                    }
                }
                mapped_cells ~= closest_cell;
            }
            if (list_mapped_cells) {
                writeln("    ghost-cell-pos=", to!string(mygc.pos[0]), 
                        " mapped-cell-pos=", to!string(mapped_cells[$-1].pos[0]));
            }
        } // end foreach mygc
    } // end set_up_cell_mapping_via_search()

    ref FVCell get_mapped_cell(size_t i)
    {
        if (i < mapped_cells.length) {
            return mapped_cells[i];
        } else {
            throw new FlowSolverException(format("Reference to requested mapped-cell[%d] is not available.", i));
        }
    }
    
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        foreach (i; 0 .. ghost_cells.length) {
            ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
            ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
            // [TODO] PJ 2018-01-14 If unstructured blocks ever get used in
            // the block-marching process, we will need a call to encode_conserved
            // at this point.  See the GhostCellFullFaceCopy class.
        }
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        foreach (i; 0 .. ghost_cells.length) {
            ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
            ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
        }
    }
} // end class GhostCellMappedCellCopy
