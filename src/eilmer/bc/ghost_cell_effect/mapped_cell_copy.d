// mapped_cell_copy.d

module bc.ghost_cell_effect.mapped_cell_copy;

import core.memory;
import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.algorithm;
import nm.complex;
import nm.number;
version(mpi_parallel) {
    import mpi;
}

import geom;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
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
    // Flow data along the boundary is stored in ghost cells.
    FVCell[] ghost_cells;
    size_t[string] ghost_cell_index_from_faceTag;
    // For each ghost-cell associated with the current boundary,
    // we will have a corresponding "mapped cell", also known as "source cell"
    // from which we will copy the flow conditions.
    // In the shared-memory flavour of the code, it is easy to get a direct
    // reference to each such mapped cell and store that for easy access.
    FVCell[] mapped_cells;
    // We may specify which source cell and block from which a particular ghost-cell
    // (a.k.a. destination cell) will copy its flow and geometry data.
    // This mapping information is prepared externally and provided in
    // a single mapped_cells file which has one line per mapped cell.
    // The first item on each line specifies the boundary face associated with
    // with the ghost cell via the faceTag.
    bool cell_mapping_from_file;
    string mapped_cells_filename;
    BlockAndCellId[string][size_t] mapped_cells_list;
    version(mpi_parallel) {
        // In the MPI-parallel code, we do not have such direct access and so
        // we store the integral ids of the source cell and block and send requests
        // to the source blocks to get the relevant geometry and flow data.
        // The particular cells and the order in which they are packed into the
        // data pipes need to be known at the source and destination ends of the pipes.
        // So, we store those cell indices in a matrix of lists with the indices
        // into the matrix being the source and destination block ids.
        size_t[][size_t][size_t] src_cell_ids;
        size_t[][size_t][size_t] ghost_cell_indices;
        //
        size_t n_incoming, n_outgoing;
        size_t[] outgoing_ncells_list, incoming_ncells_list;
        size_t[] outgoing_block_list, incoming_block_list;
        int[] outgoing_rank_list, incoming_rank_list;
        int[] outgoing_geometry_tag_list, incoming_geometry_tag_list;
        MPI_Request[] incoming_geometry_request_list;
        MPI_Status[] incoming_geometry_status_list;
        double[][] outgoing_geometry_buf_list, incoming_geometry_buf_list;
        int[] outgoing_flowstate_tag_list, incoming_flowstate_tag_list;
        MPI_Request[] incoming_flowstate_request_list;
        MPI_Status[] incoming_flowstate_status_list;
        double[][] outgoing_flowstate_buf_list, incoming_flowstate_buf_list;
        int[] outgoing_convective_gradient_tag_list, incoming_convective_gradient_tag_list;
        MPI_Request[] incoming_convective_gradient_request_list;
        MPI_Status[] incoming_convective_gradient_status_list;
        double[][] outgoing_convective_gradient_buf_list, incoming_convective_gradient_buf_list;
        int[] outgoing_viscous_gradient_tag_list, incoming_viscous_gradient_tag_list;
        MPI_Request[] incoming_viscous_gradient_request_list;
        MPI_Status[] incoming_viscous_gradient_status_list;
        double[][] outgoing_viscous_gradient_buf_list, incoming_viscous_gradient_buf_list;
    }
    //
    // Parameters for the calculation of the mapped-cell location.
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

    // not @nogc
    void set_up_cell_mapping()
    {
        if (cell_mapping_from_file) {
            if (!exists(mapped_cells_filename)) {
                string msg = format("mapped_cells file %s does not exist.", mapped_cells_filename);
                throw new FlowSolverException(msg);
            }
            final switch (blk.grid_type) {
            case Grid_t.unstructured_grid:
                // We set up the ghost-cell reference list to have the same order as
                // the list of faces that were stored in the boundary.
                // We will later confirm that the ghost cells appear in the same order
                // in the mapped_cells file.
                BoundaryCondition bc = blk.bc[which_boundary];
                foreach (i, face; bc.faces) {
                    ghost_cells ~= (bc.outsigns[i] == 1) ? face.right_cell : face.left_cell;
                    size_t[] my_vtx_list; foreach(vtx; face.vtx) { my_vtx_list ~= vtx.id; }
                    string faceTag =  makeFaceTag(my_vtx_list);
                    ghost_cell_index_from_faceTag[faceTag] = i;
                    ghost_cells[$-1].is_interior_to_domain = true;
                }
                break;
            case Grid_t.structured_grid:
                throw new Error("cell mapping from file is not implemented for structured grids");
            } // end switch grid_type
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
            //
            // Read the entire mapped_cells file.
            // The single mapped_cell file contains the indices mapped cells
            // for all ghost-cells, for all blocks.
            //
            // They are in sections labelled by the block id.
            // Each boundary face is identified by its "faceTag",
            // which is a string composed of the vertex indices, in ascending order.
            // The order of the ghost-cells is assumed the same as for each
            // grids underlying the FluidBlock.
            //
            // For the shared memory code, we only need the section for the block
            // associated with the current boundary.
            // For the MPI-parallel code, we need the mappings for all blocks,
            // so that we know what requests for data to expect from other blocks.
            //
            string txt = getHeaderContent(format("MappedBlocks in BLOCK[%d]", blk.id));
            if (!txt.length) {
                string msg = format("Did not find mapped blocks section for block id=%d.", blk.id);
                throw new FlowSolverException(msg);
            }
            size_t[] neighbour_block_id_list;
            neighbour_block_id_list ~= blk.id;
            foreach (id; txt.split()) {
                neighbour_block_id_list ~= to!int(id);
            }
            neighbour_block_id_list.sort();
            //
            foreach (dest_blk_id; neighbour_block_id_list) {
                txt = getHeaderContent(format("NMappedCells in BLOCK[%d]", dest_blk_id));
                if (!txt.length) {
                    string msg = format("Did not find mapped cells section for destination block id=%d.",
                                        dest_blk_id);
                    throw new FlowSolverException(msg);
                }
                size_t nfaces  = to!size_t(txt);
                foreach(i; 0 .. nfaces) {
                    auto lineContent = f.readln().strip();
                    auto tokens = lineContent.split();
                    string faceTag = tokens[0];
                    size_t src_blk_id = to!size_t(tokens[1]);
                    size_t src_cell_id = to!size_t(tokens[2]);
                    mapped_cells_list[dest_blk_id][faceTag] = BlockAndCellId(src_blk_id, src_cell_id);
                    version(mpi_parallel) {
                        // These lists will be used to direct data when packing and unpacking
                        // the buffers used to send data between the MPI tasks.
                        src_cell_ids[src_blk_id][dest_blk_id] ~= src_cell_id;
                        ghost_cell_indices[src_blk_id][dest_blk_id] ~= i;
                        // If we are presently reading the section for the current block,
                        // we check that the listed faces are in the same order as the
                        // underlying grid.
                        if (blk.id == dest_blk_id) {
                            if (canFind(ghost_cell_index_from_faceTag.keys(), faceTag)) {
                                if (i != ghost_cell_index_from_faceTag[faceTag]) {
                                    throw new Error(format("Oops, ghost-cell indices do not match: %d %d",
                                                           i, ghost_cell_index_from_faceTag[faceTag]));
                                }
                            } else {
                                foreach (ft; ghost_cell_index_from_faceTag.keys()) {
                                    writefln("ghost_cell_index_from_faceTag[\"%s\"] = %d",
                                             ft, ghost_cell_index_from_faceTag[ft]);
                                }
                                throw new Error(format("Oops, cannot find faceTag=\"%s\" for block id=%d", faceTag, blk.id));
                            }
                        }
                    }
                }
            } // end foreach dest_blk_id
            //
            version(mpi_parallel) {
                //
                // No communication needed just now because all MPI tasks have the full mapping,
                // however, we can prepare buffers and the like for communication of the geometry
                // and flowstate data.
                //
                // Incoming messages will carrying data from other block, to be copied into the
                // ghost cells for the current boundary.
                // N.B. We assume that there is only one such boundary per block.
                incoming_ncells_list.length = 0;
                incoming_block_list.length = 0;
                incoming_rank_list.length = 0;
                incoming_geometry_tag_list.length = 0;
                incoming_flowstate_tag_list.length = 0;
                incoming_convective_gradient_tag_list.length = 0;
                incoming_viscous_gradient_tag_list.length = 0;
                foreach (src_blk_id; neighbour_block_id_list) {
                    if (src_blk_id == blk.id) {continue;}
                    size_t nc = src_cell_ids[src_blk_id][blk.id].length;
                    if (nc > 0) {
                        incoming_ncells_list ~= nc;
                        incoming_block_list ~= src_blk_id;
                        incoming_rank_list ~= GlobalConfig.mpi_rank_for_block[src_blk_id];
                        incoming_geometry_tag_list ~= make_mpi_tag(to!int(src_blk_id), 99, 1);
                        incoming_flowstate_tag_list ~= make_mpi_tag(to!int(src_blk_id), 99, 2);
                        incoming_convective_gradient_tag_list ~= make_mpi_tag(to!int(src_blk_id), 99, 3);
                        incoming_viscous_gradient_tag_list ~= make_mpi_tag(to!int(src_blk_id), 99, 4);
                    }
                }
                n_incoming = incoming_block_list.length;
                incoming_geometry_request_list.length = n_incoming;
                incoming_geometry_status_list.length = n_incoming;
                incoming_geometry_buf_list.length = n_incoming;
                incoming_flowstate_request_list.length = n_incoming;
                incoming_flowstate_status_list.length = n_incoming;
                incoming_flowstate_buf_list.length = n_incoming;
                incoming_convective_gradient_request_list.length = n_incoming;
                incoming_convective_gradient_status_list.length = n_incoming;
                incoming_convective_gradient_buf_list.length = n_incoming;
                incoming_viscous_gradient_request_list.length = n_incoming;
                incoming_viscous_gradient_status_list.length = n_incoming;
                incoming_viscous_gradient_buf_list.length = n_incoming;
                //
                // Outgoing messages will carry data from source cells in the current block,
                // to be copied into ghost cells in another block.
                outgoing_ncells_list.length = 0;
                outgoing_block_list.length = 0;
                outgoing_rank_list.length = 0;
                outgoing_geometry_tag_list.length = 0;
                outgoing_flowstate_tag_list.length = 0;
                outgoing_convective_gradient_tag_list.length = 0;
                outgoing_viscous_gradient_tag_list.length = 0;
                foreach (dest_blk_id; neighbour_block_id_list) {
                    if (dest_blk_id == blk.id) {continue;}
                    size_t nc = src_cell_ids[blk.id][dest_blk_id].length;
                    if (nc > 0) {
                        outgoing_ncells_list ~= nc;
                        outgoing_block_list ~= dest_blk_id;
                        outgoing_rank_list ~= GlobalConfig.mpi_rank_for_block[dest_blk_id];
                        outgoing_geometry_tag_list ~= make_mpi_tag(to!int(blk.id), 99, 1);
                        outgoing_flowstate_tag_list ~= make_mpi_tag(to!int(blk.id), 99, 2);
                        outgoing_convective_gradient_tag_list ~= make_mpi_tag(to!int(blk.id), 99, 3);
                        outgoing_viscous_gradient_tag_list ~= make_mpi_tag(to!int(blk.id), 99, 4);
                    }
                }
                n_outgoing = outgoing_block_list.length;
                outgoing_geometry_buf_list.length = n_outgoing;
                outgoing_flowstate_buf_list.length = n_outgoing;
                outgoing_convective_gradient_buf_list.length = n_outgoing;
                outgoing_viscous_gradient_buf_list.length = n_outgoing;
                //
            } else { // not mpi_parallel
                // For the shared-memory code, get references to the mapped (source) cells
                // that need to be accessed for the current (destination) block.
                final switch (blk.grid_type) {
                case Grid_t.unstructured_grid:
                    BoundaryCondition bc = blk.bc[which_boundary];
                    foreach (i, face; bc.faces) {
                        size_t[] my_vtx_list; foreach(vtx; face.vtx) { my_vtx_list ~= vtx.id; }
                        string faceTag =  makeFaceTag(my_vtx_list);
                        auto src_blk_id = mapped_cells_list[blk.id][faceTag].blkId;
                        auto src_cell_id = mapped_cells_list[blk.id][faceTag].cellId;
                        if (!find(GlobalConfig.localFluidBlockIds, src_blk_id).empty) {
                            auto blk = cast(FluidBlock) globalBlocks[src_blk_id];
                            assert(blk !is null, "Oops, this should be a FluidBlock object.");
                            mapped_cells ~= blk.cells[src_cell_id];
                        } else {
                            auto msg = format("block id %d is not in localFluidBlocks", src_blk_id);
                            throw new FlowSolverException(msg);
                        }
                    } // end foreach face
                    break;
                case Grid_t.structured_grid:
                    throw new Error("cell mapping from file not implemented for structured grids");
                } // end switch grid_type
            } // end not mpi_parallel
        } else { // !cell_mapping_from_file
            set_up_cell_mapping_via_search();
        } // end if !cell_mapping_from_file
    } // end set_up_cell_mapping()

    // not @nogc
    void set_up_cell_mapping_via_search()
    {
        // For the situation when we haven't been given a file to specify
        // where to find our mapped cells.
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
            auto blk = cast(SFluidBlock) this.blk;
            assert(blk !is null, "Oops, this should be an SFluidBlock object.");
            bool nghost3 = (blk.n_ghost_cell_layers == 3);
            final switch (which_boundary) {
            case Face.north:
                size_t j = blk.njc;
                foreach (k; 0 .. blk.nkc) {
                    foreach (i; 0 .. blk.nic) {
                        auto f = blk.get_ifj(i,j,k);
                        foreach (n; 0 .. blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                    }
                }
                break;
            case Face.east:
                size_t i = blk.nic;
                foreach (k; 0 .. blk.nkc) {
                    foreach (j; 0 .. blk.njc) {
                        auto f = blk.get_ifi(i,j,k);
                        foreach (n; 0 .. blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                    }
                }
                break;
            case Face.south:
                size_t j = 0;
                foreach (k; 0 .. blk.nkc) {
                    foreach (i; 0 .. blk.nic) {
                        auto f = blk.get_ifj(i,j,k);
                        foreach (n; 0 .. blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                    }
                }
                break;
            case Face.west:
                size_t i = 0;
                foreach (k; 0 .. blk.nkc) {
                    foreach (j; 0 .. blk.njc) {
                        auto f = blk.get_ifi(i,j,k);
                        foreach (n; 0 .. blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                    }
                }
                break;
            case Face.top:
                size_t k = blk.nkc;
                foreach (i; 0 .. blk.nic) {
                    foreach (j; 0 .. blk.njc) {
                        auto f = blk.get_ifk(i,j,k);
                        foreach (n; 0 .. blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                    }
                }
                break;
            case Face.bottom:
                size_t k = 0;
                foreach (i; 0 .. blk.nic) {
                    foreach (j; 0 .. blk.njc) {
                        auto f = blk.get_ifk(i,j,k);
                        foreach (n; 0 .. blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                    }
                }
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
            mygc.is_interior_to_domain = true; // Ghost cells need to be flagged as interior
            if (transform_position) {
                Vector3 del1 = ghostpos - c0;
                Vector3 c1 = c0 + dot(n, del1) * n;
                Vector3 t1 = ghostpos - c1;
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
                double min_distance = distance_between(cellpos, mypos);
                foreach (blk; localFluidBlocks) {
                    foreach (cell; blk.cells) {
                        double distance = distance_between(cell.pos[0], mypos);
                        if (distance < min_distance) {
                            closest_cell = cell;
                            min_distance = distance;
                        }
                    }
                }
                mapped_cells ~= closest_cell;
            }
        } // end foreach mygc
	// TODO: temporarily removing the GC calls below, they are (oddly) computationally expensive - KD 26/03/2019. 
        //GC.collect();
        //GC.minimize();
    } // end set_up_cell_mapping_via_search()

    @nogc
    ref FVCell get_mapped_cell(size_t i)
    {
        if (i < mapped_cells.length) {
            return mapped_cells[i];
        } else {
            throw new FlowSolverException("Reference to requested mapped-cell is not available.");
        }
    }

    // not @nogc
    void exchange_geometry_phase0()
    {
        version(mpi_parallel) {
            // Prepare to exchange geometry data for the boundary cells.
            foreach (i; 0 .. n_incoming) {
                // To match .copy_values_from(mapped_cells[i], CopyDataOption.grid) as defined in fvcell.d.
                size_t ne = incoming_ncells_list[i] * (blk.myConfig.n_grid_time_levels * 5 + 5);
                version(complex_numbers) { ne *= 2; }
                if (incoming_geometry_buf_list[i].length < ne) { incoming_geometry_buf_list[i].length = ne; }
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the src_blk MPI process.
                MPI_Irecv(incoming_geometry_buf_list[i].ptr, to!int(ne), MPI_DOUBLE, incoming_rank_list[i],
                          incoming_geometry_tag_list[i], MPI_COMM_WORLD, &incoming_geometry_request_list[i]);
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_geometry_phase0()

    // not @nogc
    void exchange_geometry_phase1()
    {
        version(mpi_parallel) {
            foreach (i; 0 .. n_outgoing) {
                // Blocking send of this block's geometry data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                size_t ne = outgoing_ncells_list[i] * (blk.myConfig.n_grid_time_levels * 5 + 5);
                version(complex_numbers) { ne *= 2; }
                if (outgoing_geometry_buf_list[i].length < ne) { outgoing_geometry_buf_list[i].length = ne; }
                auto buf = outgoing_geometry_buf_list[i];
                size_t ii = 0;
                foreach (cid; src_cell_ids[blk.id][outgoing_block_list[i]]) {
                    auto c = blk.cells[cid];
                    foreach (j; 0 .. blk.myConfig.n_grid_time_levels) {
                        buf[ii++] = c.pos[j].x.re; version(complex_numbers) { buf[ii++] = c.pos[j].x.im; }
                        buf[ii++] = c.pos[j].y.re; version(complex_numbers) { buf[ii++] = c.pos[j].y.im; } 
                        buf[ii++] = c.pos[j].z.re; version(complex_numbers) { buf[ii++] = c.pos[j].z.im; }
                        buf[ii++] = c.volume[j].re; version(complex_numbers) { buf[ii++] = c.volume[j].im; }
                        buf[ii++] = c.areaxy[j].re; version(complex_numbers) { buf[ii++] = c.areaxy[j].im; }
                    }
                    buf[ii++] = c.iLength.re; version(complex_numbers) { buf[ii++] = c.iLength.im; }
                    buf[ii++] = c.jLength.re; version(complex_numbers) { buf[ii++] = c.jLength.im; }
                    buf[ii++] = c.kLength.re; version(complex_numbers) { buf[ii++] = c.kLength.im; }
                    buf[ii++] = c.L_min.re; version(complex_numbers) { buf[ii++] = c.L_min.im; }
                    buf[ii++] = c.L_max.re; version(complex_numbers) { buf[ii++] = c.L_max.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                              outgoing_geometry_tag_list[i], MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                             outgoing_geometry_tag_list[i], MPI_COMM_WORLD);
                }
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_geometry_phase1()

    // not @nogc
    void exchange_geometry_phase2()
    {
        version(mpi_parallel) {
            foreach (i; 0 .. n_incoming) {
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_geometry_request_list[i], &incoming_geometry_status_list[i]);
                } else {
                    MPI_Wait(&incoming_geometry_request_list[i], &incoming_geometry_status_list[i]);
                }
                auto buf = incoming_geometry_buf_list[i];
                size_t ii = 0;
                foreach (gi; ghost_cell_indices[incoming_block_list[i]][blk.id]) {
                    auto c = ghost_cells[gi];
                    foreach (j; 0 .. blk.myConfig.n_grid_time_levels) {
                        c.pos[j].x.re = buf[ii++]; version(complex_numbers) { c.pos[j].x.im = buf[ii++]; }
                        c.pos[j].y.re = buf[ii++]; version(complex_numbers) { c.pos[j].y.im = buf[ii++]; }
                        c.pos[j].z.re = buf[ii++]; version(complex_numbers) { c.pos[j].z.im = buf[ii++]; }
                        c.volume[j].re = buf[ii++]; version(complex_numbers) { c.volume[j].im = buf[ii++]; }
                        c.areaxy[j].re = buf[ii++]; version(complex_numbers) { c.areaxy[j].im = buf[ii++]; }
                    }

                    static foreach(li; 0 .. 3) {
                        c.lengths[li].re = buf[ii++]; version(complex_numbers) { c.lengths[li].im = buf[ii++]; }
                    }

                    c.L_min.re = buf[ii++]; version(complex_numbers) { c.L_min.im = buf[ii++]; }
                    c.L_max.re = buf[ii++]; version(complex_numbers) { c.L_max.im = buf[ii++]; }
                }
            }
        } else { // not mpi_parallel
            // For a single process, just access the data directly.
            foreach (i, mygc; ghost_cells) {
                mygc.copy_values_from(mapped_cells[i], CopyDataOption.grid);
            }
        }
    } // end exchange_geometry_phase2()

    @nogc
    size_t flowstate_buffer_entry_size(const LocalConfig myConfig)
    {
        /*
        Compute the amount of space needed for one flowstate in the SEND/RECV buffer 
        Previously, this code was duplicated in two places, and it was easy to 
        get bitten by not changing both at once when adding things to the flowstate variables

        Note: This routine must be kept consistent with the buffer packing in exchange_flowstate
        phases 1 and 2
        @author: Nick N. Gibbons
        */

        size_t nspecies = myConfig.n_species;
        size_t nmodes = myConfig.n_modes;

        size_t nitems = 16;
        nitems += nmodes*3; // for each of T, e and k_t
        nitems += nspecies*2; // for massf and rho_s
        version(MHD) { nitems += 5; }
        version(turbulence) { nitems += myConfig.turb_model.nturb; }

        version(complex_numbers) {
            nitems *= 2;
        }
        return nitems;
    }


    // not @nogc
    void exchange_flowstate_phase0(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            // Prepare to exchange geometry data for the boundary cells.
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            foreach (i; 0 .. n_incoming) {
                // Exchange FlowState data for the boundary cells.
                // To match the function over in flowstate.d
                // void copy_values_from(in FlowState other)
                // and over in gas_state.d
                // @nogc void copy_values_from(ref const(GasState) other)
                //size_t ne = incoming_ncells_list[i] * (nmodes*3 + nspecies + 23);
                size_t fs_size = flowstate_buffer_entry_size(blk.myConfig);
                size_t ne = incoming_ncells_list[i] * fs_size;
                if (incoming_flowstate_buf_list[i].length < ne) { incoming_flowstate_buf_list[i].length = ne; }
                // Post non-blocking receive for flowstate data that we expect to receive later
                // from the src_blk MPI process.
                MPI_Irecv(incoming_flowstate_buf_list[i].ptr, to!int(ne), MPI_DOUBLE, incoming_rank_list[i],
                          incoming_flowstate_tag_list[i], MPI_COMM_WORLD, &incoming_flowstate_request_list[i]);
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_flowstate_phase0()

    // not @nogc
    void exchange_flowstate_phase1(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach (i; 0 .. n_outgoing) {
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                //size_t nitems = 16;
                //version(MHD) { nitems += 5; }
                //version(turbulence) { nitems += 2; }
                size_t fs_size = flowstate_buffer_entry_size(blk.myConfig);
                size_t ne = outgoing_ncells_list[i] * fs_size;
                if (outgoing_flowstate_buf_list[i].length < ne) { outgoing_flowstate_buf_list[i].length = ne; }
                auto buf = outgoing_flowstate_buf_list[i];
                size_t ii = 0;
                foreach (cid; src_cell_ids[blk.id][outgoing_block_list[i]]) {
                    auto c = blk.cells[cid];
                    FlowState* fs = &(c.fs);
                    GasState* gs = &(fs.gas);
                    buf[ii++] = gs.rho.re; version(complex_numbers) { buf[ii++] = gs.rho.im; }
                    buf[ii++] = gs.p.re; version(complex_numbers) { buf[ii++] = gs.p.im; }
                    buf[ii++] = gs.T.re; version(complex_numbers) { buf[ii++] = gs.T.im; }
                    buf[ii++] = gs.u.re; version(complex_numbers) { buf[ii++] = gs.u.im; }
                    buf[ii++] = gs.p_e.re; version(complex_numbers) { buf[ii++] = gs.p_e.im; }
                    buf[ii++] = gs.a.re; version(complex_numbers) { buf[ii++] = gs.a.im; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { buf[ii++] = gs.u_modes[j].re; version(complex_numbers) { buf[ii++] = gs.u_modes[j].im; } }
                        foreach (j; 0 .. nmodes) { buf[ii++] = gs.T_modes[j].re; version(complex_numbers) { buf[ii++] = gs.T_modes[j].im; } }
                    }
                    buf[ii++] = gs.mu.re; version(complex_numbers) { buf[ii++] = gs.mu.im; }
                    buf[ii++] = gs.k.re; version(complex_numbers) { buf[ii++] = gs.k.im; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { buf[ii++] = gs.k_modes[j].re; version(complex_numbers) { buf[ii++] = gs.k_modes[j].im; } }
                    }
                    buf[ii++] = gs.sigma.re; version(complex_numbers) { buf[ii++] = gs.sigma.im; }
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) { buf[ii++] = gs.massf[j].re; version(complex_numbers) { buf[ii++] = gs.massf[j].im; } }
                        foreach (j; 0 .. nspecies) { buf[ii++] = gs.rho_s[j].re; version(complex_numbers) { buf[ii++] = gs.rho_s[j].im; } }
                    }
                    buf[ii++] = gs.quality.re; version(complex_numbers) { buf[ii++] = gs.quality.im; }
                    buf[ii++] = fs.vel.x.re; version(complex_numbers) { buf[ii++] = fs.vel.x.im; }
                    buf[ii++] = fs.vel.y.re; version(complex_numbers) { buf[ii++] = fs.vel.y.im; }
                    buf[ii++] = fs.vel.z.re; version(complex_numbers) { buf[ii++] = fs.vel.z.im; }
                    version(MHD) {
                        buf[ii++] = fs.B.x.re; version(complex_numbers) { buf[ii++] = fs.B.x.im; }
                        buf[ii++] = fs.B.y.re; version(complex_numbers) { buf[ii++] = fs.B.y.im; }
                        buf[ii++] = fs.B.z.re; version(complex_numbers) { buf[ii++] = fs.B.z.im; }
                        buf[ii++] = fs.psi.re; version(complex_numbers) { buf[ii++] = fs.psi.im; }
                        buf[ii++] = fs.divB.re; version(complex_numbers) { buf[ii++] = fs.divB.im; }
                    }
                    version(turbulence) {
                        foreach (j; 0 .. nturb) { buf[ii++] = fs.turb[j].re; version(complex_numbers) { buf[ii++] = fs.turb[j].im; } }
                    }
                    buf[ii++] = fs.mu_t.re; version(complex_numbers) { buf[ii++] = fs.mu_t.im; }
                    buf[ii++] = fs.k_t.re; version(complex_numbers) { buf[ii++] = fs.k_t.im; }
                    buf[ii++] = fs.S.re; version(complex_numbers) { buf[ii++] = fs.S.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                              outgoing_flowstate_tag_list[i], MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                             outgoing_flowstate_tag_list[i], MPI_COMM_WORLD);
                }
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_flowstate_phase1()

    // not @nogc
    void exchange_flowstate_phase2(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach (i; 0 .. n_incoming) {
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_flowstate_request_list[i], &incoming_flowstate_status_list[i]);
                } else {
                    MPI_Wait(&incoming_flowstate_request_list[i], &incoming_flowstate_status_list[i]);
                }
                auto buf = incoming_flowstate_buf_list[i];
                size_t ii = 0;
                foreach (gi; ghost_cell_indices[incoming_block_list[i]][blk.id]) {
                    auto c = ghost_cells[gi];
                    FlowState* fs = &(c.fs);
                    GasState* gs = &(fs.gas);
                    gs.rho.re = buf[ii++]; version(complex_numbers) { gs.rho.im = buf[ii++]; }
                    gs.p.re = buf[ii++]; version(complex_numbers) { gs.p.im = buf[ii++]; }
                    gs.T.re = buf[ii++]; version(complex_numbers) { gs.T.im = buf[ii++]; }
                    gs.u.re = buf[ii++]; version(complex_numbers) { gs.u.im = buf[ii++]; }
                    gs.p_e.re = buf[ii++]; version(complex_numbers) { gs.p_e.im = buf[ii++]; }
                    gs.a.re = buf[ii++]; version(complex_numbers) { gs.a.im = buf[ii++]; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { gs.u_modes[j].re = buf[ii++]; version(complex_numbers) { gs.u_modes[j].im = buf[ii++]; } }
                        foreach (j; 0 .. nmodes) { gs.T_modes[j].re = buf[ii++]; version(complex_numbers) { gs.T_modes[j].im = buf[ii++]; } }
                    }
                    gs.mu.re = buf[ii++]; version(complex_numbers) { gs.mu.im = buf[ii++]; }
                    gs.k.re = buf[ii++]; version(complex_numbers) { gs.k.im = buf[ii++]; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { gs.k_modes[j].re = buf[ii++]; version(complex_numbers) { gs.k_modes[j].im = buf[ii++]; } }
                    }
                    gs.sigma.re = buf[ii++]; version(complex_numbers) { gs.sigma.im = buf[ii++]; }
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) { gs.massf[j].re = buf[ii++]; version(complex_numbers) { gs.massf[j].im = buf[ii++]; } }
                        foreach (j; 0 .. nspecies) { gs.rho_s[j].re = buf[ii++]; version(complex_numbers) { gs.rho_s[j].im = buf[ii++]; } }
                    }
                    gs.quality.re = buf[ii++]; version(complex_numbers) { gs.quality.im = buf[ii++]; }
                    fs.vel.x.re = buf[ii++]; version(complex_numbers) { fs.vel.x.im = buf[ii++]; }
                    fs.vel.y.re = buf[ii++]; version(complex_numbers) { fs.vel.y.im = buf[ii++]; }
                    fs.vel.z.re = buf[ii++]; version(complex_numbers) { fs.vel.z.im = buf[ii++]; }
                    version(MHD) {
                        fs.B.x.re = buf[ii++]; version(complex_numbers) { fs.B.x.im = buf[ii++]; }
                        fs.B.y.re = buf[ii++]; version(complex_numbers) { fs.B.y.im = buf[ii++]; }
                        fs.B.z.re = buf[ii++]; version(complex_numbers) { fs.B.z.im = buf[ii++]; }
                        fs.psi.re = buf[ii++]; version(complex_numbers) { fs.psi.im = buf[ii++]; }
                        fs.divB.re  = buf[ii++]; version(complex_numbers) { fs.divB.im = buf[ii++]; }
                    }
                    version(turbulence) {
                        foreach (j; 0 .. nturb) { fs.turb[j].re = buf[ii++]; version(complex_numbers) { fs.turb[j].im = buf[ii++]; } }
                    }
                    fs.mu_t.re = buf[ii++]; version(complex_numbers) { fs.mu_t.im = buf[ii++]; }
                    fs.k_t.re = buf[ii++]; version(complex_numbers) { fs.k_t.im = buf[ii++]; }
                    fs.S.re = buf[ii++]; version(complex_numbers) { fs.S.im = buf[ii++]; }
                }
            }
        } else { // not mpi_parallel
            // For a single process, just access the data directly.
            foreach (i, mygc; ghost_cells) {
                mygc.fs.copy_values_from(mapped_cells[i].fs);
            }
        }
    } // end exchange_flowstate_phase2()

    void exchange_turbulent_transprops_phase0()
    {
        version(mpi_parallel) {
            foreach (i; 0 .. n_incoming) {
                size_t fs_size = 2;
                size_t ne = incoming_ncells_list[i] * fs_size;
                version(complex_numbers) { ne *= 2; }
                if (incoming_flowstate_buf_list[i].length < ne) { incoming_flowstate_buf_list[i].length = ne; }
                // Post non-blocking receive for flowstate data that we expect to receive later
                // from the src_blk MPI process.
                MPI_Irecv(incoming_flowstate_buf_list[i].ptr, to!int(ne), MPI_DOUBLE, incoming_rank_list[i],
                          incoming_flowstate_tag_list[i], MPI_COMM_WORLD, &incoming_flowstate_request_list[i]);
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    }

    void exchange_turbulent_transprops_phase1()
    {
        version(mpi_parallel) {
            foreach (i; 0 .. n_outgoing) {
                size_t fs_size = 2;
                size_t ne = outgoing_ncells_list[i] * fs_size;
                version(complex_numbers) { ne *= 2; }
                if (outgoing_flowstate_buf_list[i].length < ne) { outgoing_flowstate_buf_list[i].length = ne; }
                auto buf = outgoing_flowstate_buf_list[i];

                size_t ii = 0;
                foreach (cid; src_cell_ids[blk.id][outgoing_block_list[i]]) {
                    auto c = blk.cells[cid];
                    FlowState* fs = &(c.fs);
                    buf[ii++] = fs.mu_t.re; version(complex_numbers) { buf[ii++] = fs.mu_t.im; }
                    buf[ii++] = fs.k_t.re; version(complex_numbers) { buf[ii++] = fs.k_t.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                              outgoing_flowstate_tag_list[i], MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                             outgoing_flowstate_tag_list[i], MPI_COMM_WORLD);
                }
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    }

    void exchange_turbulent_transprops_phase2()
    {
        version(mpi_parallel) {
            foreach (i; 0 .. n_incoming) {
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_flowstate_request_list[i], &incoming_flowstate_status_list[i]);
                } else {
                    MPI_Wait(&incoming_flowstate_request_list[i], &incoming_flowstate_status_list[i]);
                }
                auto buf = incoming_flowstate_buf_list[i];
                size_t ii = 0;
                foreach (gi; ghost_cell_indices[incoming_block_list[i]][blk.id]) {
                    auto c = ghost_cells[gi];
                    FlowState* fs = &(c.fs);
                    fs.mu_t.re = buf[ii++]; version(complex_numbers) { fs.mu_t.im = buf[ii++]; }
                    fs.k_t.re = buf[ii++]; version(complex_numbers) { fs.k_t.im = buf[ii++]; }
                }
            }
        } else { // not mpi_parallel
            // For a single process, just access the data directly.
            foreach (i, mygc; ghost_cells) {
                mygc.fs.mu_t = mapped_cells[i].fs.mu_t;
                mygc.fs.k_t = mapped_cells[i].fs.k_t;
            }
        }
    }

    @nogc
    size_t convective_gradient_buffer_entry_size(const LocalConfig myConfig)
    {
        /*
        Compute the amount of space needed for one gradient in the SEND/RECV buffer

        Note: This routine must be kept consistent with the buffer packing in
        exchange_convective_gradient phases 1 and 2
        @author: Nick N. Gibbons
        */
        size_t nspecies = myConfig.n_species;
        size_t nmodes = myConfig.n_modes;
        size_t nitems = 42;
        version(MHD) { nitems += 24; }
        version(turbulence) { nitems += myConfig.turb_model.nturb*6; }
        nitems += nmodes*12;
        nitems += nspecies*6;
        version(complex_numbers) {
            nitems *= 2;
        }
        return nitems;
    }

    // not @nogc
    void exchange_convective_gradient_phase0(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            // Prepare to exchange geometry data for the boundary cells.
            //size_t nspecies = blk.myConfig.n_species;
            //size_t nmodes = blk.myConfig.n_modes;
            foreach (i; 0 .. n_incoming) {
                // Exchange cell-centered convective gradients for the boundary cells.
                // the size of the buffer should match up with that of lsqinterp.d
                //size_t nitems = 42;
                //version(MHD) { nitems += 24; }
                //version(turbulence) { nitems += 12; }
                size_t grad_size = convective_gradient_buffer_entry_size(blk.myConfig);
                size_t ne = incoming_ncells_list[i] * grad_size;
                if (incoming_convective_gradient_buf_list[i].length < ne) { incoming_convective_gradient_buf_list[i].length = ne; }
                // Post non-blocking receive for flowstate data that we expect to receive later
                // from the src_blk MPI process.
                MPI_Irecv(incoming_convective_gradient_buf_list[i].ptr, to!int(ne), MPI_DOUBLE, incoming_rank_list[i],
                          incoming_convective_gradient_tag_list[i], MPI_COMM_WORLD, &incoming_convective_gradient_request_list[i]);
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_convective_gradient_phase0()

    // not @nogc
    void exchange_convective_gradient_phase1(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach (i; 0 .. n_outgoing) {
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                //size_t nitems = 42;
                //version(MHD) { nitems += 24; }
                //version(turbulence) { nitems += 12; }
                size_t grad_size = convective_gradient_buffer_entry_size(blk.myConfig);
                size_t ne = outgoing_ncells_list[i] * grad_size;
                if (outgoing_convective_gradient_buf_list[i].length < ne) { outgoing_convective_gradient_buf_list[i].length = ne; }
                auto buf = outgoing_convective_gradient_buf_list[i];
                size_t ii = 0;
                foreach (cid; src_cell_ids[blk.id][outgoing_block_list[i]]) {
                    auto c = blk.cells[cid].gradients;
                    // velocity
                    buf[ii++] = c.velx[0].re; version(complex_numbers) { buf[ii++] = c.velx[0].im; }
                    buf[ii++] = c.velx[1].re; version(complex_numbers) { buf[ii++] = c.velx[1].im; }
                    buf[ii++] = c.velx[2].re; version(complex_numbers) { buf[ii++] = c.velx[2].im; }
                    buf[ii++] = c.velxPhi.re; version(complex_numbers) { buf[ii++] = c.velxPhi.im; }
                    buf[ii++] = c.velxMin.re; version(complex_numbers) { buf[ii++] = c.velxMin.im; }
                    buf[ii++] = c.velxMax.re; version(complex_numbers) { buf[ii++] = c.velxMax.im; }
                    buf[ii++] = c.vely[0].re; version(complex_numbers) { buf[ii++] = c.vely[0].im; }
                    buf[ii++] = c.vely[1].re; version(complex_numbers) { buf[ii++] = c.vely[1].im; }
                    buf[ii++] = c.vely[2].re; version(complex_numbers) { buf[ii++] = c.vely[2].im; }
                    buf[ii++] = c.velyPhi.re; version(complex_numbers) { buf[ii++] = c.velyPhi.im; }
                    buf[ii++] = c.velyMin.re; version(complex_numbers) { buf[ii++] = c.velyMin.im; }
                    buf[ii++] = c.velyMax.re; version(complex_numbers) { buf[ii++] = c.velyMax.im; }
                    buf[ii++] = c.velz[0].re; version(complex_numbers) { buf[ii++] = c.velz[0].im; }
                    buf[ii++] = c.velz[1].re; version(complex_numbers) { buf[ii++] = c.velz[1].im; }
                    buf[ii++] = c.velz[2].re; version(complex_numbers) { buf[ii++] = c.velz[2].im; }
                    buf[ii++] = c.velzPhi.re; version(complex_numbers) { buf[ii++] = c.velzPhi.im; }
                    buf[ii++] = c.velzMin.re; version(complex_numbers) { buf[ii++] = c.velzMin.im; }
                    buf[ii++] = c.velzMax.re; version(complex_numbers) { buf[ii++] = c.velzMax.im; }
                    // rho, p, T, u
                    buf[ii++] = c.rho[0].re; version(complex_numbers) { buf[ii++] = c.rho[0].im; }
                    buf[ii++] = c.rho[1].re; version(complex_numbers) { buf[ii++] = c.rho[1].im; }
                    buf[ii++] = c.rho[2].re; version(complex_numbers) { buf[ii++] = c.rho[2].im; }
                    buf[ii++] = c.rhoPhi.re; version(complex_numbers) { buf[ii++] = c.rhoPhi.im; }
                    buf[ii++] = c.rhoMin.re; version(complex_numbers) { buf[ii++] = c.rhoMin.im; }
                    buf[ii++] = c.rhoMax.re; version(complex_numbers) { buf[ii++] = c.rhoMax.im; }
                    buf[ii++] = c.p[0].re; version(complex_numbers) { buf[ii++] = c.p[0].im; }
                    buf[ii++] = c.p[1].re; version(complex_numbers) { buf[ii++] = c.p[1].im; }
                    buf[ii++] = c.p[2].re; version(complex_numbers) { buf[ii++] = c.p[2].im; }
                    buf[ii++] = c.pPhi.re; version(complex_numbers) { buf[ii++] = c.pPhi.im; }
                    buf[ii++] = c.pMin.re; version(complex_numbers) { buf[ii++] = c.pMin.im; }
                    buf[ii++] = c.pMax.re; version(complex_numbers) { buf[ii++] = c.pMax.im; }
                    buf[ii++] = c.T[0].re; version(complex_numbers) { buf[ii++] = c.T[0].im; }
                    buf[ii++] = c.T[1].re; version(complex_numbers) { buf[ii++] = c.T[1].im; }
                    buf[ii++] = c.T[2].re; version(complex_numbers) { buf[ii++] = c.T[2].im; }
                    buf[ii++] = c.TPhi.re; version(complex_numbers) { buf[ii++] = c.TPhi.im; }
                    buf[ii++] = c.TMin.re; version(complex_numbers) { buf[ii++] = c.TMin.im; }
                    buf[ii++] = c.TMax.re; version(complex_numbers) { buf[ii++] = c.TMax.im; }
                    buf[ii++] = c.u[0].re; version(complex_numbers) { buf[ii++] = c.u[0].im; }
                    buf[ii++] = c.u[1].re; version(complex_numbers) { buf[ii++] = c.u[1].im; }
                    buf[ii++] = c.u[2].re; version(complex_numbers) { buf[ii++] = c.u[2].im; }
                    buf[ii++] = c.uPhi.re; version(complex_numbers) { buf[ii++] = c.uPhi.im; }
                    buf[ii++] = c.uMin.re; version(complex_numbers) { buf[ii++] = c.uMin.im; }
                    buf[ii++] = c.uMax.re; version(complex_numbers) { buf[ii++] = c.uMax.im; }
                    // formerly tke, omega
                    version(turbulence) {
                        foreach(j; 0 .. nturb) {
                            buf[ii++] = c.turb[j][0].re; version(complex_numbers) { buf[ii++] = c.turb[j][0].im; }
                            buf[ii++] = c.turb[j][1].re; version(complex_numbers) { buf[ii++] = c.turb[j][1].im; }
                            buf[ii++] = c.turb[j][2].re; version(complex_numbers) { buf[ii++] = c.turb[j][2].im; }
                            buf[ii++] = c.turbPhi[j].re; version(complex_numbers) { buf[ii++] = c.turbPhi[j].im; }
                            buf[ii++] = c.turbMin[j].re; version(complex_numbers) { buf[ii++] = c.turbMin[j].im; }
                            buf[ii++] = c.turbMax[j].re; version(complex_numbers) { buf[ii++] = c.turbMax[j].im; }
                        }
                    }
                    // MHD
                    version(MHD) {
                        buf[ii++] = c.Bx[0].re; version(complex_numbers) { buf[ii++] = c.Bx[0].im; }
                        buf[ii++] = c.Bx[1].re; version(complex_numbers) { buf[ii++] = c.Bx[1].im; }
                        buf[ii++] = c.Bx[2].re; version(complex_numbers) { buf[ii++] = c.Bx[2].im; }
                        buf[ii++] = c.BxPhi.re; version(complex_numbers) { buf[ii++] = c.BxPhi.im; }
                        buf[ii++] = c.BxMin.re; version(complex_numbers) { buf[ii++] = c.BxMin.im; }
                        buf[ii++] = c.BxMax.re; version(complex_numbers) { buf[ii++] = c.BxMax.im; }
                        buf[ii++] = c.By[0].re; version(complex_numbers) { buf[ii++] = c.By[0].im; }
                        buf[ii++] = c.By[1].re; version(complex_numbers) { buf[ii++] = c.By[1].im; }
                        buf[ii++] = c.By[2].re; version(complex_numbers) { buf[ii++] = c.By[2].im; }
                        buf[ii++] = c.ByPhi.re; version(complex_numbers) { buf[ii++] = c.ByPhi.im; }
                        buf[ii++] = c.ByMin.re; version(complex_numbers) { buf[ii++] = c.ByMin.im; }
                        buf[ii++] = c.ByMax.re; version(complex_numbers) { buf[ii++] = c.ByMax.im; }
                        buf[ii++] = c.Bz[0].re; version(complex_numbers) { buf[ii++] = c.Bz[0].im; }
                        buf[ii++] = c.Bz[1].re; version(complex_numbers) { buf[ii++] = c.Bz[1].im; }
                        buf[ii++] = c.Bz[2].re; version(complex_numbers) { buf[ii++] = c.Bz[2].im; }
                        buf[ii++] = c.BzPhi.re; version(complex_numbers) { buf[ii++] = c.BzPhi.im; }
                        buf[ii++] = c.BzMin.re; version(complex_numbers) { buf[ii++] = c.BzMin.im; }
                        buf[ii++] = c.BzMax.re; version(complex_numbers) { buf[ii++] = c.BzMax.im; }
                        buf[ii++] = c.psi[0].re; version(complex_numbers) { buf[ii++] = c.psi[0].im; }
                        buf[ii++] = c.psi[1].re; version(complex_numbers) { buf[ii++] = c.psi[1].im; }
                        buf[ii++] = c.psi[2].re; version(complex_numbers) { buf[ii++] = c.psi[2].im; }
                        buf[ii++] = c.psiPhi.re; version(complex_numbers) { buf[ii++] = c.psiPhi.im; }
                        buf[ii++] = c.psiMin.re; version(complex_numbers) { buf[ii++] = c.psiMin.im; }
                        buf[ii++] = c.psiMax.re; version(complex_numbers) { buf[ii++] = c.psiMax.im; }
                    }
                    // multi-species
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) {
                            buf[ii++] = c.rho_s[j][0].re; version(complex_numbers) { buf[ii++] = c.rho_s[j][0].im; }
                            buf[ii++] = c.rho_s[j][1].re; version(complex_numbers) { buf[ii++] = c.rho_s[j][1].im; }
                            buf[ii++] = c.rho_s[j][2].re; version(complex_numbers) { buf[ii++] = c.rho_s[j][2].im; }
                            buf[ii++] = c.rho_sPhi[j].re; version(complex_numbers) { buf[ii++] = c.rho_sPhi[j].im; }
                            buf[ii++] = c.rho_sMin[j].re; version(complex_numbers) { buf[ii++] = c.rho_sMin[j].im; }
                            buf[ii++] = c.rho_sMax[j].re; version(complex_numbers) { buf[ii++] = c.rho_sMax[j].im; }
                        }
                    }
                    // multi-T
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) {
                            buf[ii++] = c.T_modes[j][0].re; version(complex_numbers) { buf[ii++] = c.T_modes[j][0].im; }
                            buf[ii++] = c.T_modes[j][1].re; version(complex_numbers) { buf[ii++] = c.T_modes[j][1].im; }
                            buf[ii++] = c.T_modes[j][2].re; version(complex_numbers) { buf[ii++] = c.T_modes[j][2].im; }
                            buf[ii++] = c.T_modesPhi[j].re; version(complex_numbers) { buf[ii++] = c.T_modesPhi[j].im; }
                            buf[ii++] = c.T_modesMin[j].re; version(complex_numbers) { buf[ii++] = c.T_modesMin[j].im; }
                            buf[ii++] = c.T_modesMax[j].re; version(complex_numbers) { buf[ii++] = c.T_modesMax[j].im; }
                        }
                        foreach (j; 0 .. nmodes) {
                            buf[ii++] = c.u_modes[j][0].re; version(complex_numbers) { buf[ii++] = c.u_modes[j][0].im; }
                            buf[ii++] = c.u_modes[j][1].re; version(complex_numbers) { buf[ii++] = c.u_modes[j][1].im; }
                            buf[ii++] = c.u_modes[j][2].re; version(complex_numbers) { buf[ii++] = c.u_modes[j][2].im; }
                            buf[ii++] = c.u_modesPhi[j].re; version(complex_numbers) { buf[ii++] = c.u_modesPhi[j].im; }
                            buf[ii++] = c.u_modesMin[j].re; version(complex_numbers) { buf[ii++] = c.u_modesMin[j].im; }
                            buf[ii++] = c.u_modesMax[j].re; version(complex_numbers) { buf[ii++] = c.u_modesMax[j].im; }
                        }
                    }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                              outgoing_convective_gradient_tag_list[i], MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                             outgoing_convective_gradient_tag_list[i], MPI_COMM_WORLD);
                }
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_convective_gradient_phase1()

    // not @nogc
    void exchange_convective_gradient_phase2(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach (i; 0 .. n_incoming) {
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_convective_gradient_request_list[i], &incoming_convective_gradient_status_list[i]);
                } else {
                    MPI_Wait(&incoming_convective_gradient_request_list[i], &incoming_convective_gradient_status_list[i]);
                }
                auto buf = incoming_convective_gradient_buf_list[i];
                size_t ii = 0;
                foreach (gi; ghost_cell_indices[incoming_block_list[i]][blk.id]) {
                    auto c = ghost_cells[gi].gradients;
                    // velocity
                    c.velx[0].re = buf[ii++]; version(complex_numbers) { c.velx[0].im = buf[ii++]; }
                    c.velx[1].re = buf[ii++]; version(complex_numbers) { c.velx[1].im = buf[ii++]; }
                    c.velx[2].re = buf[ii++]; version(complex_numbers) { c.velx[2].im = buf[ii++]; }
                    c.velxPhi.re = buf[ii++]; version(complex_numbers) { c.velxPhi.im = buf[ii++]; }
                    c.velxMin.re = buf[ii++]; version(complex_numbers) { c.velxMin.im = buf[ii++]; }
                    c.velxMax.re = buf[ii++]; version(complex_numbers) { c.velxMax.im = buf[ii++]; }
                    c.vely[0].re = buf[ii++]; version(complex_numbers) { c.vely[0].im = buf[ii++]; }
                    c.vely[1].re = buf[ii++]; version(complex_numbers) { c.vely[1].im = buf[ii++]; }
                    c.vely[2].re = buf[ii++]; version(complex_numbers) { c.vely[2].im = buf[ii++]; }
                    c.velyPhi.re = buf[ii++]; version(complex_numbers) { c.velyPhi.im = buf[ii++]; }
                    c.velyMin.re = buf[ii++]; version(complex_numbers) { c.velyMin.im = buf[ii++]; }
                    c.velyMax.re = buf[ii++]; version(complex_numbers) { c.velyMax.im = buf[ii++]; }
                    c.velz[0].re = buf[ii++]; version(complex_numbers) { c.velz[0].im = buf[ii++]; }
                    c.velz[1].re = buf[ii++]; version(complex_numbers) { c.velz[1].im = buf[ii++]; }
                    c.velz[2].re = buf[ii++]; version(complex_numbers) { c.velz[2].im = buf[ii++]; }
                    c.velzPhi.re = buf[ii++]; version(complex_numbers) { c.velzPhi.im = buf[ii++]; }
                    c.velzMin.re = buf[ii++]; version(complex_numbers) { c.velzMin.im = buf[ii++]; }
                    c.velzMax.re = buf[ii++]; version(complex_numbers) { c.velzMax.im = buf[ii++]; }
                    // rho, p, T, u
                    c.rho[0].re = buf[ii++]; version(complex_numbers) { c.rho[0].im = buf[ii++]; }
                    c.rho[1].re = buf[ii++]; version(complex_numbers) { c.rho[1].im = buf[ii++]; }
                    c.rho[2].re = buf[ii++]; version(complex_numbers) { c.rho[2].im = buf[ii++]; }
                    c.rhoPhi.re = buf[ii++]; version(complex_numbers) { c.rhoPhi.im = buf[ii++]; }
                    c.rhoMin.re = buf[ii++]; version(complex_numbers) { c.rhoMin.im = buf[ii++]; }
                    c.rhoMax.re = buf[ii++]; version(complex_numbers) { c.rhoMax.im = buf[ii++]; }
                    c.p[0].re = buf[ii++]; version(complex_numbers) { c.p[0].im = buf[ii++]; }
                    c.p[1].re = buf[ii++]; version(complex_numbers) { c.p[1].im = buf[ii++]; }
                    c.p[2].re = buf[ii++]; version(complex_numbers) { c.p[2].im = buf[ii++]; }
                    c.pPhi.re = buf[ii++]; version(complex_numbers) { c.pPhi.im = buf[ii++]; }
                    c.pMin.re = buf[ii++]; version(complex_numbers) { c.pMin.im = buf[ii++]; }
                    c.pMax.re = buf[ii++]; version(complex_numbers) { c.pMax.im = buf[ii++]; }
                    c.T[0].re = buf[ii++]; version(complex_numbers) { c.T[0].im = buf[ii++]; }
                    c.T[1].re = buf[ii++]; version(complex_numbers) { c.T[1].im = buf[ii++]; }
                    c.T[2].re = buf[ii++]; version(complex_numbers) { c.T[2].im = buf[ii++]; }
                    c.TPhi.re = buf[ii++]; version(complex_numbers) { c.TPhi.im = buf[ii++]; }
                    c.TMin.re = buf[ii++]; version(complex_numbers) { c.TMin.im = buf[ii++]; }
                    c.TMax.re = buf[ii++]; version(complex_numbers) { c.TMax.im = buf[ii++]; }
                    c.u[0].re = buf[ii++]; version(complex_numbers) { c.u[0].im = buf[ii++]; }
                    c.u[1].re = buf[ii++]; version(complex_numbers) { c.u[1].im = buf[ii++]; }
                    c.u[2].re = buf[ii++]; version(complex_numbers) { c.u[2].im = buf[ii++]; }
                    c.uPhi.re = buf[ii++]; version(complex_numbers) { c.uPhi.im = buf[ii++]; }
                    c.uMin.re = buf[ii++]; version(complex_numbers) { c.uMin.im = buf[ii++]; }
                    c.uMax.re = buf[ii++]; version(complex_numbers) { c.uMax.im = buf[ii++]; }
                    // tke, omega
                    version(turbulence) {
                        foreach (j; 0 .. nturb){
                            c.turb[j][0].re = buf[ii++]; version(complex_numbers) { c.turb[j][0].im = buf[ii++]; }
                            c.turb[j][1].re = buf[ii++]; version(complex_numbers) { c.turb[j][1].im = buf[ii++]; }
                            c.turb[j][2].re = buf[ii++]; version(complex_numbers) { c.turb[j][2].im = buf[ii++]; }
                            c.turbPhi[j].re = buf[ii++]; version(complex_numbers) { c.turbPhi[j].im = buf[ii++]; }
                            c.turbMin[j].re = buf[ii++]; version(complex_numbers) { c.turbMin[j].im = buf[ii++]; }
                            c.turbMax[j].re = buf[ii++]; version(complex_numbers) { c.turbMax[j].im = buf[ii++]; }
                        }
                    }
                    // MHD
                    version(MHD) {
                        c.Bx[0].re = buf[ii++]; version(complex_numbers) { c.Bx[0].im = buf[ii++]; }
                        c.Bx[1].re = buf[ii++]; version(complex_numbers) { c.Bx[1].im = buf[ii++]; }
                        c.Bx[2].re = buf[ii++]; version(complex_numbers) { c.Bx[2].im = buf[ii++]; }
                        c.BxPhi.re = buf[ii++]; version(complex_numbers) { c.BxPhi.im = buf[ii++]; }
                        c.BxMin.re = buf[ii++]; version(complex_numbers) { c.BxMin.im = buf[ii++]; }
                        c.BxMax.re = buf[ii++]; version(complex_numbers) { c.BxMax.im = buf[ii++]; }
                        c.By[0].re = buf[ii++]; version(complex_numbers) { c.By[0].im = buf[ii++]; }
                        c.By[1].re = buf[ii++]; version(complex_numbers) { c.By[1].im = buf[ii++]; }
                        c.By[2].re = buf[ii++]; version(complex_numbers) { c.By[2].im = buf[ii++]; }
                        c.ByPhi.re = buf[ii++]; version(complex_numbers) { c.ByPhi.im = buf[ii++]; }
                        c.ByMin.re = buf[ii++]; version(complex_numbers) { c.ByMin.im = buf[ii++]; }
                        c.ByMax.re = buf[ii++]; version(complex_numbers) { c.ByMax.im = buf[ii++]; }
                        c.Bz[0].re = buf[ii++]; version(complex_numbers) { c.Bz[0].im = buf[ii++]; }
                        c.Bz[1].re = buf[ii++]; version(complex_numbers) { c.Bz[1].im = buf[ii++]; }
                        c.Bz[2].re = buf[ii++]; version(complex_numbers) { c.Bz[2].im = buf[ii++]; }
                        c.BzPhi.re = buf[ii++]; version(complex_numbers) { c.BzPhi.im = buf[ii++]; }
                        c.BzMin.re = buf[ii++]; version(complex_numbers) { c.BzMin.im = buf[ii++]; }
                        c.BzMax.re = buf[ii++]; version(complex_numbers) { c.BzMax.im = buf[ii++]; }
                        c.psi[0].re = buf[ii++]; version(complex_numbers) {c.psi[0].im = buf[ii++]; }
                        c.psi[1].re = buf[ii++]; version(complex_numbers) {c.psi[1].im = buf[ii++]; }
                        c.psi[2].re = buf[ii++]; version(complex_numbers) {c.psi[2].im = buf[ii++]; }
                        c.psiPhi.re = buf[ii++]; version(complex_numbers) {c.psiPhi.im = buf[ii++]; }
                        c.psiMin.re = buf[ii++]; version(complex_numbers) {c.psiMin.im = buf[ii++]; }
                        c.psiMax.re = buf[ii++]; version(complex_numbers) {c.psiMax.im = buf[ii++]; }
                    }
                    // multi-species
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) {
                            c.rho_s[j][0].re = buf[ii++]; version(complex_numbers) { c.rho_s[j][0].im = buf[ii++]; }
                            c.rho_s[j][1].re = buf[ii++]; version(complex_numbers) { c.rho_s[j][1].im = buf[ii++]; }
                            c.rho_s[j][2].re = buf[ii++]; version(complex_numbers) { c.rho_s[j][2].im = buf[ii++]; }
                            c.rho_sPhi[j].re = buf[ii++]; version(complex_numbers) { c.rho_sPhi[j].im = buf[ii++]; }
                            c.rho_sMin[j].re = buf[ii++]; version(complex_numbers) { c.rho_sMin[j].im = buf[ii++]; }
                            c.rho_sMax[j].re = buf[ii++]; version(complex_numbers) { c.rho_sMax[j].im = buf[ii++]; }
                        }
                    }
                    // multi-T
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) {
                            c.T_modes[j][0].re = buf[ii++]; version(complex_numbers) { c.T_modes[j][0].im = buf[ii++]; }
                            c.T_modes[j][1].re = buf[ii++]; version(complex_numbers) { c.T_modes[j][1].im = buf[ii++]; }
                            c.T_modes[j][2].re = buf[ii++]; version(complex_numbers) { c.T_modes[j][2].im = buf[ii++]; }
                            c.T_modesPhi[j].re = buf[ii++]; version(complex_numbers) { c.T_modesPhi[j].im = buf[ii++]; }
                            c.T_modesMin[j].re = buf[ii++]; version(complex_numbers) { c.T_modesMin[j].im = buf[ii++]; }
                            c.T_modesMax[j].re = buf[ii++]; version(complex_numbers) { c.T_modesMax[j].im = buf[ii++]; }
                        }
                        foreach (j; 0 .. nmodes) {
                            c.u_modes[j][0].re = buf[ii++]; version(complex_numbers) { c.u_modes[j][0].im = buf[ii++]; }
                            c.u_modes[j][1].re = buf[ii++]; version(complex_numbers) { c.u_modes[j][1].im = buf[ii++]; }
                            c.u_modes[j][2].re = buf[ii++]; version(complex_numbers) { c.u_modes[j][2].im = buf[ii++]; }
                            c.u_modesPhi[j].re = buf[ii++]; version(complex_numbers) { c.u_modesPhi[j].im = buf[ii++]; }
                            c.u_modesMin[j].re = buf[ii++]; version(complex_numbers) { c.u_modesMin[j].im = buf[ii++]; }
                            c.u_modesMax[j].re = buf[ii++]; version(complex_numbers) { c.u_modesMax[j].im = buf[ii++]; }
                        }
                    }
                }
            }
        } else { // not mpi_parallel
            // For a single process, just access the data directly.
            foreach (i, mygc; ghost_cells) {
                mygc.gradients.copy_values_from(*(mapped_cells[i].gradients));
            }
        }
    } // end exchange_convective_gradient_phase2()

    @nogc
    size_t viscous_gradient_buffer_entry_size(const LocalConfig myConfig)
    {
        /*
        Compute the amount of space needed for one gradient in the SEND/RECV buffer

        Note: This routine must be kept consistent with the buffer packing in exchange_viscous_gradient
        phases 1 and 2
        @author: Nick N. Gibbons
        */
        size_t nspecies = myConfig.n_species;
        size_t nmodes = myConfig.n_modes;
        size_t nitems = 12;
        version(turbulence) { nitems += myConfig.turb_model.nturb*3; }
        nitems += nmodes*3;
        nitems += nspecies*3;
        version(complex_numbers) {
            nitems *= 2;
        }
        return nitems;
    }

    // not @nogc
    void exchange_viscous_gradient_phase0(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            // Prepare to exchange geometry data for the boundary cells.
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            foreach (i; 0 .. n_incoming) {
                // Exchange cell-centered viscous gradients for the boundary cells.
                // the size of the buffer should match up with that of lsqinterp.d
                //size_t nitems = 12;
                //version(turbulence) { nitems += 6; }
                size_t grad_size = viscous_gradient_buffer_entry_size(blk.myConfig);
                size_t ne = incoming_ncells_list[i] * grad_size;
                if (incoming_viscous_gradient_buf_list[i].length < ne) { incoming_viscous_gradient_buf_list[i].length = ne; }
                // Post non-blocking receive for flowstate data that we expect to receive later
                // from the src_blk MPI process.
                MPI_Irecv(incoming_viscous_gradient_buf_list[i].ptr, to!int(ne), MPI_DOUBLE, incoming_rank_list[i],
                          incoming_viscous_gradient_tag_list[i], MPI_COMM_WORLD, &incoming_viscous_gradient_request_list[i]);
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_viscous_gradient_phase0()

    // not @nogc
    void exchange_viscous_gradient_phase1(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach (i; 0 .. n_outgoing) {
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                //size_t nitems = 12;
                //version(turbulence) { nitems += 6; }
                size_t grad_size = viscous_gradient_buffer_entry_size(blk.myConfig);
                size_t ne = outgoing_ncells_list[i] * grad_size;
                if (outgoing_viscous_gradient_buf_list[i].length < ne) { outgoing_viscous_gradient_buf_list[i].length = ne; }
                auto buf = outgoing_viscous_gradient_buf_list[i];
                size_t ii = 0;
                foreach (cid; src_cell_ids[blk.id][outgoing_block_list[i]]) {
                    auto c = blk.cells[cid].grad;
                    // velocity
                    buf[ii++] = c.vel[0][0].re; version(complex_numbers) { buf[ii++] = c.vel[0][0].im; }
                    buf[ii++] = c.vel[0][1].re; version(complex_numbers) { buf[ii++] = c.vel[0][1].im; }
                    buf[ii++] = c.vel[0][2].re; version(complex_numbers) { buf[ii++] = c.vel[0][2].im; }
                    buf[ii++] = c.vel[1][0].re; version(complex_numbers) { buf[ii++] = c.vel[1][0].im; }
                    buf[ii++] = c.vel[1][1].re; version(complex_numbers) { buf[ii++] = c.vel[1][1].im; }
                    buf[ii++] = c.vel[1][2].re; version(complex_numbers) { buf[ii++] = c.vel[1][2].im; }
                    buf[ii++] = c.vel[2][0].re; version(complex_numbers) { buf[ii++] = c.vel[2][0].im; }
                    buf[ii++] = c.vel[2][1].re; version(complex_numbers) { buf[ii++] = c.vel[2][1].im; }
                    buf[ii++] = c.vel[2][2].re; version(complex_numbers) { buf[ii++] = c.vel[2][2].im; }
                    // T
                    buf[ii++] = c.T[0].re; version(complex_numbers) { buf[ii++] = c.T[0].im; }
                    buf[ii++] = c.T[1].re; version(complex_numbers) { buf[ii++] = c.T[1].im; }
                    buf[ii++] = c.T[2].re; version(complex_numbers) { buf[ii++] = c.T[2].im; }
                    // tke, omega
                    version(turbulence) {
                        foreach (j; 0 .. nturb) {
                            buf[ii++] = c.turb[j][0].re; version(complex_numbers) { buf[ii++] = c.turb[j][0].im; }
                            buf[ii++] = c.turb[j][1].re; version(complex_numbers) { buf[ii++] = c.turb[j][1].im; }
                            buf[ii++] = c.turb[j][2].re; version(complex_numbers) { buf[ii++] = c.turb[j][2].im; }
                        }
                    }
                    // multi-species
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) {
                            buf[ii++] = c.massf[j][0].re; version(complex_numbers) { buf[ii++] = c.massf[j][0].im; }
                            buf[ii++] = c.massf[j][1].re; version(complex_numbers) { buf[ii++] = c.massf[j][1].im; }
                            buf[ii++] = c.massf[j][2].re; version(complex_numbers) { buf[ii++] = c.massf[j][2].im; }
                        }
                    }
                    // multi-T
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) {
                            buf[ii++] = c.T_modes[j][0].re; version(complex_numbers) { buf[ii++] = c.T_modes[j][0].im; }
                            buf[ii++] = c.T_modes[j][1].re; version(complex_numbers) { buf[ii++] = c.T_modes[j][1].im; }
                            buf[ii++] = c.T_modes[j][2].re; version(complex_numbers) { buf[ii++] = c.T_modes[j][2].im; }
                        }
                    }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                              outgoing_viscous_gradient_tag_list[i], MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(buf.ptr, to!int(ne), MPI_DOUBLE, outgoing_rank_list[i],
                             outgoing_viscous_gradient_tag_list[i], MPI_COMM_WORLD);
                }
            }
        } else { // not mpi_parallel
            // For a single process, nothing to be done because
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_viscous_gradient_phase1()

    // not @nogc
    void exchange_viscous_gradient_phase2(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            size_t nspecies = blk.myConfig.n_species;
            size_t nmodes = blk.myConfig.n_modes;
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach (i; 0 .. n_incoming) {
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_viscous_gradient_request_list[i], &incoming_viscous_gradient_status_list[i]);
                } else {
                    MPI_Wait(&incoming_viscous_gradient_request_list[i], &incoming_viscous_gradient_status_list[i]);
                }
                auto buf = incoming_viscous_gradient_buf_list[i];
                size_t ii = 0;
                foreach (gi; ghost_cell_indices[incoming_block_list[i]][blk.id]) {
                    auto c = ghost_cells[gi].grad;
                    // velocity
                    c.vel[0][0].re = buf[ii++]; version(complex_numbers) { c.vel[0][0].im = buf[ii++]; }
                    c.vel[0][1].re = buf[ii++]; version(complex_numbers) { c.vel[0][1].im = buf[ii++]; }
                    c.vel[0][2].re = buf[ii++]; version(complex_numbers) { c.vel[0][2].im = buf[ii++]; }
                    c.vel[1][0].re = buf[ii++]; version(complex_numbers) { c.vel[1][0].im = buf[ii++]; }
                    c.vel[1][1].re = buf[ii++]; version(complex_numbers) { c.vel[1][1].im = buf[ii++]; }
                    c.vel[1][2].re = buf[ii++]; version(complex_numbers) { c.vel[1][2].im = buf[ii++]; }
                    c.vel[2][0].re = buf[ii++]; version(complex_numbers) { c.vel[2][0].im = buf[ii++]; }
                    c.vel[2][1].re = buf[ii++]; version(complex_numbers) { c.vel[2][1].im = buf[ii++]; }
                    c.vel[2][2].re = buf[ii++]; version(complex_numbers) { c.vel[2][2].im = buf[ii++]; }
                    // T
                    c.T[0].re = buf[ii++]; version(complex_numbers) { c.T[0].im = buf[ii++]; }
                    c.T[1].re = buf[ii++]; version(complex_numbers) { c.T[1].im = buf[ii++]; }
                    c.T[2].re = buf[ii++]; version(complex_numbers) { c.T[2].im = buf[ii++]; }
                    // tke, omega
                    version(turbulence) {
                        foreach (j; 0 .. nturb){
                            c.turb[j][0].re = buf[ii++]; version(complex_numbers) { c.turb[j][0].im = buf[ii++]; }
                            c.turb[j][1].re = buf[ii++]; version(complex_numbers) { c.turb[j][1].im = buf[ii++]; }
                            c.turb[j][2].re = buf[ii++]; version(complex_numbers) { c.turb[j][2].im = buf[ii++]; }
                        }
                    }
                    // multi-species
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) {
                            c.massf[j][0].re = buf[ii++]; version(complex_numbers) { c.massf[j][0].im = buf[ii++]; }
                            c.massf[j][1].re = buf[ii++]; version(complex_numbers) { c.massf[j][1].im = buf[ii++]; }
                            c.massf[j][2].re = buf[ii++]; version(complex_numbers) { c.massf[j][2].im = buf[ii++]; }
                        }
                    }
                    // multi-T
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) {
                            c.T_modes[j][0].re = buf[ii++]; version(complex_numbers) { c.T_modes[j][0].im = buf[ii++]; }
                            c.T_modes[j][1].re = buf[ii++]; version(complex_numbers) { c.T_modes[j][1].im = buf[ii++]; }
                            c.T_modes[j][2].re = buf[ii++]; version(complex_numbers) { c.T_modes[j][2].im = buf[ii++]; }
                        }
                    }
                }
            }
        } else { // not mpi_parallel
            // For a single process, just access the data directly.
            foreach (i, mygc; ghost_cells) {
                mygc.grad.copy_values_from(*(mapped_cells[i].grad));
            }
        }
    } // end exchange_viscous_gradient_phase2()

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	assert(0, "apply_for_interface_unstructured_grid not implemented for this BC.");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        // We presume that all of the exchange of data happened earlier,
        // and that the ghost cells have been filled with flow state data
        // from their respective source cells.
        foreach (i, mygc; ghost_cells) {
            if (reorient_vector_quantities) {
                mygc.fs.reorient_vector_quantities(Rmatrix);
            }
            // [TODO] PJ 2018-01-14 If unstructured blocks ever get used in
            // the block-marching process, we will need a call to encode_conserved
            // at this point.  See the GhostCellFullFaceCopy class.
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	assert(0, "apply_for_interface_structured_grid not implemented for this BC.");
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        foreach (i, mygc; ghost_cells) {
            if (reorient_vector_quantities) {
                mygc.fs.reorient_vector_quantities(Rmatrix);
            }
        }
    } // end apply_unstructured_grid()
} // end class GhostCellMappedCellCopy
