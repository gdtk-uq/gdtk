module lmr.bc.boundary_vertex_full_face_copy;

import std.algorithm;
import std.conv;
import std.math;
import std.string;

import lmr.fvvertex;
import lmr.fluidblock;
import lmr.sfluidblock;
import lmr.globaldata;
import lmr.globalconfig;
import lmr.bc;
import geom;
import lmr.bc.ghost_cell_effect.full_face_copy : make_mpi_tag;

version(mpi_parallel) {
    import mpi;
    import lmr.bc.ghost_cell_effect.full_face_copy : MPI_Wait_a_while;
}

// This class manages exchanging vertex positions along block boundaries.
// In steady-state shock fitting simulations, without exchanging this
// information the blocks may drift apart from one another.
// Unlike the other full-face exchanges, the communication here is one way.
// The vertex positions in one block are considered the truth, and they
// get communicated to the other both. This lets the steady-state solver
// only solve for the vertex position of the vertices once.
// 
// Author: Robert Watt (copied heavily from full_face_copy.d)
// Date: 15/11/2024
class BoundaryVertexFullFaceCopy {
public:
    FVVertex[] mapped_vertices;
    size_t[] mapped_vertex_ids;
    FVVertex[] this_vertices;

    // Keeps track of if this boundary is receiving data or 
    // sending data, since vertex positions are only sent one way
    bool incoming;

    int which_boundary;
    int other_boundary;

    string desc;
    
    SFluidBlock this_blk;
    SFluidBlock other_blk;

    FluidBlock blk;

    version(mpi_parallel) {
        size_t[] outgoing_mapped_vertex_ids;
        FVVertex[] outgoing_mapped_vertices;
        int other_blk_rank;
        int outgoing_vertex_ids_tag, incoming_vertex_ids_tag;
        MPI_Request incoming_vertex_ids_request;
        MPI_Status incoming_vertex_ids_status;
        int[] outgoing_vertex_ids_buf, incoming_vertex_ids_buf;

        int outgoing_pos_tag, incoming_pos_tag;
        MPI_Request incoming_pos_request;
        MPI_Status incoming_pos_status;
        double[] outgoing_pos_buf, incoming_pos_buf;
    }


    this(int id, int boundary, int otherBlock, int otherFace) {
        which_boundary = boundary;
        other_boundary = otherFace;
        desc = "BoundaryVertexFullFaceCopy";

        blk = cast(FluidBlock) globalBlocks[id];
        if (!blk) { throw new Error("Destination Block must be FluidBlock"); }
        
        this_blk = cast(SFluidBlock) globalBlocks[id];
        if (!this_blk) { throw new Error("Destination FluidBlock must be a structured-grid block"); }

        other_blk = cast(SFluidBlock) globalBlocks[otherBlock];
        if (!other_blk) { throw new Error("Source FluidBlock must be a structured-grid block"); }

        // data is sent from the block with the lower id to the block with the higher id
        incoming = (otherBlock < id);
    }

    void setup_vertex_mapping_phase0() {
        size_t i_src, j_src, k_src;
        size_t i_dest, j_dest, k_dest;
        this_blk = cast(SFluidBlock) blk;

        version(mpi_parallel) {
            other_blk_rank = GlobalConfig.mpi_rank_for_block[other_blk.id];
        }

        if (this_blk.myConfig.dimensions == 2) {
            switch (which_boundary) {
                case Face.north:
                    j_dest = this_blk.njv - 1;
                    foreach (i; 0 .. this_blk.niv) {
                        i_dest = i;
                        this_vertices ~= this_blk.get_vtx(i_dest, j_dest);

                        switch (other_boundary) {
                            case Face.north:
                                i_src = other_blk.niv - 1 - i;
                                j_src = other_blk.njv - 1;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.east:
                                i_src = other_blk.niv - 1;
                                j_src = i;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.south:
                                i_src = i;
                                j_src = 0;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.west:
                                i_src = 0;
                                j_src = other_blk.njv - i - 1;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            default:
                                throw new FlowSolverException("Incorrect boundary connection, source face.");
                        } // end switch other_boundary
                    } // i loop
                    break;

                case Face.east:
                    i_dest = this_blk.niv - 1;
                    foreach (j; 0 .. this_blk.njv) {
                        j_dest = j;
                        this_vertices ~= this_blk.get_vtx(i_dest, j_dest);

                        switch (other_boundary) {
                            case Face.north:
                                i_src = other_blk.njv - 1;
                                j_src = j;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.east:
                                i_src = other_blk.niv - 1;
                                j_src = other_blk.njv - 1 - j;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.south:
                                i_src = other_blk.niv - j - 1;
                                j_src = 0;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.west:
                                i_src = 0;
                                j_src = j;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            default:
                                throw new FlowSolverException("Incorrect boundary connection, source face.");
                        } // end switch other_boundary
                    } // j loop
                    break;
                case Face.south:
                    j_dest = 0;
                    foreach (i; 0 .. this_blk.niv) {
                        i_dest = i;
                        this_vertices ~= this_blk.get_vtx(i_dest, j_dest);

                        switch (other_boundary) {
                            case Face.north:
                                i_src = i;
                                j_src = other_blk.njv - 1;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.east:
                                i_src = other_blk.niv - 1;
                                j_src = other_blk.njv - i - 1;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.south:
                                i_src = other_blk.niv - i - 1;
                                j_src = 0;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.west:
                                i_src = 0;
                                j_src = i;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            default:
                                throw new FlowSolverException("Incorrect boundary connection, source face.");
                        } // end switch other_boundary
                    } // i loop
                    break;
                case Face.west:
                    i_dest = 0;
                    foreach (j; 0 .. this_blk.njv) {
                        j_dest = j;
                        this_vertices ~= this_blk.get_vtx(i_dest, j_dest);

                        switch (other_boundary) {
                            case Face.north:
                                i_src = other_blk.niv - j - 1;
                                j_src = other_blk.njv - 1;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.east:
                                i_src = other_blk.niv - 1;
                                j_src = j;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.south:
                                i_src = j;
                                j_src = 0;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            case Face.west:
                                i_src = 0;
                                j_src = other_blk.njv - j - 1;
                                mapped_vertex_ids ~= other_blk.vertex_index(i_src, j_src);
                                break;
                            default:
                                throw new FlowSolverException("Incorrect boundary connection, source face.");
                        } // end switch other_boundary
                    } // j loop
                    break;
                default:
                    throw new FlowSolverException("Incorrect boundary connection, destination face.");
            }          
        } else {
            throw new Error("BoundaryVertexFullFaceCopy not implemented in 3D yet");
        }
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, so we get the data via messages.
                size_t ne = this_vertices.length;
                if (incoming_vertex_ids_buf.length < ne) { incoming_vertex_ids_buf.length = ne; }

                incoming_vertex_ids_tag = make_mpi_tag(other_blk.id, other_boundary, 0);
                MPI_Irecv(incoming_vertex_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                          incoming_vertex_ids_tag, MPI_COMM_WORLD, &incoming_vertex_ids_request);
            } else {
                // The other block happens to be in this MPI process so
                // we can directly access the memory in the final phase
            }
        } else {
            // shared memory -- we can directly acces the memory in the last step
        }

        // flag vertices not to be solved
        if (incoming) {
            foreach (ref v; this_vertices) v.solve_position = false;
        }
    }

    void setup_vertex_mapping_phase1() {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, so we get the data via MPI
                size_t ne = this_vertices.length;
                if (outgoing_vertex_ids_buf.length < ne) { outgoing_vertex_ids_buf.length = ne; }
                assert(ne == mapped_vertex_ids.length, "wrong length");
                outgoing_vertex_ids_tag = make_mpi_tag(blk.id, which_boundary, 0);
                foreach (i; 0 .. ne) {
                    outgoing_vertex_ids_buf[i] = to!int(mapped_vertex_ids[i]);
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_vertex_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                              outgoing_vertex_ids_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_vertex_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                             outgoing_vertex_ids_tag, MPI_COMM_WORLD);
                }
            } else {
                // The other block happens to be in this MPI process,
                // so we can directly access the memory in the final phase
            }
        } else {
            // shared memory -- we can directly acces the memory in the last step
        }
    }

    void setup_vertex_mapping_phase2() {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI, so we communicate with it over MPI
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_vertex_ids_request, &incoming_vertex_ids_status);
                } else {
                    MPI_Wait(&incoming_vertex_ids_request, &incoming_vertex_ids_status);
                }
                size_t ne = this_vertices.length;
                outgoing_mapped_vertex_ids.length = ne;
                foreach (i; 0 .. ne) {
                    outgoing_mapped_vertex_ids[i] = to!size_t(incoming_vertex_ids_buf[i]);
                }
                outgoing_mapped_vertices.length = 0;
                foreach (id; outgoing_mapped_vertex_ids) {
                    outgoing_mapped_vertices ~= this_blk.vertices[id];
                }
                assert(outgoing_mapped_vertices.length == this_vertices.length,
                       "mismatch in outgoing_mapped_vertices and this_vertices");
            } else {
                // The other block happens to be in this MPI process so we can
                // directly access the memory
                foreach (i; 0 .. this_vertices.length) {
                    mapped_vertices ~= other_blk.vertices[mapped_vertex_ids[i]];
                }
            }
        } else {
            // shared memory -- we can directly acces the memory
            foreach (i; 0 .. this_vertices.length) {
                mapped_vertices ~= other_blk.vertices[mapped_vertex_ids[i]];
            }
        }
    }

    void exchange_vertex_pos_phase0(int gtl) {
        if (!incoming) return; // we don't need to receive and data

        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, so communicate with it via MPI
                size_t ne = this_vertices.length * 3;
                version(complex_numbers) { ne *= 2; }
                if (incoming_pos_buf.length < ne) { incoming_pos_buf.length = ne; }

                incoming_pos_tag = make_mpi_tag(other_blk.id, other_boundary, 1);
                MPI_Irecv(incoming_pos_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_pos_tag, MPI_COMM_WORLD, &incoming_pos_request);
            } else {
                // The other block happens to be in this MPI process, so we can
                // directly access the memory in the final phase
            }
        } else {
            // shared memory -- we can directly access the memory in the last step
        }
    }

    void exchange_vertex_pos_phase1(int gtl) {
        if (incoming) return; // we don't need to send any data
        version(mpi_parallel) {

            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                outgoing_pos_tag = make_mpi_tag(blk.id, which_boundary, 1);
                size_t ne = this_vertices.length * 3;
                version(complex_numbers) { ne *= 2; }
                if (outgoing_pos_buf.length < ne) { outgoing_pos_buf.length = ne; }
                size_t ii = 0;
                foreach (v; outgoing_mapped_vertices) {
                    outgoing_pos_buf[ii++] = v.pos[gtl].x.re;
                    version(complex_numbers) { outgoing_pos_buf[ii++] = v.pos[gtl].x.im; }
                    outgoing_pos_buf[ii++] = v.pos[gtl].y.re;
                    version(complex_numbers) { outgoing_pos_buf[ii++] = v.pos[gtl].y.im; }
                    outgoing_pos_buf[ii++] = v.pos[gtl].z.re;
                    version(complex_numbers) { outgoing_pos_buf[ii++] = v.pos[gtl].z.im; }
                }

                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_pos_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_pos_tag, MPI_COMM_WORLD, &send_requent);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_pos_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_pos_tag, MPI_COMM_WORLD);
                }
            } else {
                // The other block is in this MPI process -- we can access the data directly
                // in the final phase
            }
        } else {
            // shared memory -- we can directly access the memory in the last step
        }
    }

    void exchange_vertex_pos_phase2(int gtl) {
        if (!incoming) return; // we don't need to receive any data

        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_pos_request, &incoming_pos_status);
                } else {
                    MPI_Wait(&incoming_pos_request, &incoming_pos_status);
                }
                size_t ii = 0;
                foreach (v; this_vertices) {
                    v.pos[gtl].x.re = incoming_pos_buf[ii++];
                    version(complex_numbers) { v.pos[gtl].x.im = incoming_pos_buf[ii++]; }
                    v.pos[gtl].y.re = incoming_pos_buf[ii++];
                    version(complex_numbers) { v.pos[gtl].y.im = incoming_pos_buf[ii++]; }
                    v.pos[gtl].z.re = incoming_pos_buf[ii++];
                    version(complex_numbers) { v.pos[gtl].z.im = incoming_pos_buf[ii++]; }
                }
            } else {
                // The other block happens to be in this MPI process so we
                // directly access the memory
                foreach (i; 0 .. this_vertices.length) {
                    this_vertices[i].pos[gtl].set(mapped_vertices[i].pos[gtl]);
                }
            }
        } else {
            // shared memory -- we can access the data directly
            foreach (i; 0 .. this_vertices.length) {
                this_vertices[i].pos[gtl].set(mapped_vertices[i].pos[gtl]);
            }
        }
    }
}
