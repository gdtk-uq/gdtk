module bc.boundary_vertex_full_face_copy;

import fvvertex;
import fluidblock;
import sfluidblock;
import globaldata;
import globalconfig;
import bc;
import geom;

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
            // TODO
        } else {
            // shared memory -- we can directly acces the memory in the last step
        }
    }

    void setup_vertex_mapping_phase1() {
        version(mpi_parallel) {
            // TODO
        } else {
            // shared memory -- we can directly acces the memory in the last step
        }
    }

    void setup_vertex_mapping_phase2() {
        version(mpi_parallel) {
            // TODO
        } else {
            // shared memory -- we can directly acces the memory
            foreach (i; 0 .. this_vertices.length) {
                mapped_vertices ~= other_blk.vertices[mapped_vertex_ids[i]];
            }

            // flag destination vertices not to be solved
            if (incoming) {
                foreach (v; this_vertices) v.solve_position = false;
            }
        }
    }

    void exchange_vertex_pos_phase0(int gtl) {
        version(mpi_parallel) {
            // TODO
        } else {
            // shared memory -- we can directly access the memory in the last step
        }
    }

    void exchange_vertex_pos_phase1(int gtl) {
        version(mpi_parallel) {
            // TODO
        } else {
            // shared memory -- we can directly access the memory in the last step
        }
    }

    void exchange_vertex_pos_phase2(int gtl) {
        version(mpi_parallel) {
            // TODO
        } else {
            // shared memory -- we can access the data directly
            if (incoming) {
                foreach (i; 0 .. this_vertices.length) {
                    this_vertices[i].pos[gtl].set(mapped_vertices[i].pos[gtl]);
                }
            }
        }
    }
}
