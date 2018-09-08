// full_face_copy.d

module bc.ghost_cell_effect.full_face_copy;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.algorithm;
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

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

// ----------------------------------------------------------------------------------
// MPI-specific services.

int make_mpi_tag(int blk_id, int bndry_id, int seq)
{
    // This function provides a somewhat unique integer
    // that can be used to encode message intent.
    assert(seq < 10, "mpi-oops, too many different messages");
    assert(bndry_id < 100, "mpi-oops, too many boundaries in a block");
    assert(blk_id < 100000, "mpi_oops, too many blocks in a simulation");
    return blk_id*1000 + bndry_id*10 + seq; // less than 2^^31
}

version(mpi_parallel) {
    int MPI_Wait_a_while(MPI_Request* request, MPI_Status *status)
    {
        int ierr = 0;
        int flag = 0;
        while (!flag) {
            ierr = MPI_Test(request, &flag, status);
        }
        return ierr;
    } // end MPI_Wait_a_while()
}

// ----------------------------------------------------------------------------------


class GhostCellFullFaceCopy : GhostCellEffect {
public:
    FluidBlock neighbourBlock;
    int neighbourFace;
    int neighbourOrientation;
    bool reorient_vector_quantities;
    double[] Rmatrix;
    // For each ghost cell associated with the boundary,
    // we will have a corresponding "mapped" or "source" cell
    // from which we will copy the flow conditions.
    FVCell[] ghost_cells;
    FVCell[] mapped_cells;
    size_t[] mapped_cell_ids;
    // Later, it is convenient to use a different notation for the data exchange.
    // Also, note that we require structured-grid blocks.
    SFluidBlock this_blk;
    SFluidBlock other_blk;
    int other_face;
    int other_orientation;
    version(mpi_parallel) {
        // This GhostCellEffect is somewhat symmetric in that for each ghost-cell
        // source-cell mapping, there should be a corresponding mapping over in
        // the other source block so these the cells in the current block
        // for which data should be sent to the source block.
        size_t[] outgoing_mapped_cell_ids;
        FVCell[] outgoing_mapped_cells;
        int other_blk_rank;
        int outgoing_cell_ids_tag, incoming_cell_ids_tag;
        MPI_Request incoming_cell_ids_request;
        MPI_Status incoming_cell_ids_status;
        int[] outgoing_cell_ids_buf, incoming_cell_ids_buf;
        int outgoing_geometry_tag, incoming_geometry_tag;
        MPI_Request incoming_geometry_request;
        MPI_Status incoming_geometry_status;
        double[] outgoing_geometry_buf, incoming_geometry_buf;
        int outgoing_flowstate_tag, incoming_flowstate_tag;
        MPI_Request incoming_flowstate_request;
        MPI_Status incoming_flowstate_status;
        double[] outgoing_flowstate_buf, incoming_flowstate_buf;
    }

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient,
         bool reorient_vector_quantities,
         ref const(double[]) Rmatrix)
    {
        super(id, boundary, "FullFaceCopy");
        neighbourBlock = globalFluidBlocks[otherBlock];
        neighbourFace = otherFace;
        neighbourOrientation = orient;
        this.reorient_vector_quantities = reorient_vector_quantities;
        this.Rmatrix = Rmatrix.dup();
    }

    override string toString() const
    { 
        string str = "FullFaceCopy(otherBlock=" ~ to!string(neighbourBlock.id) ~ 
            ", otherFace=" ~ to!string(neighbourFace) ~ 
            ", orient=" ~ to!string(neighbourOrientation) ~
            ", reorient_vector_quantities=" ~ to!string(reorient_vector_quantities) ~
            ", Rmatrix=[";
        foreach(i, v; Rmatrix) {
            str ~= to!string(v);
            str ~= (i < Rmatrix.length-1) ? ", " : "]";
        }
        str ~= ")";
        return str;
    }

    void set_up_cell_mapping_phase0()
    {
        // We call this function only after all blocks have been constructed.
        //
        // The following references and names will be handy for the data exchange
        // that occurs for this ghost-cell effect.
        this_blk = cast(SFluidBlock) blk;
        assert(this_blk, "Destination FlowBlock must be a structured-grid block.");
        other_blk = cast(SFluidBlock) neighbourBlock;
        assert(other_blk, "Source FlowBlock must be a structured-grid block.");
        other_face = neighbourFace;
        other_orientation = neighbourOrientation;
        version(mpi_parallel) {
            other_blk_rank = GlobalConfig.mpi_rank_for_block[other_blk.id];
        }
        //
        // For the source cells, we use indices into the hypothetical block of active cells. 
        size_t i_src, j_src, k_src;
        // For ghost-cell indices into the destination block, we use the raw indices into the
        // larger underlying block array that includes the surrounding layers of ghost cells.
        size_t i_dest, j_dest, k_dest;
        //
        if (blk.myConfig.dimensions == 2) {
            // Handle the 2D case separately.
            switch (which_boundary) {
            case Face.north:
                j_dest = this_blk.jmax;  // index of the north-most plane of active cells
                foreach (i; 0 .. this_blk.nicell) {
                    i_dest = i + this_blk.imin;
                    ghost_cells ~= this_blk.get_cell(i_dest,j_dest+1);
                    ghost_cells ~= this_blk.get_cell(i_dest,j_dest+2);
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = other_blk.nicell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = other_blk.njcell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // i loop
                break;
            case Face.east:
                i_dest = this_blk.imax;  // index of the east-most plane of active cells
                foreach (j; 0 .. this_blk.njcell) {
                    j_dest = j + this_blk.jmin;
                    ghost_cells ~= this_blk.get_cell(i_dest+1,j_dest);
                    ghost_cells ~= this_blk.get_cell(i_dest+2,j_dest);
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = other_blk.njcell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = other_blk.nicell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // j loop
                break;
            case Face.south:
                j_dest = this_blk.jmin;  // index of the south-most plane of active cells
                foreach (i; 0 .. this_blk.nicell) {
                    i_dest = i + this_blk.imin;
                    ghost_cells ~= this_blk.get_cell(i_dest,j_dest-1);
                    ghost_cells ~= this_blk.get_cell(i_dest,j_dest-2);
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = other_blk.njcell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = other_blk.nicell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // i loop
                break;
            case Face.west:
                i_dest = this_blk.imin;  // index of the west-most plane of active cells
                foreach (j; 0 .. this_blk.njcell) {
                    j_dest = j + this_blk.jmin;
                    ghost_cells ~= this_blk.get_cell(i_dest-1,j_dest);
                    ghost_cells ~= this_blk.get_cell(i_dest-2,j_dest);
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = other_blk.nicell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = other_blk.njcell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // j loop
                break;
            default:
                assert(false, "Incorrect boundary connection, which_boundary.");
            } // end switch which_boundary
        } else {
            // presume dimensions == 3
            // Continue on with 3D work...
            final switch (which_boundary) {
            case Face.north:
                j_dest = this_blk.jmax;  // index of the north-most plane of active cells
                foreach (i; 0 .. this_blk.nicell) {
                    i_dest = i + this_blk.imin;
                    foreach (k; 0 .. this_blk.nkcell) {
                        k_dest = k + this_blk.kmin;
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest+1,k_dest);
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest+2,k_dest);
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; k_src = k; break;
                            case 1: i_src = k; k_src = i; break;
                            case 2: i_src = i; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = other_blk.nicell - 1; 
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = other_blk.njcell - k - 1; k_src = i; break;
                            case 2: j_src = other_blk.njcell - i - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; k_src = i; break;
                            case 2: i_src = other_blk.nicell - i - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njcell - i - 1; k_src = k; break;
                            case 1: j_src = k; k_src = i; break;
                            case 2: j_src = i; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = other_blk.njcell - k - 1; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = other_blk.nkcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; j_src = i; break;
                            case 2: i_src = other_blk.nicell - i - 1; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = other_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; j_src = k; break;
                            case 1: i_src = k; j_src = i; break;
                            case 2: i_src = i; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; j_src = other_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (other_face)
                    } // k loop
                } // i loop
                break;
            case Face.east:
                i_dest = this_blk.imax;  // index of the east-most plane of active cells
                foreach (j; 0 .. this_blk.njcell) {
                    j_dest = j + this_blk.jmin;
                    foreach (k; 0 .. this_blk.nkcell) {
                        k_dest = k + this_blk.kmin;
                        ghost_cells ~= this_blk.get_cell(i_dest+1,j_dest,k_dest);
                        ghost_cells ~= this_blk.get_cell(i_dest+2,j_dest,k_dest);
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = k; k_src = other_blk.nkcell - j - 1; break;
                            case 2: i_src = other_blk.nicell - j - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; k_src = j;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = other_blk.nicell - 1; 
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njcell - j - 1; k_src = k; break;
                            case 1: j_src = other_blk.njcell - k - 1; k_src = other_blk.nkcell - j - 1; break;
                            case 2: j_src = j; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = j;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - j - 1; k_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; k_src = other_blk.nkcell - j - 1; break;
                            case 2: i_src = j; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = j;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (other_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = k; k_src = other_blk.nkcell - j - 1; break;
                            case 2: j_src = other_blk.njcell - j - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = other_blk.njcell - k - 1; k_src = j;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = other_blk.nkcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - j - 1; j_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; j_src = other_blk.njcell - j - 1; break;
                            case 2: i_src = j; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = j;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = k; j_src = other_blk.njcell - j - 1; break;
                            case 2: i_src = other_blk.nicell - j - 1; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; j_src = j;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (other_face)
                    } // k loop
                } // j loop
                break;
            case Face.south:
                j_dest = this_blk.jmin;  // index of the south-most plane of active cells
                foreach (i; 0 .. this_blk.nicell) {
                    i_dest = i + this_blk.imin;
                    foreach (k; 0 .. this_blk.nkcell) {
                        k_dest = k + this_blk.kmin;
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest-1,k_dest);
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest-2,k_dest);
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = k; k_src = other_blk.nkcell - i - 1; break;
                            case 2: i_src = other_blk.nicell - i - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = other_blk.nicell - 1; 
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njcell - i - 1; k_src = k; break;
                            case 1: j_src = other_blk.njcell - k - 1; k_src = other_blk.nkcell - i - 1; break;
                            case 2: j_src = i; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; k_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; k_src = other_blk.nkcell - i - 1; break;
                            case 2: i_src = i; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = k; k_src = other_blk.nkcell - i - 1; break;
                            case 2: j_src = other_blk.njcell - i - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = other_blk.njcell - k - 1; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = other_blk.nkcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; j_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; j_src = other_blk.njcell - i - 1; break;
                            case 2: i_src = i; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = k; j_src = other_blk.njcell - i - 1; break;
                            case 2: i_src = other_blk.nicell - i - 1; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; j_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (other_face)
                    } // k loop
                } // i loop
                break;
            case Face.west:
                i_dest = this_blk.imin;  // index of the west-most plane of active cells
                foreach (j; 0 .. this_blk.njcell) {
                    j_dest = j + this_blk.jmin;
                    foreach (k; 0 .. this_blk.nkcell) {
                        k_dest = k + this_blk.kmin;
                        ghost_cells ~= this_blk.get_cell(i_dest-1,j_dest,k_dest);
                        ghost_cells ~= this_blk.get_cell(i_dest-2,j_dest,k_dest);
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - j - 1; k_src = k; break;
                            case 1: i_src = k; k_src = j; break;
                            case 2: i_src = j; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; k_src = other_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = other_blk.nicell - 1; 
                            final switch (other_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = other_blk.njcell - k - 1; k_src = j; break;
                            case 2: j_src = other_blk.njcell - j - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = other_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; k_src = j; break;
                            case 2: i_src = other_blk.nicell - j - 1; k_src = other_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = other_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njcell - j - 1; k_src = k; break;
                            case 1: j_src = k; k_src = j; break;
                            case 2: j_src = j; k_src = other_blk.nkcell - k - 1; break;
                            case 3: j_src = other_blk.njcell - k - 1; k_src = other_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = other_blk.nkcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = other_blk.nicell - k - 1; j_src = j; break;
                            case 2: i_src = other_blk.nicell - j - 1; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = other_blk.njcell - j - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - j - 1; j_src = k; break;
                            case 1: i_src = k; j_src = j; break;
                            case 2: i_src = j; j_src = other_blk.njcell - k - 1; break;
                            case 3: i_src = other_blk.nicell - k - 1; j_src = other_blk.njcell - j - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (other_face)
                    } // k loop
                } // j loop
                break;
            case Face.top:
                k_dest = this_blk.kmax;  // index of the top-most plane of active cells
                foreach (j; 0 .. this_blk.njcell) {
                    j_dest = j + this_blk.jmin;
                    foreach (i; 0 .. this_blk.nicell) {
                        i_dest = i + this_blk.imin;
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest,k_dest+1);
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest,k_dest+2);
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = j; k_src = other_blk.nkcell - i - 1; break;
                            case 2: i_src = other_blk.nicell - i - 1; k_src = other_blk.nkcell - j - 1; break;
                            case 3: i_src = other_blk.nicell - j - 1; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = other_blk.nicell - 1; 
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njcell - i - 1; k_src = j; break;
                            case 1: j_src = other_blk.njcell - j - 1; k_src = other_blk.nkcell - i - 1; break;
                            case 2: j_src = i; k_src = other_blk.nkcell - j - 1; break;
                            case 3: j_src = j; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; k_src = j; break;
                            case 1: i_src = other_blk.nicell - j - 1; k_src = other_blk.nkcell - i - 1; break;
                            case 2: i_src = i; k_src = other_blk.nkcell - j - 1; break;
                            case 3: i_src = j; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = j; k_src = other_blk.nkcell - i - 1; break;
                            case 2: j_src = other_blk.njcell - i - 1; k_src = other_blk.nkcell - j - 1; break;
                            case 3: j_src = other_blk.njcell - j - 1; k_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = other_blk.nkcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; j_src = j; break;
                            case 1: i_src = other_blk.nicell - j - 1; j_src = other_blk.njcell - i - 1; break;
                            case 2: i_src = i; j_src = other_blk.njcell - j - 1; break;
                            case 3: i_src = j; j_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = j; j_src = other_blk.njcell - i - 1; break;
                            case 2: i_src = other_blk.nicell - i - 1; j_src = other_blk.njcell - j - 1; break;
                            case 3: i_src = other_blk.nicell - j - 1; j_src = i;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (other_face)
                    } // i loop
                } // j loop
                break;
            case Face.bottom:
                k_dest = this_blk.kmin;  // index of the bottom-most plane of active cells
                foreach (j; 0 .. this_blk.njcell) {
                    j_dest = j + this_blk.jmin;
                    foreach (i; 0 .. this_blk.nicell) {
                        i_dest = i + this_blk.imin;
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest,k_dest-1);
                        ghost_cells ~= this_blk.get_cell(i_dest,j_dest,k_dest-2);
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; k_src = j; break;
                            case 1: i_src = j; k_src = i; break;
                            case 2: i_src = i; k_src = other_blk.nkcell - j - 1; break;
                            case 3: i_src = other_blk.nicell - j - 1; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = other_blk.nicell - 1; 
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = other_blk.njcell - j - 1; k_src = i; break;
                            case 2: j_src = other_blk.njcell - i - 1; k_src = other_blk.nkcell - j - 1; break;
                            case 3: j_src = j; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = other_blk.nicell - j - 1; k_src = i; break;
                            case 2: i_src = other_blk.nicell - i - 1; k_src = other_blk.nkcell - j - 1; break;
                            case 3: i_src = j; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njcell - i - 1; k_src = j; break;
                            case 1: j_src = j; k_src = i; break;
                            case 2: j_src = i; k_src = other_blk.nkcell - j - 1; break;
                            case 3: j_src = other_blk.njcell - j - 1; k_src = other_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = other_blk.nkcell - 1; 
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = other_blk.nicell - j - 1; j_src = i; break;
                            case 2: i_src = other_blk.nicell - i - 1; j_src = other_blk.njcell - j - 1; break;
                            case 3: i_src = j; j_src = other_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nicell - i - 1; j_src = j; break;
                            case 1: i_src = j; j_src = i; break;
                            case 2: i_src = i; j_src = other_blk.njcell - j - 1; break;
                            case 3: i_src = other_blk.nicell - j - 1; j_src = other_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch other_face
                    } // i loop
                } // j loop
            } // end switch which_boundary
        } // end if dimensions == ...
        //
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // For this particular GhostCellEffect, we are expecting somewhat symmetric
                // interaction with the other MPI process.  This block has to receive
                // a list of mapped source cells from the other block and it has to send its own
                // list of mapped source cells from which it wants geometry and flow data.
                //
                size_t ne = ghost_cells.length;
                if (incoming_cell_ids_buf.length < ne) { incoming_cell_ids_buf.length = ne; }
                //
                // Post non-blocking receive for data that we expect to receive later
                // from the other_blk MPI process.
                incoming_cell_ids_tag = make_mpi_tag(other_blk.id, other_face, 0);
                MPI_Irecv(incoming_cell_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                          incoming_cell_ids_tag, MPI_COMM_WORLD, &incoming_cell_ids_request);
            } else {
                // The source block happens to be in this MPI process so
                // we know that we can just access the cell data directly
                // in the final phase.
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly,
            // in the final phase.
        }
    } // end set_up_cell_mapping_phase0()

    void set_up_cell_mapping_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // Blocking send to corresponding non-blocking receive that was posted
                // at in other_blk MPI process.
                size_t ne = ghost_cells.length;
                if (outgoing_cell_ids_buf.length < ne) { outgoing_cell_ids_buf.length = ne; }
                assert(ne == mapped_cell_ids.length, "oops, wrong length");
                outgoing_cell_ids_tag = make_mpi_tag(blk.id, which_boundary, 0);
                foreach (i; 0 .. ne) { outgoing_cell_ids_buf[i] = to!int(mapped_cell_ids[i]); }
                MPI_Send(outgoing_cell_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                         outgoing_cell_ids_tag, MPI_COMM_WORLD);
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly
                // in the final phase.
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end set_up_cell_mapping_send()

    void set_up_cell_mapping_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // Once complete, copy the data back into the local context.
                MPI_Wait_a_while(&incoming_cell_ids_request, &incoming_cell_ids_status);
                size_t ne = ghost_cells.length;
                outgoing_mapped_cell_ids.length = ne;
                foreach (i; 0 .. ne) { outgoing_mapped_cell_ids[i] = to!size_t(incoming_cell_ids_buf[i]); }
                outgoing_mapped_cells.length = 0;
                foreach (id; outgoing_mapped_cell_ids) { outgoing_mapped_cells ~= this_blk.cells[id]; }
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    mapped_cells ~= other_blk.cells[mapped_cell_ids[i]];
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                mapped_cells ~= other_blk.cells[mapped_cell_ids[i]];
            }
        }
    } // end set_up_cell_mapping_phase2()
    
    ref FVCell get_mapped_cell(size_t i)
    {
        if (i < mapped_cells.length) {
            return mapped_cells[i];
        } else {
            throw new FlowSolverException(format("Reference to requested mapped-cell[%d] is not available.", i));
        }
    }

    void exchange_geometry_phase0()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the geometry data via messages.
                //
                // Prepare to exchange geometry data for the boundary cells.
                // To match .copy_values_from(mapped_cells[i], CopyDataOption.grid) as defined in fvcell.d.
                //
                size_t ne = ghost_cells.length * (this_blk.myConfig.n_grid_time_levels * 5 + 4);
                if (incoming_geometry_buf.length < ne) { incoming_geometry_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the other_blk MPI process.
                incoming_geometry_tag = make_mpi_tag(other_blk.id, other_face, 1);
                MPI_Irecv(incoming_geometry_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_geometry_tag, MPI_COMM_WORLD, &incoming_geometry_request);
            } else {
                // The source block happens to be in this MPI process so
                // we know that we can just access the cell data directly
                // in the final phase.
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly,
            // in the final phase.
        }
    } // end exchange_geometry_phase0()

    void exchange_geometry_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // Blocking send of this block's geometry data
                // to the corresponding non-blocking receive that was posted
                // in the other MPI process.
                outgoing_geometry_tag = make_mpi_tag(blk.id, which_boundary, 1);
                size_t ne = ghost_cells.length * (this_blk.myConfig.n_grid_time_levels * 5 + 4);
                if (outgoing_geometry_buf.length < ne) { outgoing_geometry_buf.length = ne; }
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    foreach (j; 0 .. this_blk.myConfig.n_grid_time_levels) {
                        outgoing_geometry_buf[ii++] = c.pos[j].x;
                        outgoing_geometry_buf[ii++] = c.pos[j].y;
                        outgoing_geometry_buf[ii++] = c.pos[j].z;
                        outgoing_geometry_buf[ii++] = c.volume[j];
                        outgoing_geometry_buf[ii++] = c.areaxy[j];
                    }
                    outgoing_geometry_buf[ii++] = c.iLength;
                    outgoing_geometry_buf[ii++] = c.jLength;
                    outgoing_geometry_buf[ii++] = c.kLength;
                    outgoing_geometry_buf[ii++] = c.L_min;
                }
                MPI_Send(outgoing_geometry_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                         outgoing_geometry_tag, MPI_COMM_WORLD);
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly
                // in the final phase.
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly
            // in the final phase.
        }
    } // end exchange_geometry_phase1()

    void exchange_geometry_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                MPI_Wait_a_while(&incoming_geometry_request, &incoming_geometry_status);
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    foreach (j; 0 .. this_blk.myConfig.n_grid_time_levels) {
                        c.pos[j].refx = incoming_geometry_buf[ii++];
                        c.pos[j].refy = incoming_geometry_buf[ii++];
                        c.pos[j].refz = incoming_geometry_buf[ii++];
                        c.volume[j] = incoming_geometry_buf[ii++];
                        c.areaxy[j] = incoming_geometry_buf[ii++];
                    }
                    c.iLength = incoming_geometry_buf[ii++];
                    c.jLength = incoming_geometry_buf[ii++];
                    c.kLength = incoming_geometry_buf[ii++];
                    c.L_min = incoming_geometry_buf[ii++];
                }
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            }
        }
    } // end exchange_geometry_phase2()

    void exchange_flowstate_phase0(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // For this particular GhostCellEffect, we are expecting somewhat symmetric
                // interaction with the other MPI process.
                //
                // Exchange FlowState data for the boundary cells.
                // To match the function over in flowstate.d
                // void copy_values_from(in FlowState other)
                // and over in gas_state.d
                // @nogc void copy_values_from(ref const(GasState) other) 
                //
                auto gmodel = this_blk.myConfig.gmodel;
                size_t nspecies = gmodel.n_species();
                size_t nmodes = gmodel.n_modes();
                size_t ne = ghost_cells.length * (nmodes*3 + nspecies + 23);
                if (incoming_flowstate_buf.length < ne) { incoming_flowstate_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the other_blk MPI process.
                incoming_flowstate_tag = make_mpi_tag(other_blk.id, other_face, 0);
                MPI_Irecv(incoming_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_flowstate_tag, MPI_COMM_WORLD, &incoming_flowstate_request);
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly
                // in the final phase.
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly
            // in the final phase.
        }
        // Done with setting up all non-blocking reads for MPI.
    } // end exchange_flowstate_phase0()

    void exchange_flowstate_phase1(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // For this particular GhostCellEffect, we are expecting somewhat symmetric
                // interaction with the other MPI process.
                //
                // Exchange FlowState data for the boundary cells.
                // To match the function over in flowstate.d
                // void copy_values_from(in FlowState other)
                // {
                //     gas.copy_values_from(other.gas);
                //     vel.set(other.vel);
                //     B.set(other.B);
                //     psi = other.psi;
                //     divB = other.divB;
                //     tke = other.tke;
                //     omega = other.omega;
                //     mu_t = other.mu_t;
                //     k_t = other.k_t;
                //     S = other.S;
                // }
                // and over in gas_state.d
                // @nogc void copy_values_from(ref const(GasState) other) 
                // {
                //     rho = other.rho;
                //     p = other.p;
                //     T = other.T;
                //     u = other.u;
                //     p_e = other.p_e;
                //     a = other.a;
                //     foreach (i; 0 .. u_modes.length) { u_modes[i] = other.u_modes[i]; }
                //     foreach (i; 0 .. T_modes.length) { T_modes[i] = other.T_modes[i]; }
                //     mu = other.mu;
                //     k = other.k;
                //     foreach (i; 0 .. k_modes.length) { k_modes[i] = other.k_modes[i]; }
                //     sigma = other.sigma;
                //     foreach (i; 0 .. massf.length) { massf[i] = other.massf[i]; }
                //     quality = other.quality;
                // }
                //
                auto gmodel = this_blk.myConfig.gmodel;
                size_t nspecies = gmodel.n_species();
                size_t nmodes = gmodel.n_modes();
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // in the other MPI process.
                size_t ne = ghost_cells.length * (nmodes*3 + nspecies + 23);
                if (outgoing_flowstate_buf.length < ne) { outgoing_flowstate_buf.length = ne; }
                outgoing_flowstate_tag = make_mpi_tag(blk.id, which_boundary, 0);
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    FlowState fs = c.fs;
                    GasState gs = fs.gas;
                    outgoing_flowstate_buf[ii++] = gs.rho;
                    outgoing_flowstate_buf[ii++] = gs.p;
                    outgoing_flowstate_buf[ii++] = gs.T;
                    outgoing_flowstate_buf[ii++] = gs.u;
                    outgoing_flowstate_buf[ii++] = gs.p_e;
                    outgoing_flowstate_buf[ii++] = gs.a;
                    foreach (j; 0 .. nmodes) { outgoing_flowstate_buf[ii++] = gs.u_modes[j]; }
                    foreach (j; 0 .. nmodes) { outgoing_flowstate_buf[ii++] = gs.T_modes[j]; }
                    outgoing_flowstate_buf[ii++] = gs.mu;
                    outgoing_flowstate_buf[ii++] = gs.k;
                    foreach (j; 0 .. nmodes) { outgoing_flowstate_buf[ii++] = gs.k_modes[j]; }
                    outgoing_flowstate_buf[ii++] = gs.sigma;
                    foreach (j; 0 .. nspecies) { outgoing_flowstate_buf[ii++] = gs.massf[j]; }
                    outgoing_flowstate_buf[ii++] = gs.quality;
                    outgoing_flowstate_buf[ii++] = fs.vel.x;
                    outgoing_flowstate_buf[ii++] = fs.vel.y;
                    outgoing_flowstate_buf[ii++] = fs.vel.z;
                    outgoing_flowstate_buf[ii++] = fs.B.x;
                    outgoing_flowstate_buf[ii++] = fs.B.y;
                    outgoing_flowstate_buf[ii++] = fs.B.z;
                    outgoing_flowstate_buf[ii++] = fs.psi;
                    outgoing_flowstate_buf[ii++] = fs.divB;
                    outgoing_flowstate_buf[ii++] = fs.tke;
                    outgoing_flowstate_buf[ii++] = fs.omega;
                    outgoing_flowstate_buf[ii++] = fs.mu_t;
                    outgoing_flowstate_buf[ii++] = fs.k_t;
                    outgoing_flowstate_buf[ii++] = to!double(fs.S);
                }
                MPI_Send(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                         outgoing_flowstate_tag, MPI_COMM_WORLD);
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly
                // in the final phase.
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly
            // in the final phase.
        }
        // Done with copying from source cells.
    } // end exchange_flowstate_phase1()
    
    void exchange_flowstate_phase2(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, other_blk.id).empty) {
                // The source block is in another MPI process, go fetch the data via messages.
                //
                auto gmodel = this_blk.myConfig.gmodel;
                size_t nspecies = gmodel.n_species();
                size_t nmodes = gmodel.n_modes();
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                MPI_Wait_a_while(&incoming_flowstate_request, &incoming_flowstate_status);
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    FlowState fs = c.fs;
                    GasState gs = fs.gas;
                    gs.rho = incoming_flowstate_buf[ii++];
                    gs.p = incoming_flowstate_buf[ii++];
                    gs.T = incoming_flowstate_buf[ii++];
                    gs.u = incoming_flowstate_buf[ii++];
                    gs.p_e = incoming_flowstate_buf[ii++];
                    gs.a = incoming_flowstate_buf[ii++];
                    foreach (j; 0 .. nmodes) { gs.u_modes[j] = incoming_flowstate_buf[ii++]; }
                    foreach (j; 0 .. nmodes) { gs.T_modes[j] = incoming_flowstate_buf[ii++]; }
                    gs.mu = incoming_flowstate_buf[ii++];
                    gs.k = incoming_flowstate_buf[ii++];
                    foreach (j; 0 .. nmodes) { gs.k_modes[j] = incoming_flowstate_buf[ii++]; }
                    gs.sigma = incoming_flowstate_buf[ii++];
                    foreach (j; 0 .. nspecies) { gs.massf[j] = incoming_flowstate_buf[ii++]; }
                    gs.quality = incoming_flowstate_buf[ii++];
                    fs.vel.refx = incoming_flowstate_buf[ii++];
                    fs.vel.refy = incoming_flowstate_buf[ii++];
                    fs.vel.refz = incoming_flowstate_buf[ii++];
                    fs.B.refx = incoming_flowstate_buf[ii++];
                    fs.B.refy = incoming_flowstate_buf[ii++];
                    fs.B.refz = incoming_flowstate_buf[ii++];
                    fs.psi = incoming_flowstate_buf[ii++];
                    fs.divB = incoming_flowstate_buf[ii++];
                    fs.tke = incoming_flowstate_buf[ii++];
                    fs.omega = incoming_flowstate_buf[ii++];
                    fs.mu_t = incoming_flowstate_buf[ii++];
                    fs.k_t = incoming_flowstate_buf[ii++];
                    fs.S = to!int(incoming_flowstate_buf[ii++]);
                }
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
            }
        }
        // Done with copying from source cells.
    } // end exchange_flowstate_phase2()

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("GhostCellFullFaceCopy.apply_unstructured_grid() not implemented");
    }
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        // We presume that all of the exchange of data happened earlier,
        // and that the ghost cells have been filled with flow state data
        // from their respective source cells.
        foreach (i; 0 .. ghost_cells.length) {
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
            // The following call to encode_conserved is needed
            // for the block-marching process.
            ghost_cells[i].encode_conserved(gtl, ftl, blk.omegaz);
        }
    } // end apply_structured_grid()
    
} // end class GhostCellFullFaceCopy
