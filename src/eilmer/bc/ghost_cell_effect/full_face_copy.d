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
    assert(seq < 10, "mpi-oops, too many different messages");
    assert(bndry_id < 100, "mpi-oops, too many boundaries in a block");
    assert(blk_id < 100000, "mpi_oops, too many blocks in a simulation");
    return blk_id*1000 + bndry_id*10 + seq; // less than 2^^31
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
    version(mpi_parallel) {
        // This GhostCellEffect is somewhat symmetric in that for each ghost-cell
        // source-cell mapping, there should be a corresponding mapping over in
        // the other source block so these the cells in the current block
        // for which data should be sent to the source block.
        size_t[] outgoing_mapped_cell_ids;
        FVCell[] outgoing_mapped_cells;
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

    void set_up_cell_mapping()
    {
        auto dest_blk = cast(SFluidBlock) blk;
        assert(dest_blk, "Destination FlowBlock must be a structured-grid block.");
        int destination_face = which_boundary;
        auto src_blk = cast(SFluidBlock) neighbourBlock;
        assert(src_blk, "Source FlowBlock must be a structured-grid block.");
        int src_face = neighbourFace;
        int src_orientation = neighbourOrientation;
        //
        // For the source cells, we use indices into the hypothetical block of active cells. 
        size_t i_src, j_src, k_src;
        // For ghost-cell indices into the destination block, we use the raw indices into the
        // larger underlying block array that includes the surrounding layers of ghost cells.
        size_t i_dest, j_dest, k_dest;
        //
        if (blk.myConfig.dimensions == 2) {
            // Handle the 2D case separately.
            switch (destination_face) {
            case Face.north:
                j_dest = dest_blk.jmax;  // index of the north-most plane of active cells
                foreach (i; 0 .. dest_blk.nicell) {
                    i_dest = i + dest_blk.imin;
                    ghost_cells ~= dest_blk.get_cell(i_dest,j_dest+1);
                    ghost_cells ~= dest_blk.get_cell(i_dest,j_dest+2);
                    switch (src_face) {
                    case Face.north:
                        j_src = src_blk.njcell - 1; 
                        i_src = src_blk.nicell - i - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = src_blk.nicell - 1; 
                        j_src = i;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = i;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = src_blk.njcell - i - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                } // i loop
                break;
            case Face.east:
                i_dest = dest_blk.imax;  // index of the east-most plane of active cells
                foreach (j; 0 .. dest_blk.njcell) {
                    j_dest = j + dest_blk.jmin;
                    ghost_cells ~= dest_blk.get_cell(i_dest+1,j_dest);
                    ghost_cells ~= dest_blk.get_cell(i_dest+2,j_dest);
                    switch (src_face) {
                    case Face.north:
                        j_src = src_blk.njcell - 1; 
                        i_src = j;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = src_blk.nicell - 1; 
                        j_src = src_blk.njcell - j - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = src_blk.nicell - j - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = j;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                } // j loop
                break;
            case Face.south:
                j_dest = dest_blk.jmin;  // index of the south-most plane of active cells
                foreach (i; 0 .. dest_blk.nicell) {
                    i_dest = i + dest_blk.imin;
                    ghost_cells ~= dest_blk.get_cell(i_dest,j_dest-1);
                    ghost_cells ~= dest_blk.get_cell(i_dest,j_dest-2);
                    switch (src_face) {
                    case Face.north:
                        j_src = src_blk.njcell - 1; 
                        i_src = i;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = src_blk.nicell - 1; 
                        j_src = src_blk.njcell - i - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = src_blk.nicell - i - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = i;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                } // i loop
                break;
            case Face.west:
                i_dest = dest_blk.imin;  // index of the west-most plane of active cells
                foreach (j; 0 .. dest_blk.njcell) {
                    j_dest = j + dest_blk.jmin;
                    ghost_cells ~= dest_blk.get_cell(i_dest-1,j_dest);
                    ghost_cells ~= dest_blk.get_cell(i_dest-2,j_dest);
                    switch (src_face) {
                    case Face.north:
                        j_src = src_blk.njcell - 1; 
                        i_src = src_blk.nicell - j - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        i_src = src_blk.nicell - 1; 
                        j_src = j;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = j;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = src_blk.njcell - j - 1;
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                } // j loop
                break;
            default:
                assert(false, "Incorrect boundary connection, destination face.");
            } // end switch destination_face
        } else {
            // presume dimensions == 3
            // Continue on with 3D work...
            final switch (destination_face) {
            case Face.north:
                j_dest = dest_blk.jmax;  // index of the north-most plane of active cells
                foreach (i; 0 .. dest_blk.nicell) {
                    i_dest = i + dest_blk.imin;
                    foreach (k; 0 .. dest_blk.nkcell) {
                        k_dest = k + dest_blk.kmin;
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest+1,k_dest);
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest+2,k_dest);
                        final switch (src_face) {
                        case Face.north:
                            j_src = src_blk.njcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = k; break;
                            case 1: i_src = k; k_src = i; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = src_blk.nicell - 1; 
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = i; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = k; break;
                            case 1: j_src = k; k_src = i; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = src_blk.nkcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = src_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = k; break;
                            case 1: i_src = k; j_src = i; break;
                            case 2: i_src = i; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                    } // k loop
                } // i loop
                break;
            case Face.east:
                i_dest = dest_blk.imax;  // index of the east-most plane of active cells
                foreach (j; 0 .. dest_blk.njcell) {
                    j_dest = j + dest_blk.jmin;
                    foreach (k; 0 .. dest_blk.nkcell) {
                        k_dest = k + dest_blk.kmin;
                        ghost_cells ~= dest_blk.get_cell(i_dest+1,j_dest,k_dest);
                        ghost_cells ~= dest_blk.get_cell(i_dest+2,j_dest,k_dest);
                        final switch (src_face) {
                        case Face.north:
                            j_src = src_blk.njcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = k; k_src = src_blk.nkcell - j - 1; break;
                            case 2: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = j;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = src_blk.nicell - 1; 
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - j - 1; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 2: j_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = j;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 2: i_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = j;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (src_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = k; k_src = src_blk.nkcell - j - 1; break;
                            case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = j;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = src_blk.nkcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - j - 1; break;
                            case 2: i_src = j; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = j;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = k; j_src = src_blk.njcell - j - 1; break;
                            case 2: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = j;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                    } // k loop
                } // j loop
                break;
            case Face.south:
                j_dest = dest_blk.jmin;  // index of the south-most plane of active cells
                foreach (i; 0 .. dest_blk.nicell) {
                    i_dest = i + dest_blk.imin;
                    foreach (k; 0 .. dest_blk.nkcell) {
                        k_dest = k + dest_blk.kmin;
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest-1,k_dest);
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest-2,k_dest);
                        final switch (src_face) {
                        case Face.north:
                            j_src = src_blk.njcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = k; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = src_blk.nicell - 1; 
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = k; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = src_blk.nkcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = i; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = k; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                    } // k loop
                } // i loop
                break;
            case Face.west:
                i_dest = dest_blk.imin;  // index of the west-most plane of active cells
                foreach (j; 0 .. dest_blk.njcell) {
                    j_dest = j + dest_blk.jmin;
                    foreach (k; 0 .. dest_blk.nkcell) {
                        k_dest = k + dest_blk.kmin;
                        ghost_cells ~= dest_blk.get_cell(i_dest-1,j_dest,k_dest);
                        ghost_cells ~= dest_blk.get_cell(i_dest-2,j_dest,k_dest);
                        final switch (src_face) {
                        case Face.north:
                            j_src = src_blk.njcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; k_src = k; break;
                            case 1: i_src = k; k_src = j; break;
                            case 2: i_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = src_blk.nicell - 1; 
                            final switch (src_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = j; break;
                            case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = src_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = j; break;
                            case 2: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = src_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - j - 1; k_src = k; break;
                            case 1: j_src = k; k_src = j; break;
                            case 2: j_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - j - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = src_blk.nkcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = j; break;
                            case 2: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = src_blk.njcell - j - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; j_src = k; break;
                            case 1: i_src = k; j_src = j; break;
                            case 2: i_src = j; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - j - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                    } // k loop
                } // j loop
                break;
            case Face.top:
                k_dest = dest_blk.kmax;  // index of the top-most plane of active cells
                foreach (j; 0 .. dest_blk.njcell) {
                    j_dest = j + dest_blk.jmin;
                    foreach (i; 0 .. dest_blk.nicell) {
                        i_dest = i + dest_blk.imin;
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest,k_dest+1);
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest,k_dest+2);
                        final switch (src_face) {
                        case Face.north:
                            j_src = src_blk.njcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = j; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = src_blk.nicell - 1; 
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = j; break;
                            case 1: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = j; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = j; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = j; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = src_blk.njcell - j - 1; k_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = src_blk.nkcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = i; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = j; j_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = j; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; j_src = i;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                    } // i loop
                } // j loop
                break;
            case Face.bottom:
                k_dest = dest_blk.kmin;  // index of the bottom-most plane of active cells
                foreach (j; 0 .. dest_blk.njcell) {
                    j_dest = j + dest_blk.jmin;
                    foreach (i; 0 .. dest_blk.nicell) {
                        i_dest = i + dest_blk.imin;
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest,k_dest-1);
                        ghost_cells ~= dest_blk.get_cell(i_dest,j_dest,k_dest-2);
                        final switch (src_face) {
                        case Face.north:
                            j_src = src_blk.njcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
                            case 1: i_src = j; k_src = i; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            i_src = src_blk.nicell - 1; 
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = src_blk.njcell - j - 1; k_src = i; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = j; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            j_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; k_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = j; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            i_src = 0; 
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = j; break;
                            case 1: j_src = j; k_src = i; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            k_src = src_blk.nkcell - 1; 
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; j_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = j; j_src = src_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            k_src = 0; 
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = j; break;
                            case 1: i_src = j; j_src = i; break;
                            case 2: i_src = i; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - i - 1;
                            }
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                            mapped_cell_ids ~= src_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch src_face
                    } // i loop
                } // j loop
            } // end switch destination_face
        } // end if dimensions == ...
        //
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, src_blk.id).empty) {
                // The source block is in another MPI process, go fetch the data via messages.
                //
                // For this particular GhostCellEffect, we are expecting somewhat symmetric
                // interaction with the other MPI process.  This block has to receive
                // a list of mapped source cells from the other block and it has to send its own
                // list of mapped source cells from which it wants geometry and flow data.
                //
                size_t ne = ghost_cells.length;
                if (incoming_int_buf.length < ne) { incoming_int_buf.length = ne; }
                if (outgoing_int_buf.length < ne) { outgoing_int_buf.length = ne; }
                //
                // Post non-blocking receive for data that we expect to receive later
                // from the src_blk MPI process.
                int tag = make_mpi_tag(src_blk.id, src_face, 0);
                int src_blk_rank = GlobalConfig.mpi_rank_for_block[src_blk.id];
                int dest_blk_rank = GlobalConfig.mpi_rank_for_local_task;
                MPI_Request my_request;
                MPI_Status my_status;
                MPI_Irecv(incoming_int_buf.ptr, to!int(ne), MPI_INT, src_blk_rank, tag,
                          MPI_COMM_WORLD, &my_request);
                //
                // Blocking send to corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                assert(ne == mapped_cell_ids.length, "oops, wrong length");
                tag = make_mpi_tag(dest_blk.id, destination_face, 0);
                foreach (i; 0 .. ne) { outgoing_int_buf[i] = to!int(mapped_cell_ids[i]); }
                MPI_Send(outgoing_int_buf.ptr, to!int(ne), MPI_INT, src_blk_rank, tag, MPI_COMM_WORLD);
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
		MPI_Wait(&my_request, &my_status);
                outgoing_mapped_cell_ids.length = ne;
                foreach (i; 0 .. ne) { outgoing_mapped_cell_ids[i] = to!size_t(incoming_int_buf[i]); }
                outgoing_mapped_cells.length = 0;
                foreach (i; outgoing_mapped_cell_ids) { outgoing_mapped_cells ~= dest_blk.cells[i]; }
                //
                // Exchange geometry data for the boundary cells.
                // To match .copy_values_from(mapped_cells[i], CopyDataOption.grid) as defined in fvcell.d.
                //
                ne = ghost_cells.length * (dest_blk.myConfig.n_grid_time_levels * 5 + 4);
                if (incoming_double_buf.length < ne) { incoming_double_buf.length = ne; }
                if (outgoing_double_buf.length < ne) { outgoing_double_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the src_blk MPI process.
                tag = make_mpi_tag(src_blk.id, src_face, 1);
                MPI_Irecv(incoming_double_buf.ptr, to!int(ne), MPI_DOUBLE, src_blk_rank, tag,
                          MPI_COMM_WORLD, &my_request);
                // Blocking send of this block's geometry data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                tag = make_mpi_tag(dest_blk.id, destination_face, 1);
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    foreach (j; 0 .. dest_blk.myConfig.n_grid_time_levels) {
                        outgoing_double_buf[ii++] = c.pos[j].x;
                        outgoing_double_buf[ii++] = c.pos[j].y;
                        outgoing_double_buf[ii++] = c.pos[j].z;
                        outgoing_double_buf[ii++] = c.volume[j];
                        outgoing_double_buf[ii++] = c.areaxy[j];
                    }
                    outgoing_double_buf[ii++] = c.iLength;
                    outgoing_double_buf[ii++] = c.jLength;
                    outgoing_double_buf[ii++] = c.kLength;
                    outgoing_double_buf[ii++] = c.L_min;
                }
                MPI_Send(outgoing_double_buf.ptr, to!int(ne), MPI_DOUBLE, src_blk_rank, tag, MPI_COMM_WORLD);
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
		MPI_Wait(&my_request, &my_status);
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                ii = 0;
                foreach (c; ghost_cells) {
                    foreach (j; 0 .. dest_blk.myConfig.n_grid_time_levels) {
                        c.pos[j].refx = incoming_double_buf[ii++];
                        c.pos[j].refy = incoming_double_buf[ii++];
                        c.pos[j].refz = incoming_double_buf[ii++];
                        c.volume[j] = incoming_double_buf[ii++];
                        c.areaxy[j] = incoming_double_buf[ii++];
                    }
                    c.iLength = incoming_double_buf[ii++];
                    c.jLength = incoming_double_buf[ii++];
                    c.kLength = incoming_double_buf[ii++];
                    c.L_min = incoming_double_buf[ii++];
                }
            } else {
                // The source block happens to be in this MPI process so
                // we know that we can just access the source-cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    mapped_cells ~= src_blk.cells[mapped_cell_ids[i]];
                }
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                mapped_cells ~= src_blk.cells[mapped_cell_ids[i]];
            }
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            }
        }
    } // end set_up_cell_mapping()
    
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
        throw new Error("GhostCellFullFaceCopy.apply_unstructured_grid() not implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto dest_blk = cast(SFluidBlock) blk;
        assert(dest_blk, "Destination FlowBlock must be a structured-grid block.");
        int destination_face = which_boundary;
        auto src_blk = cast(SFluidBlock) neighbourBlock;
        assert(src_blk, "Source FlowBlock must be a structured-grid block.");
        int src_face = neighbourFace;
        version(mpi_parallel) {
            if (find(GlobalConfig.localBlockIds, src_blk.id).empty) {
                // The source block is in another MPI process, go fetch the data via messages.
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
                auto gmodel = dest_blk.myConfig.gmodel;
                size_t nspecies = gmodel.n_species();
                size_t nmodes = gmodel.n_modes();
                size_t ne = ghost_cells.length * (nmodes*3 + nspecies + 23);
                if (incoming_double_buf.length < ne) { incoming_double_buf.length = ne; }
                if (outgoing_double_buf.length < ne) { outgoing_double_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the src_blk MPI process.
                int tag = make_mpi_tag(src_blk.id, src_face, 0);
                int src_blk_rank = GlobalConfig.mpi_rank_for_block[src_blk.id];
                int dest_blk_rank = GlobalConfig.mpi_rank_for_local_task;
                MPI_Request my_request;
                MPI_Status my_status;
                MPI_Irecv(incoming_double_buf.ptr, to!int(ne), MPI_DOUBLE, src_blk_rank, tag,
                          MPI_COMM_WORLD, &my_request);
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                tag = make_mpi_tag(dest_blk.id, destination_face, 0);
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    FlowState fs = c.fs;
                    GasState gs = fs.gas;
                    outgoing_double_buf[ii++] = gs.rho;
                    outgoing_double_buf[ii++] = gs.p;
                    outgoing_double_buf[ii++] = gs.T;
                    outgoing_double_buf[ii++] = gs.u;
                    outgoing_double_buf[ii++] = gs.p_e;
                    outgoing_double_buf[ii++] = gs.a;
                    foreach (j; 0 .. nmodes) { outgoing_double_buf[ii++] = gs.u_modes[j]; }
                    foreach (j; 0 .. nmodes) { outgoing_double_buf[ii++] = gs.T_modes[j]; }
                    outgoing_double_buf[ii++] = gs.mu;
                    outgoing_double_buf[ii++] = gs.k;
                    foreach (j; 0 .. nmodes) { outgoing_double_buf[ii++] = gs.k_modes[j]; }
                    outgoing_double_buf[ii++] = gs.sigma;
                    foreach (j; 0 .. nspecies) { outgoing_double_buf[ii++] = gs.massf[j]; }
                    outgoing_double_buf[ii++] = gs.quality;
                    outgoing_double_buf[ii++] = fs.vel.x;
                    outgoing_double_buf[ii++] = fs.vel.y;
                    outgoing_double_buf[ii++] = fs.vel.z;
                    outgoing_double_buf[ii++] = fs.B.x;
                    outgoing_double_buf[ii++] = fs.B.y;
                    outgoing_double_buf[ii++] = fs.B.z;
                    outgoing_double_buf[ii++] = fs.psi;
                    outgoing_double_buf[ii++] = fs.divB;
                    outgoing_double_buf[ii++] = fs.tke;
                    outgoing_double_buf[ii++] = fs.omega;
                    outgoing_double_buf[ii++] = fs.mu_t;
                    outgoing_double_buf[ii++] = fs.k_t;
                    outgoing_double_buf[ii++] = to!double(fs.S);
                }
                MPI_Send(outgoing_double_buf.ptr, to!int(ne), MPI_DOUBLE, src_blk_rank, tag, MPI_COMM_WORLD);
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
		MPI_Wait(&my_request, &my_status);
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                ii = 0;
                foreach (c; ghost_cells) {
                    FlowState fs = c.fs;
                    GasState gs = fs.gas;
                    gs.rho = incoming_double_buf[ii++];
                    gs.p = incoming_double_buf[ii++];
                    gs.T = incoming_double_buf[ii++];
                    gs.u = incoming_double_buf[ii++];
                    gs.p_e = incoming_double_buf[ii++];
                    gs.a = incoming_double_buf[ii++];
                    foreach (j; 0 .. nmodes) { gs.u_modes[j] = incoming_double_buf[ii++]; }
                    foreach (j; 0 .. nmodes) { gs.T_modes[j] = incoming_double_buf[ii++]; }
                    gs.mu = incoming_double_buf[ii++];
                    gs.k = incoming_double_buf[ii++];
                    foreach (j; 0 .. nmodes) { gs.k_modes[j] = incoming_double_buf[ii++]; }
                    gs.sigma = incoming_double_buf[ii++];
                    foreach (j; 0 .. nspecies) { gs.massf[j] = incoming_double_buf[ii++]; }
                    gs.quality = incoming_double_buf[ii++];
                    fs.vel.refx = incoming_double_buf[ii++];
                    fs.vel.refy = incoming_double_buf[ii++];
                    fs.vel.refz = incoming_double_buf[ii++];
                    fs.B.refx = incoming_double_buf[ii++];
                    fs.B.refy = incoming_double_buf[ii++];
                    fs.B.refz = incoming_double_buf[ii++];
                    fs.psi = incoming_double_buf[ii++];
                    fs.divB = incoming_double_buf[ii++];
                    fs.tke = incoming_double_buf[ii++];
                    fs.omega = incoming_double_buf[ii++];
                    fs.mu_t = incoming_double_buf[ii++];
                    fs.k_t = incoming_double_buf[ii++];
                    fs.S = to!int(incoming_double_buf[ii++]);
                }
                //
                // Exchange geometry data for the boundary cells.
                // To match .copy_values_from(mapped_cells[i], CopyDataOption.grid) as defined in fvcell.d.
                //
                ne = ghost_cells.length * (dest_blk.myConfig.n_grid_time_levels * 5 + 4);
                if (incoming_double_buf.length < ne) { incoming_double_buf.length = ne; }
                if (outgoing_double_buf.length < ne) { outgoing_double_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the src_blk MPI process.
                tag = make_mpi_tag(src_blk.id, src_face, 1);
                MPI_Irecv(incoming_double_buf.ptr, to!int(ne), MPI_DOUBLE, src_blk_rank, tag,
                          MPI_COMM_WORLD, &my_request);
                // Blocking send of this block's geometry data
                // to the corresponding non-blocking receive that was posted
                // at in src_blk MPI process.
                tag = make_mpi_tag(dest_blk.id, destination_face, 1);
                ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    foreach (j; 0 .. dest_blk.myConfig.n_grid_time_levels) {
                        outgoing_double_buf[ii++] = c.pos[j].x;
                        outgoing_double_buf[ii++] = c.pos[j].y;
                        outgoing_double_buf[ii++] = c.pos[j].z;
                        outgoing_double_buf[ii++] = c.volume[j];
                        outgoing_double_buf[ii++] = c.areaxy[j];
                    }
                    outgoing_double_buf[ii++] = c.iLength;
                    outgoing_double_buf[ii++] = c.jLength;
                    outgoing_double_buf[ii++] = c.kLength;
                    outgoing_double_buf[ii++] = c.L_min;
                }
                MPI_Send(outgoing_double_buf.ptr, to!int(ne), MPI_DOUBLE, src_blk_rank, tag, MPI_COMM_WORLD);
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
		MPI_Wait(&my_request, &my_status);
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                ii = 0;
                foreach (c; ghost_cells) {
                    foreach (j; 0 .. dest_blk.myConfig.n_grid_time_levels) {
                        c.pos[j].refx = incoming_double_buf[ii++];
                        c.pos[j].refy = incoming_double_buf[ii++];
                        c.pos[j].refz = incoming_double_buf[ii++];
                        c.volume[j] = incoming_double_buf[ii++];
                        c.areaxy[j] = incoming_double_buf[ii++];
                    }
                    c.iLength = incoming_double_buf[ii++];
                    c.jLength = incoming_double_buf[ii++];
                    c.kLength = incoming_double_buf[ii++];
                    c.L_min = incoming_double_buf[ii++];
                }
            } else {
                // The source block happens to be in this MPI process so
                // we know that we can just access the source-cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
                    ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
                ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            }
        }
        // Done with copying from source cells.
        foreach (i; 0 .. ghost_cells.length) {
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
            // The following call to encode_conserved is needed
            // for the block-marching process.
            ghost_cells[i].encode_conserved(gtl, ftl, blk.omegaz);
        }
    }
} // end class GhostCellFullFaceCopy
