// solid_full_face_copy.d

module solid_full_face_copy;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.algorithm;
import std.datetime;
version(mpi_parallel) {
    import mpi;
}

import simcore;
import util.json_helper;
import geom;
import globaldata;
import globalconfig;
import solidfvinterface;
import ssolidblock;
import solidfvcell;
import solidbc;
import solid_ghost_cell;

// ----------------------------------------------------------------------------------
// MPI-specific services.

@nogc
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
    // not @nogc
    int MPI_Wait_a_while(MPI_Request* request, MPI_Status *status)
    {
        int ierr = 0; // Normally we return the MPI_Test result.
        long timeOut_msecs = 10000; // Surely 10 seconds will be enough.
        SysTime startTime = Clock.currTime();
        int flag = 0;
        while (!flag) {
            ierr = MPI_Test(request, &flag, status);
            long elapsedTime_msecs = (Clock.currTime() - startTime).total!"msecs"();
            if (elapsedTime_msecs > timeOut_msecs) {
                // [TODO] Maybe we should consider cleaning up MPI resources, however,
                // we do not expect our job to recover gracefully from this point.
                throw new Exception("MPI_wait_a_while time-out has occurred.");
            }
        }
        return ierr;
    } // end MPI_Wait_a_while()
}

// ----------------------------------------------------------------------------------


class SolidGCE_SolidGhostCellFullFaceCopy : SolidGhostCellEffect {
public:
    SSolidBlock neighbourBlock;
    int neighbourFace;
    int neighbourOrientation;
    // For each ghost cell associated with the boundary,
    // we will have a corresponding "mapped" or "source" cell
    // from which we will copy the flow conditions.
    SolidFVCell[] ghost_cells;
    SolidFVCell[] mapped_cells;
    size_t[] mapped_cell_ids;
    // Later, it is convenient to use a different notation for the data exchange.
    // Also, note that we require structured-grid blocks.
    SSolidBlock this_blk;
    SSolidBlock other_blk;
    int other_face;
    int other_orientation;
    version(mpi_parallel) {
        // This GhostCellEffect is somewhat symmetric in that for each ghost-cell
        // source-cell mapping, there should be a corresponding mapping over in
        // the other source block so these the cells in the current block
        // for which data should be sent to the source block.
        size_t[] outgoing_mapped_cell_ids;
        SolidFVCell[] outgoing_mapped_cells;
        int other_blk_rank;
        int outgoing_cell_ids_tag, incoming_cell_ids_tag;
        MPI_Request incoming_cell_ids_request;
        MPI_Status incoming_cell_ids_status;
        int[] outgoing_cell_ids_buf, incoming_cell_ids_buf;
        int outgoing_solidstate_tag, incoming_solidstate_tag;
        MPI_Request incoming_solidstate_request;
        MPI_Status incoming_solidstate_status;
        double[] outgoing_solidstate_buf, incoming_solidstate_buf;
    }

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient)
    {
        super(id, boundary, "SolidFullFaceCopy");
        neighbourBlock = cast(SSolidBlock) globalBlocks[otherBlock];
        assert(blk !is null, "Oops, this should be a SSolidBlock object.");
        neighbourFace = otherFace;
        neighbourOrientation = orient;
    }

    override string toString() const
    { 
        string str = "FullFaceCopy(otherBlock=" ~ to!string(neighbourBlock.id) ~ 
            ", otherFace=" ~ to!string(neighbourFace) ~ 
            ", orient=" ~ to!string(neighbourOrientation); 
        str ~= ")";
        return str;
    }

    void set_up_cell_mapping_phase0()
    {
        // We call this function only after all blocks have been constructed.
        //
        // The following references and names will be handy for the data exchange
        // that occurs for this ghost-cell effect.
        this_blk = cast(SSolidBlock) blk;
        if (!this_blk) { throw new Error("Destination SolidBlock must be a structured-grid block."); }
        bool nghost3 = (this_blk.n_ghost_cell_layers == 1);
        other_blk = cast(SSolidBlock) neighbourBlock;
        if (!other_blk) { throw new Error("Source SolidBlock must be a structured-grid block."); }
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
                    ghost_cells ~= this_blk.getCell(i_dest,j_dest+1);
                    ghost_cells ~= this_blk.getCell(i_dest,j_dest+2);
                    if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest,j_dest+3); }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = other_blk.nicell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2); }
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src); }
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2); }
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = other_blk.njcell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src,k_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                    ghost_cells ~= this_blk.getCell(i_dest+1,j_dest);
                    ghost_cells ~= this_blk.getCell(i_dest+2,j_dest);
                    if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest+3,j_dest); }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2); }
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = other_blk.njcell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src); }
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = other_blk.nicell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2); }
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src); }
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
                    ghost_cells ~= this_blk.getCell(i_dest,j_dest-1);
                    ghost_cells ~= this_blk.getCell(i_dest,j_dest-2);
                    if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest,j_dest-3); }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2); }
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = other_blk.njcell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src); }
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = other_blk.nicell - i - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2); }
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = i;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src); }
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
                    ghost_cells ~= this_blk.getCell(i_dest-1,j_dest);
                    ghost_cells ~= this_blk.getCell(i_dest-2,j_dest);
                    if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest-3,j_dest); }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njcell - 1; 
                        i_src = other_blk.nicell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2); }
                        break;
                    case Face.east:
                        i_src = other_blk.nicell - 1; 
                        j_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src); }
                        break;
                    case Face.south:
                        j_src = 0; 
                        i_src = j;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+1);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2); }
                        break;
                    case Face.west:
                        i_src = 0; 
                        j_src = other_blk.njcell - j - 1;
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src);
                        mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+1,j_src);
                        if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src); }
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
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest+1,k_dest);
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest+2,k_dest);
                        if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest,j_dest+3,k_dest); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-2); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+2); }
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
                        ghost_cells ~= this_blk.getCell(i_dest+1,j_dest,k_dest);
                        ghost_cells ~= this_blk.getCell(i_dest+2,j_dest,k_dest);
                        if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest+3,j_dest,k_dest); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-2); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+2); }
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
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest-1,k_dest);
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest-2,k_dest);
                        if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest,j_dest-3,k_dest); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-2); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+2); }
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
                        ghost_cells ~= this_blk.getCell(i_dest-1,j_dest,k_dest);
                        ghost_cells ~= this_blk.getCell(i_dest-2,j_dest,k_dest);
                        if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest-3,j_dest,k_dest); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-2); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+2); }
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
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest,k_dest+1);
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest,k_dest+2);
                        if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest,j_dest,k_dest+3); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-2); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+2); }
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
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest,k_dest-1);
                        ghost_cells ~= this_blk.getCell(i_dest,j_dest,k_dest-2);
                        if (nghost3) { ghost_cells ~= this_blk.getCell(i_dest,j_dest,k_dest-3); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src-2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src-2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src+2,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src+2,j_src,k_src); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src-2); }
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
                            if (nghost3) { mapped_cell_ids ~= other_blk.ijk_0n_indices_to_cell_id(i_src,j_src,k_src+2); }
                        } // end switch other_face
                    } // i loop
                } // j loop
            } // end switch which_boundary
        } // end if dimensions == ...
        //
        version(mpi_parallel) {
            if (find(GlobalConfig.localSolidBlockIds, other_blk.id).empty) {
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
                incoming_cell_ids_tag = make_mpi_tag(other_blk.id, other_face, 5);
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

    // not @nogc because we set array length
    void set_up_cell_mapping_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localSolidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // Blocking send to corresponding non-blocking receive that was posted
                // at in other_blk MPI process.
                size_t ne = ghost_cells.length;
                if (outgoing_cell_ids_buf.length < ne) { outgoing_cell_ids_buf.length = ne; }
                assert(ne == mapped_cell_ids.length, "oops, wrong length");
                outgoing_cell_ids_tag = make_mpi_tag(blk.id, which_boundary, 5);
                foreach (i; 0 .. ne) { outgoing_cell_ids_buf[i] = to!int(mapped_cell_ids[i]); }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_cell_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                              outgoing_cell_ids_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_cell_ids_buf.ptr, to!int(ne), MPI_INT, other_blk_rank,
                             outgoing_cell_ids_tag, MPI_COMM_WORLD);
                }
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

    // not @nogc
    void set_up_cell_mapping_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localSolidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_cell_ids_request, &incoming_cell_ids_status);
                } else {
                    MPI_Wait(&incoming_cell_ids_request, &incoming_cell_ids_status);
                }
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

    @nogc
    ref SolidFVCell get_mapped_cell(size_t i)
    {
        if (i < mapped_cells.length) {
            return mapped_cells[i];
        } else {
            string msg = "Reference to requested mapped-cell is not available.";
            debug { msg ~= format(" cell[%d]", i); }
            throw new FlowSolverException(msg);
        }
    }
         
    // not @nogc because we may set length and use MPI_Irecv
    void exchange_solid_data_phase0()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localSolidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the geometry data via messages.
                //
                // Prepare to exchange solid data for the boundary cells.
                // To match .copy_values_from(mapped_cells[i]) as defined in solidfvcell.d.
                //
                size_t ne = ghost_cells.length * (this_blk.myConfig.n_flow_time_levels * 2 + 24);
                version(complex_numbers) {
                    ne += ghost_cells.length * (this_blk.myConfig.n_flow_time_levels * 2 + 11);
                }
                if (incoming_solidstate_buf.length < ne) { incoming_solidstate_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the other_blk MPI process.
                incoming_solidstate_tag = make_mpi_tag(other_blk.id, other_face, 5);
                MPI_Irecv(incoming_solidstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_solidstate_tag, MPI_COMM_WORLD, &incoming_solidstate_request);
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

    // not @nogc
    void exchange_solid_data_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localSolidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // Blocking send of this block's geometry data
                // to the corresponding non-blocking receive that was posted
                // in the other MPI process.
                outgoing_solidstate_tag = make_mpi_tag(blk.id, which_boundary, 5);
                size_t ne = ghost_cells.length * (this_blk.myConfig.n_flow_time_levels * 2 + 24);
                version(complex_numbers) {
                    ne += ghost_cells.length * (this_blk.myConfig.n_flow_time_levels * 2 + 11);
                }
                if (outgoing_solidstate_buf.length < ne) { outgoing_solidstate_buf.length = ne; }
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    foreach (j; 0 .. this_blk.myConfig.n_flow_time_levels) {
                        outgoing_solidstate_buf[ii++] = c.e[j].re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.e[j].im; }
                        outgoing_solidstate_buf[ii++] = c.dedt[j].re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.dedt[j].im; }
                    }
                    outgoing_solidstate_buf[ii++] = c.volume.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.volume.im; }
                    outgoing_solidstate_buf[ii++] = c.areaxy.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.areaxy.im; }
                    outgoing_solidstate_buf[ii++] = c.pos.x.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.pos.x.im; }
                    outgoing_solidstate_buf[ii++] = c.pos.y.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.pos.y.im; }
                    outgoing_solidstate_buf[ii++] = c.pos.z.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.pos.z.im; }
                    outgoing_solidstate_buf[ii++] = c.sp.rho.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k.re;
                    outgoing_solidstate_buf[ii++] = c.sp.Cp.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k11.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k12.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k13.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k21.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k22.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k23.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k31.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k32.re;
                    outgoing_solidstate_buf[ii++] = c.sp.k33.re;
                    outgoing_solidstate_buf[ii++] = c.T.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.T.im; }
                    outgoing_solidstate_buf[ii++] = c.de_prev.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.de_prev.im; }
                    outgoing_solidstate_buf[ii++] = c.Q.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.Q.im; }
                    outgoing_solidstate_buf[ii++] = c.dTdx.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.dTdx.im; }
                    outgoing_solidstate_buf[ii++] = c.dTdy.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.dTdy.im; }
                    outgoing_solidstate_buf[ii++] = c.dTdz.re; version(complex_numbers) { outgoing_solidstate_buf[ii++] = c.dTdz.im; }
                    if (c.is_ghost) { outgoing_solidstate_buf[ii++] = 1.0; }
                    else { outgoing_solidstate_buf[ii++] = -1.0; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_solidstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_solidstate_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_solidstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_solidstate_tag, MPI_COMM_WORLD);
                }
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

    // not @nogc
    void exchange_solid_data_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localSolidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_solidstate_request, &incoming_solidstate_status);
                } else {
                    MPI_Wait(&incoming_solidstate_request, &incoming_solidstate_status);
                }
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    foreach (j; 0 .. this_blk.myConfig.n_flow_time_levels) {
                        c.e[j].re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.e[j].im = incoming_solidstate_buf[ii++]; }
                        c.dedt[j].re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.dedt[j].im = incoming_solidstate_buf[ii++]; }
                    }
                    c.volume.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.volume.im = incoming_solidstate_buf[ii++]; }
                    c.areaxy.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.areaxy.im = incoming_solidstate_buf[ii++]; }
                    c.pos.x.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.pos.x.im = incoming_solidstate_buf[ii++]; }
                    c.pos.y.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.pos.y.im = incoming_solidstate_buf[ii++]; }
                    c.pos.z.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.pos.z.im = incoming_solidstate_buf[ii++]; }
                    c.sp.rho.re = incoming_solidstate_buf[ii++];
                    c.sp.k.re = incoming_solidstate_buf[ii++];
                    c.sp.Cp.re = incoming_solidstate_buf[ii++];
                    c.sp.k11.re = incoming_solidstate_buf[ii++];
                    c.sp.k12.re = incoming_solidstate_buf[ii++];
                    c.sp.k13.re = incoming_solidstate_buf[ii++];
                    c.sp.k21.re = incoming_solidstate_buf[ii++];
                    c.sp.k22.re = incoming_solidstate_buf[ii++];
                    c.sp.k23.re = incoming_solidstate_buf[ii++];
                    c.sp.k31.re = incoming_solidstate_buf[ii++];
                    c.sp.k32.re = incoming_solidstate_buf[ii++];
                    c.sp.k33.re = incoming_solidstate_buf[ii++];
                    c.T.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.T.im = incoming_solidstate_buf[ii++]; }
                    c.de_prev.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.de_prev.im = incoming_solidstate_buf[ii++]; }
                    c.Q.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.Q.im = incoming_solidstate_buf[ii++]; }
                    c.dTdx.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.dTdx.im = incoming_solidstate_buf[ii++]; }
                    c.dTdy.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.dTdy.im = incoming_solidstate_buf[ii++]; }
                    c.dTdz.re = incoming_solidstate_buf[ii++]; version(complex_numbers) { c.dTdz.im = incoming_solidstate_buf[ii++]; }
                    if ( incoming_solidstate_buf[ii++] > 0.0)  { c.is_ghost = true; }
                    else { c.is_ghost = false; }
                }
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].copy_values_from(mapped_cells[i]);
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].copy_values_from(mapped_cells[i]);
            }
        }
    } // end exchange_geometry_phase2()

    @nogc
    override void apply_for_interface(double t, int tLevel, SolidFVInterface f)
    {
	throw new Error("SolidGCE_SolidGhostCellFullFaceCopy.apply_for_interface() not implemented");
    }

    @nogc
    override void apply(double t, int ftl)
    {
        // We presume that all of the exchange of data happened earlier,
        // and that the ghost cells have been filled with flow state data
        // from their respective source cells.
    } // end apply_structured_grid()
    
} // end class GhostCellFullFaceCopy
