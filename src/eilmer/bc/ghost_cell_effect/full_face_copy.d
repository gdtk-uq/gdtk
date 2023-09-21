// full_face_copy.d

module bc.ghost_cell_effect.full_face_copy;

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
        int outgoing_viscous_gradient_tag, incoming_viscous_gradient_tag;
        MPI_Request incoming_viscous_gradient_request;
        MPI_Status incoming_viscous_gradient_status;
        double[] outgoing_viscous_gradient_buf, incoming_viscous_gradient_buf;
    }

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient,
         bool reorient_vector_quantities,
         ref const(double[]) Rmatrix)
    {
        super(id, boundary, "FullFaceCopy");
        neighbourBlock = cast(FluidBlock) globalBlocks[otherBlock];
        assert(neighbourBlock !is null, "Oops, this should be a FluidBlock object.");
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

    int check_cell_mapping()
    {
        // Returns 0 if the numbers of cells along exchange block boundaries are consistent, -1 otherwise.
        // The checking process stops at the first mismatch.
        //
        // We call this function only after all blocks have been constructed.
        //
        this_blk = cast(SFluidBlock) blk;
        if (!this_blk) { writeln("Destination FlowBlock must be a structured-grid block."); return -4; }
        bool nghost3 = (this_blk.n_ghost_cell_layers == 3);
        other_blk = cast(SFluidBlock) neighbourBlock;
        if (!other_blk) { writeln("Source FlowBlock must be a structured-grid block.");  return -4; }
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
            // First, check consistency of numbers of cells along the exchange boundaries.
            string fstr = "Block[%d] %s exchange with Block[%d] %s, inconsistent number of cells.";
            switch (which_boundary) {
            case Face.north:
                switch (other_face) {
                case Face.north:
                    if (this_blk.nic != other_blk.nic) {
                        writefln(fstr, this_blk.id, "north", other_blk.id, "north");
                        return -1;
                    }
                    break;
                case Face.east:
                    if (this_blk.nic != other_blk.njc) {
                        writefln(fstr, this_blk.id, "north", other_blk.id, "east");
                        return -1;
                    }
                    break;
                case Face.south:
                    if (this_blk.nic != other_blk.nic) {
                        writefln(fstr, this_blk.id, "north", other_blk.id, "south");
                        return -1;
                    }
                    break;
                case Face.west:
                    if (this_blk.nic != other_blk.njc) {
                        writefln(fstr, this_blk.id, "north", other_blk.id, "west");
                        return -1;
                    }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            case Face.east:
                switch (other_face) {
                case Face.north:
                    if (this_blk.njc != other_blk.nic) {
                        writefln(fstr, this_blk.id, "east", other_blk.id, "north");
                        return -1;
                    }
                    break;
                case Face.east:
                    if (this_blk.njc != other_blk.njc) {
                        writefln(fstr, this_blk.id, "east", other_blk.id, "east");
                        return -1;
                    }
                    break;
                case Face.south:
                    if (this_blk.njc != other_blk.nic) {
                        writefln(fstr, this_blk.id, "east", other_blk.id, "south");
                        return -1;
                    }
                    break;
                case Face.west:
                    if (this_blk.njc != other_blk.njc) {
                        writefln(fstr, this_blk.id, "east", other_blk.id, "west");
                        return -1;
                    }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            case Face.south:
                switch (other_face) {
                case Face.north:
                    if (this_blk.nic != other_blk.nic) {
                        writefln(fstr, this_blk.id, "south", other_blk.id, "north");
                        return -1;
                    }
                    break;
                case Face.east:
                    if (this_blk.nic != other_blk.njc) {
                        writefln(fstr, this_blk.id, "south", other_blk.id, "east");
                        return -1;
                    }
                    break;
                case Face.south:
                    if (this_blk.nic != other_blk.nic) {
                        writefln(fstr, this_blk.id, "south", other_blk.id, "south");
                        return -1;
                    }
                    break;
                case Face.west:
                    if (this_blk.nic != other_blk.njc) {
                        writefln(fstr, this_blk.id, "south", other_blk.id, "west");
                        return -1;
                    }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            case Face.west:
                switch (other_face) {
                case Face.north:
                    if (this_blk.njc != other_blk.nic) {
                        writefln(fstr, this_blk.id, "west", other_blk.id, "north");
                        return -1;
                    }
                    break;
                case Face.east:
                    if (this_blk.njc != other_blk.njc) {
                        writefln(fstr, this_blk.id, "west", other_blk.id, "east");
                        return -1;
                    }
                    break;
                case Face.south:
                    if (this_blk.njc != other_blk.nic) {
                        writefln(fstr, this_blk.id, "west", other_blk.id, "south");
                        return -1;
                    }
                    break;
                case Face.west:
                    if (this_blk.njc != other_blk.njc) {
                        writefln(fstr, this_blk.id, "west", other_blk.id, "west");
                        return -1;
                    }
                    break;
                default:
                    writefln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            default:
                writeln("Incorrect boundary connection, invalid which_boundary id.");
                return -3;
            } // end switch which_boundary
        } else {
            // presume dimensions == 3
            // Continue on with 3D work...
            // For checking orientations, see the code in bc_conn.lua.
            string fstr = "Block[%d] %s exchange with Block[%d] %s, orientation=%d, inconsistent number of cells.";
            bool ok = false;
            switch (which_boundary) {
            case Face.north:
                switch (other_face) {
                case Face.north:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "north", other_blk.id, "north", other_orientation);  return -1; }
                    break;
                case Face.east:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "north", other_blk.id, "east", other_orientation); return -1; }
                    break;
                case Face.south:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "north", other_blk.id, "south", other_orientation); return -1;  }
                    break;
                case Face.west:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "north", other_blk.id, "west", other_orientation); return -1; }
                    break;
                case Face.top:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "north", other_blk.id, "top", other_orientation); return -1; }
                    break;
                case Face.bottom:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "north", other_blk.id, "bottom", other_orientation); return -1; }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            case Face.east:
                switch (other_face) {
                case Face.north:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "east", other_blk.id, "north", other_orientation); return -1; }
                    break;
                case Face.east:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "east", other_blk.id, "east", other_orientation); return -1; }
                    break;
                case Face.south:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "east", other_blk.id, "south", other_orientation); return -1; }
                    break;
                case Face.west:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "east", other_blk.id, "west", other_orientation); return -1; }
                    break;
                case Face.top:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "east", other_blk.id, "top", other_orientation); return -1; }
                    break;
                case Face.bottom:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "east", other_blk.id, "bottom", other_orientation); return -1; }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            case Face.south:
                switch (other_face) {
                case Face.north:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "south", other_blk.id, "north", other_orientation); return -1; }
                    break;
                case Face.east:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "south", other_blk.id, "east", other_orientation); return -1; }
                    break;
                case Face.south:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "south", other_blk.id, "south", other_orientation); return -1; }
                    break;
                case Face.west:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "south", other_blk.id, "west", other_orientation); return -1; }
                    break;
                case Face.top:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "south", other_blk.id, "top", other_orientation); return -1; }
                    break;
                case Face.bottom:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "south", other_blk.id, "bottom", other_orientation); return -1; }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            case Face.west:
                switch (other_face) {
                case Face.north:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "west", other_blk.id, "north", other_orientation); return -1; }
                    break;
                case Face.east:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "west", other_blk.id, "east", other_orientation); return -1; }
                    break;
                case Face.south:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "west", other_blk.id, "south", other_orientation); return -1; }
                    break;
                case Face.west:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 1: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc); break;
                    case 2: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nkc); break;
                    case 3: ok = (this_blk.njc == other_blk.nkc) && (this_blk.nkc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "west", other_blk.id, "west", other_orientation); return -1; }
                    break;
                case Face.top:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "west", other_blk.id, "top", other_orientation); return -1; }
                    break;
                case Face.bottom:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 1: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic); break;
                    case 2: ok = (this_blk.njc == other_blk.nic) && (this_blk.nkc == other_blk.njc); break;
                    case 3: ok = (this_blk.njc == other_blk.njc) && (this_blk.nkc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "west", other_blk.id, "bottom", other_orientation); return -1; }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -1;
                } // end switch other_face
                break;
            case Face.top:
                switch (other_face) {
                case Face.north:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "top", other_blk.id, "north", other_orientation); return -1; }
                    break;
                case Face.east:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "top", other_blk.id, "east", other_orientation); return -1; }
                    break;
                case Face.south:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "top", other_blk.id, "south", other_orientation); return -1; }
                    break;
                case Face.west:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "top", other_blk.id, "west", other_orientation); return -1; }
                    break;
                case Face.top:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "top", other_blk.id, "top", other_orientation); return -1; }
                    break;
                case Face.bottom:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "top", other_blk.id, "bottom", other_orientation); return -1; }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -1;
                } // end switch other_face
                break;
            case Face.bottom:
                switch (other_face) {
                case Face.north:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "bottom", other_blk.id, "north", other_orientation); return -1; }
                    break;
                case Face.east:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "bottom", other_blk.id, "east", other_orientation); return -1; }
                    break;
                case Face.south:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "bottom", other_blk.id, "south", other_orientation); return -1; }
                    break;
                case Face.west:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 1: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc); break;
                    case 2: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nkc); break;
                    case 3: ok = (this_blk.nic == other_blk.nkc) && (this_blk.njc == other_blk.njc);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "bottom", other_blk.id, "west", other_orientation); return -1; }
                    break;
                case Face.top:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "bottom", other_blk.id, "top", other_orientation); return -1; }
                    break;
                case Face.bottom:
                    final switch (other_orientation) {
                    case 0: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 1: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic); break;
                    case 2: ok = (this_blk.nic == other_blk.nic) && (this_blk.njc == other_blk.njc); break;
                    case 3: ok = (this_blk.nic == other_blk.njc) && (this_blk.njc == other_blk.nic);
                    }
                    if (!ok) { writefln(fstr, this_blk.id, "bottom", other_blk.id, "bottom", other_orientation); return -1; }
                    break;
                default:
                    writeln("Incorrect boundary connection, invalid source face id.");
                    return -2;
                } // end switch other_face
                break;
            default:
                writeln("Incorrect boundary connection, invalid which_boundary id.");
                return -3;
            } // end switch which_boundary
        } // end if dimensions == ...
        //
        return 0; // to arrive here, the numbers of cells along the exchange boundaries are consistent.
    } // end check_cell_mapping()

    void set_up_cell_mapping_phase0()
    {
        // We call this function only after all blocks have been constructed.
        //
        // The following references and names will be handy for the data exchange
        // that occurs for this ghost-cell effect.
        this_blk = cast(SFluidBlock) blk;
        if (!this_blk) { throw new Error("Destination FlowBlock must be a structured-grid block."); }
        bool nghost3 = (this_blk.n_ghost_cell_layers == 3);
        other_blk = cast(SFluidBlock) neighbourBlock;
        if (!other_blk) { throw new Error("Source FlowBlock must be a structured-grid block."); }
        other_face = neighbourFace;
        other_orientation = neighbourOrientation;
        version(mpi_parallel) {
            other_blk_rank = GlobalConfig.mpi_rank_for_block[other_blk.id];
        }
        //
        // Flow data is copied from the other (source) block into this (destination) block.
        // For the source cells, we use indices into the other block of active cells.
        size_t i_src, j_src, k_src;
        // To find ghost-cells in this destination block, we use the boundary face that holds
        // arrays of references that include the attached ghost cells.
        size_t i_dest, j_dest, k_dest;
        //
        if (blk.myConfig.dimensions == 2) {
            // Handle the 2D case separately.
            switch (which_boundary) {
            case Face.north:
                j_dest = this_blk.njc;  // index of the north-most plane of faces
                foreach (i; 0 .. this_blk.nic) {
                    i_dest = i;
                    auto f = this_blk.get_ifj(i_dest,j_dest);
                    foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njc - 1;
                        i_src = other_blk.nic - i - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n);
                        }
                        break;
                    case Face.east:
                        i_src = other_blk.nic - 1;
                        j_src = i;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src);
                        }
                        break;
                    case Face.south:
                        j_src = 0;
                        i_src = i;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n);
                        }
                        break;
                    case Face.west:
                        i_src = 0;
                        j_src = other_blk.njc - i - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src);
                        }
                        break;
                    default:
                        throw new FlowSolverException("Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // i loop
                break;
            case Face.east:
                i_dest = this_blk.nic;  // index of the east-most plane of faces
                foreach (j; 0 .. this_blk.njc) {
                    j_dest = j;
                    auto f = this_blk.get_ifi(i_dest,j_dest);
                    foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njc - 1;
                        i_src = j;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n);
                        }
                        break;
                    case Face.east:
                        i_src = other_blk.nic - 1;
                        j_src = other_blk.njc - j - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src);
                        }
                        break;
                    case Face.south:
                        j_src = 0;
                        i_src = other_blk.nic - j - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n);
                        }
                        break;
                    case Face.west:
                        i_src = 0;
                        j_src = j;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src);
                        }
                        break;
                    default:
                        throw new FlowSolverException("Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // j loop
                break;
            case Face.south:
                j_dest = 0;  // index of the south-most plane of faces
                foreach (i; 0 .. this_blk.nic) {
                    i_dest = i;
                    auto f = this_blk.get_ifj(i_dest,j_dest);
                    foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njc - 1;
                        i_src = i;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n);
                        }
                        break;
                    case Face.east:
                        i_src = other_blk.nic - 1;
                        j_src = other_blk.njc - i - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src);
                        }
                        break;
                    case Face.south:
                        j_src = 0;
                        i_src = other_blk.nic - i - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n);
                        }
                        break;
                    case Face.west:
                        i_src = 0;
                        j_src = i;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src);
                        }
                        break;
                    default:
                        throw new FlowSolverException("Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // i loop
                break;
            case Face.west:
                i_dest = 0;  // index of the west-most plane of faces
                foreach (j; 0 .. this_blk.njc) {
                    j_dest = j;
                    auto f = this_blk.get_ifi(i_dest,j_dest);
                    foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                    switch (other_face) {
                    case Face.north:
                        j_src = other_blk.njc - 1;
                        i_src = other_blk.nic - j - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n);
                        }
                        break;
                    case Face.east:
                        i_src = other_blk.nic - 1;
                        j_src = j;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src);
                        }
                        break;
                    case Face.south:
                        j_src = 0;
                        i_src = j;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n);
                        }
                        break;
                    case Face.west:
                        i_src = 0;
                        j_src = other_blk.njc - j - 1;
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                            mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src);
                        }
                        break;
                    default:
                        throw new FlowSolverException("Incorrect boundary connection, source face.");
                    } // end switch other_face
                } // j loop
                break;
            default:
                throw new FlowSolverException("Incorrect boundary connection, which_boundary.");
            } // end switch which_boundary
        } else {
            // presume dimensions == 3
            // Continue on with 3D work...
            final switch (which_boundary) {
            case Face.north:
                j_dest = this_blk.njc;  // index of the north-most plane of faces
                foreach (i; 0 .. this_blk.nic) {
                    i_dest = i;
                    foreach (k; 0 .. this_blk.nkc) {
                        k_dest = k;
                        auto f = this_blk.get_ifj(i_dest,j_dest,k_dest);
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; k_src = k; break;
                            case 1: i_src = k; k_src = i; break;
                            case 2: i_src = i; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n,k_src);
                            }
                            break;
                        case Face.east:
                            i_src = other_blk.nic - 1;
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = other_blk.njc - k - 1; k_src = i; break;
                            case 2: j_src = other_blk.njc - i - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = k; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src,k_src);
                            }
                            break;
                        case Face.south:
                            j_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; k_src = i; break;
                            case 2: i_src = other_blk.nic - i - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = k; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n,k_src);
                            }
                            break;
                        case Face.west:
                            i_src = 0;
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njc - i - 1; k_src = k; break;
                            case 1: j_src = k; k_src = i; break;
                            case 2: j_src = i; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = other_blk.njc - k - 1; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src,k_src);
                            }
                            break;
                        case Face.top:
                            k_src = other_blk.nkc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; j_src = i; break;
                            case 2: i_src = other_blk.nic - i - 1; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = k; j_src = other_blk.njc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src-n);
                            }
                            break;
                        case Face.bottom:
                            k_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; j_src = k; break;
                            case 1: i_src = k; j_src = i; break;
                            case 2: i_src = i; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; j_src = other_blk.njc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src+n);
                            }
                        } // end switch (other_face)
                    } // k loop
                } // i loop
                break;
            case Face.east:
                i_dest = this_blk.nic;  // index of the east-most plane of faces
                foreach (j; 0 .. this_blk.njc) {
                    j_dest = j;
                    foreach (k; 0 .. this_blk.nkc) {
                        k_dest = k;
                        auto f = this_blk.get_ifi(i_dest,j_dest,k_dest);
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = k; k_src = other_blk.nkc - j - 1; break;
                            case 2: i_src = other_blk.nic - j - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; k_src = j;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n,k_src);
                            }
                            break;
                        case Face.east:
                            i_src = other_blk.nic - 1;
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njc - j - 1; k_src = k; break;
                            case 1: j_src = other_blk.njc - k - 1; k_src = other_blk.nkc - j - 1; break;
                            case 2: j_src = j; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = k; k_src = j;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src,k_src);
                            }
                            break;
                        case Face.south:
                            j_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - j - 1; k_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; k_src = other_blk.nkc - j - 1; break;
                            case 2: i_src = j; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = k; k_src = j;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n,k_src);
                            }
                            break;
                        case Face.west:
                            i_src = 0;
                            final switch (other_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = k; k_src = other_blk.nkc - j - 1; break;
                            case 2: j_src = other_blk.njc - j - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = other_blk.njc - k - 1; k_src = j;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src,k_src);
                            }
                            break;
                        case Face.top:
                            k_src = other_blk.nkc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - j - 1; j_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; j_src = other_blk.njc - j - 1; break;
                            case 2: i_src = j; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = k; j_src = j;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src-n);
                            }
                            break;
                        case Face.bottom:
                            k_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = k; j_src = other_blk.njc - j - 1; break;
                            case 2: i_src = other_blk.nic - j - 1; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; j_src = j;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src+n);
                            }
                        } // end switch (other_face)
                    } // k loop
                } // j loop
                break;
            case Face.south:
                j_dest = 0;  // index of the south-most plane of faces
                foreach (i; 0 .. this_blk.nic) {
                    i_dest = i;
                    foreach (k; 0 .. this_blk.nkc) {
                        k_dest = k;
                        auto f = this_blk.get_ifj(i_dest,j_dest,k_dest);
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = k; k_src = other_blk.nkc - i - 1; break;
                            case 2: i_src = other_blk.nic - i - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n,k_src);
                            }
                            break;
                        case Face.east:
                            i_src = other_blk.nic - 1;
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njc - i - 1; k_src = k; break;
                            case 1: j_src = other_blk.njc - k - 1; k_src = other_blk.nkc - i - 1; break;
                            case 2: j_src = i; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = k; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src,k_src);
                            }
                            break;
                        case Face.south:
                            j_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; k_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; k_src = other_blk.nkc - i - 1; break;
                            case 2: i_src = i; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = k; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n,k_src);
                            }
                            break;
                        case Face.west:
                            i_src = 0;
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = k; k_src = other_blk.nkc - i - 1; break;
                            case 2: j_src = other_blk.njc - i - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = other_blk.njc - k - 1; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src,k_src);
                            }
                            break;
                        case Face.top:
                            k_src = other_blk.nkc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; j_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; j_src = other_blk.njc - i - 1; break;
                            case 2: i_src = i; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = k; j_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src-n);
                            }
                            break;
                        case Face.bottom:
                            k_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = k; j_src = other_blk.njc - i - 1; break;
                            case 2: i_src = other_blk.nic - i - 1; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; j_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src+n);
                            }
                        } // end switch (other_face)
                    } // k loop
                } // i loop
                break;
            case Face.west:
                i_dest = 0;  // index of the west-most plane of faces
                foreach (j; 0 .. this_blk.njc) {
                    j_dest = j;
                    foreach (k; 0 .. this_blk.nkc) {
                        k_dest = k;
                        auto f = this_blk.get_ifi(i_dest,j_dest,k_dest);
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - j - 1; k_src = k; break;
                            case 1: i_src = k; k_src = j; break;
                            case 2: i_src = j; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; k_src = other_blk.nkc - j - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n,k_src);
                            }
                            break;
                        case Face.east:
                            i_src = other_blk.nic - 1;
                            final switch (other_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = other_blk.njc - k - 1; k_src = j; break;
                            case 2: j_src = other_blk.njc - j - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = k; k_src = other_blk.nkc - j - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src,k_src);
                            }
                            break;
                        case Face.south:
                            j_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; k_src = j; break;
                            case 2: i_src = other_blk.nic - j - 1; k_src = other_blk.nkc - k - 1; break;
                            case 3: i_src = k; k_src = other_blk.nkc - j - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n,k_src);
                            }
                            break;
                        case Face.west:
                            i_src = 0;
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njc - j - 1; k_src = k; break;
                            case 1: j_src = k; k_src = j; break;
                            case 2: j_src = j; k_src = other_blk.nkc - k - 1; break;
                            case 3: j_src = other_blk.njc - k - 1; k_src = other_blk.nkc - j - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src,k_src);
                            }
                            break;
                        case Face.top:
                            k_src = other_blk.nkc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = other_blk.nic - k - 1; j_src = j; break;
                            case 2: i_src = other_blk.nic - j - 1; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = k; j_src = other_blk.njc - j - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src-n);
                            }
                            break;
                        case Face.bottom:
                            k_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - j - 1; j_src = k; break;
                            case 1: i_src = k; j_src = j; break;
                            case 2: i_src = j; j_src = other_blk.njc - k - 1; break;
                            case 3: i_src = other_blk.nic - k - 1; j_src = other_blk.njc - j - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src+n);
                            }
                        } // end switch (other_face)
                    } // k loop
                } // j loop
                break;
            case Face.top:
                k_dest = this_blk.nkc;  // index of the top-most plane of faces
                foreach (j; 0 .. this_blk.njc) {
                    j_dest = j;
                    foreach (i; 0 .. this_blk.nic) {
                        i_dest = i;
                        auto f = this_blk.get_ifk(i_dest,j_dest,k_dest);
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.right_cells[n]; }
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = j; k_src = other_blk.nkc - i - 1; break;
                            case 2: i_src = other_blk.nic - i - 1; k_src = other_blk.nkc - j - 1; break;
                            case 3: i_src = other_blk.nic - j - 1; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n,k_src);
                            }
                            break;
                        case Face.east:
                            i_src = other_blk.nic - 1;
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njc - i - 1; k_src = j; break;
                            case 1: j_src = other_blk.njc - j - 1; k_src = other_blk.nkc - i - 1; break;
                            case 2: j_src = i; k_src = other_blk.nkc - j - 1; break;
                            case 3: j_src = j; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src,k_src);
                            }
                            break;
                        case Face.south:
                            j_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; k_src = j; break;
                            case 1: i_src = other_blk.nic - j - 1; k_src = other_blk.nkc - i - 1; break;
                            case 2: i_src = i; k_src = other_blk.nkc - j - 1; break;
                            case 3: i_src = j; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n,k_src);
                            }
                            break;
                        case Face.west:
                            i_src = 0;
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = j; k_src = other_blk.nkc - i - 1; break;
                            case 2: j_src = other_blk.njc - i - 1; k_src = other_blk.nkc - j - 1; break;
                            case 3: j_src = other_blk.njc - j - 1; k_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src,k_src);
                            }
                            break;
                        case Face.top:
                            k_src = other_blk.nkc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; j_src = j; break;
                            case 1: i_src = other_blk.nic - j - 1; j_src = other_blk.njc - i - 1; break;
                            case 2: i_src = i; j_src = other_blk.njc - j - 1; break;
                            case 3: i_src = j; j_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src-n);
                            }
                            break;
                        case Face.bottom:
                            k_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = j; j_src = other_blk.njc - i - 1; break;
                            case 2: i_src = other_blk.nic - i - 1; j_src = other_blk.njc - j - 1; break;
                            case 3: i_src = other_blk.nic - j - 1; j_src = i;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src+n);
                            }
                        } // end switch (other_face)
                    } // i loop
                } // j loop
                break;
            case Face.bottom:
                k_dest = 0;  // index of the bottom-most plane of faces
                foreach (j; 0 .. this_blk.njc) {
                    j_dest = j;
                    foreach (i; 0 .. this_blk.nic) {
                        i_dest = i;
                        auto f = this_blk.get_ifk(i_dest,j_dest,k_dest);
                        foreach (n; 0 .. this_blk.n_ghost_cell_layers) { ghost_cells ~= f.left_cells[n]; }
                        final switch (other_face) {
                        case Face.north:
                            j_src = other_blk.njc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; k_src = j; break;
                            case 1: i_src = j; k_src = i; break;
                            case 2: i_src = i; k_src = other_blk.nkc - j - 1; break;
                            case 3: i_src = other_blk.nic - j - 1; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src-n,k_src);
                            }
                            break;
                        case Face.east:
                            i_src = other_blk.nic - 1;
                            final switch (other_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = other_blk.njc - j - 1; k_src = i; break;
                            case 2: j_src = other_blk.njc - i - 1; k_src = other_blk.nkc - j - 1; break;
                            case 3: j_src = j; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src-n,j_src,k_src);
                            }
                            break;
                        case Face.south:
                            j_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = other_blk.nic - j - 1; k_src = i; break;
                            case 2: i_src = other_blk.nic - i - 1; k_src = other_blk.nkc - j - 1; break;
                            case 3: i_src = j; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src+n,k_src);
                            }
                            break;
                        case Face.west:
                            i_src = 0;
                            final switch (other_orientation) {
                            case 0: j_src = other_blk.njc - i - 1; k_src = j; break;
                            case 1: j_src = j; k_src = i; break;
                            case 2: j_src = i; k_src = other_blk.nkc - j - 1; break;
                            case 3: j_src = other_blk.njc - j - 1; k_src = other_blk.nkc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src+n,j_src,k_src);
                            }
                            break;
                        case Face.top:
                            k_src = other_blk.nkc - 1;
                            final switch (other_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = other_blk.nic - j - 1; j_src = i; break;
                            case 2: i_src = other_blk.nic - i - 1; j_src = other_blk.njc - j - 1; break;
                            case 3: i_src = j; j_src = other_blk.njc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src-n);
                            }
                            break;
                        case Face.bottom:
                            k_src = 0;
                            final switch (other_orientation) {
                            case 0: i_src = other_blk.nic - i - 1; j_src = j; break;
                            case 1: i_src = j; j_src = i; break;
                            case 2: i_src = i; j_src = other_blk.njc - j - 1; break;
                            case 3: i_src = other_blk.nic - j - 1; j_src = other_blk.njc - i - 1;
                            }
                            foreach (n; 0 .. this_blk.n_ghost_cell_layers) {
                                mapped_cell_ids ~= other_blk.cell_index(i_src,j_src,k_src+n);
                            }
                        } // end switch other_face
                    } // i loop
                } // j loop
            } // end switch which_boundary
        } // end if dimensions == ...
        //
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
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

    // not @nogc because we set array length
    void set_up_cell_mapping_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                // Blocking send to corresponding non-blocking receive that was posted
                // at in other_blk MPI process.
                size_t ne = ghost_cells.length;
                if (outgoing_cell_ids_buf.length < ne) { outgoing_cell_ids_buf.length = ne; }
                assert(ne == mapped_cell_ids.length, "oops, wrong length");
                outgoing_cell_ids_tag = make_mpi_tag(blk.id, which_boundary, 0);
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
    } // end set_up_cell_mapping_phase1()

    // not @nogc
    void set_up_cell_mapping_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
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
        // Ghost cells need to be flagged as interior, for later (NNG 23/08/22)
        foreach (c; ghost_cells) c.is_interior_to_domain = true;
    } // end set_up_cell_mapping_phase2()

    @nogc
    ref FVCell get_mapped_cell(size_t i)
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
    void exchange_geometry_phase0()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the geometry data via messages.
                //
                // Prepare to exchange geometry data for the boundary cells.
                // To match .copy_values_from(mapped_cells[i], CopyDataOption.grid) as defined in fvcell.d.
                //
                size_t ne = ghost_cells.length * (this_blk.myConfig.n_grid_time_levels * 5 + 5);
                version(complex_numbers) { ne *= 2; }
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

    // not @nogc
    void exchange_geometry_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // Blocking send of this block's geometry data
                // to the corresponding non-blocking receive that was posted
                // in the other MPI process.
                outgoing_geometry_tag = make_mpi_tag(blk.id, which_boundary, 1);
                size_t ne = ghost_cells.length * (this_blk.myConfig.n_grid_time_levels * 5 + 5);
                version(complex_numbers) { ne *= 2; }
                if (outgoing_geometry_buf.length < ne) { outgoing_geometry_buf.length = ne; }
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    foreach (j; 0 .. this_blk.myConfig.n_grid_time_levels) {
                        outgoing_geometry_buf[ii++] = c.pos[j].x.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.pos[j].x.im; }
                        outgoing_geometry_buf[ii++] = c.pos[j].y.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.pos[j].y.im; }
                        outgoing_geometry_buf[ii++] = c.pos[j].z.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.pos[j].z.im; }
                        outgoing_geometry_buf[ii++] = c.volume[j].re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.volume[j].im; }
                        outgoing_geometry_buf[ii++] = c.areaxy[j].re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.areaxy[j].im; }
                    }
                    outgoing_geometry_buf[ii++] = c.iLength.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.iLength.im; }
                    outgoing_geometry_buf[ii++] = c.jLength.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.jLength.im; }
                    outgoing_geometry_buf[ii++] = c.kLength.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.kLength.im; }
                    outgoing_geometry_buf[ii++] = c.L_min.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.L_min.im; }
                    outgoing_geometry_buf[ii++] = c.L_max.re; version(complex_numbers) { outgoing_geometry_buf[ii++] = c.L_max.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_geometry_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_geometry_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_geometry_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_geometry_tag, MPI_COMM_WORLD);
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
    void exchange_geometry_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The other block is in another MPI process, go fetch the data via messages.
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_geometry_request, &incoming_geometry_status);
                } else {
                    MPI_Wait(&incoming_geometry_request, &incoming_geometry_status);
                }
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    foreach (j; 0 .. this_blk.myConfig.n_grid_time_levels) {
                        c.pos[j].x.re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.pos[j].x.im = incoming_geometry_buf[ii++]; }
                        c.pos[j].y.re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.pos[j].y.im = incoming_geometry_buf[ii++]; }
                        c.pos[j].z.re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.pos[j].z.im = incoming_geometry_buf[ii++]; }
                        c.volume[j].re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.volume[j].im = incoming_geometry_buf[ii++]; }
                        c.areaxy[j].re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.areaxy[j].im = incoming_geometry_buf[ii++]; }
                    }

                    static foreach(li; 0 .. 3) {
                        c.lengths[li].re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.lengths[li].im = incoming_geometry_buf[ii++]; }
                    }

                    c.L_min.re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.L_min.im = incoming_geometry_buf[ii++]; }
                    c.L_max.re = incoming_geometry_buf[ii++]; version(complex_numbers) { c.L_max.im = incoming_geometry_buf[ii++]; }
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
        // Need to have the ghost cells to have the correct position for the viscous calc
        if (reorient_vector_quantities) {
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].pos[0].apply_matrix_transform(Rmatrix);
            }
        }
    } // end exchange_geometry_phase2()

    @nogc
    size_t flowstate_buffer_entry_size(const LocalConfig myConfig){
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
        nitems += nmodes*3;
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
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
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
                //size_t nspecies = this_blk.myConfig.n_species;
                //size_t nmodes = this_blk.myConfig.n_modes;
                //size_t ne = ghost_cells.length * (nmodes*3 + nspecies + 23);
                size_t fs_size = flowstate_buffer_entry_size(this_blk.myConfig);
                size_t ne = ghost_cells.length * fs_size;
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

    // not @nogc
    void exchange_flowstate_phase1(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
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
                size_t nspecies = this_blk.myConfig.n_species;
                size_t nmodes = this_blk.myConfig.n_modes;
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // in the other MPI process.
                //size_t nitems = 16;
                //version(MHD) { nitems += 5; }
                //version(turbulence) { nitems += 2; }
                //size_t ne = ghost_cells.length * (nmodes*3 + nspecies + nitems);

                size_t fs_size = flowstate_buffer_entry_size(this_blk.myConfig);
                size_t ne = ghost_cells.length * fs_size;
                size_t nturb = this_blk.myConfig.turb_model.nturb;
                if (outgoing_flowstate_buf.length < ne) { outgoing_flowstate_buf.length = ne; }
                outgoing_flowstate_tag = make_mpi_tag(blk.id, which_boundary, 0);
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    FlowState* fs = &(c.fs);
                    GasState* gs = &(fs.gas);
                    outgoing_flowstate_buf[ii++] = gs.rho.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.rho.im; }
                    outgoing_flowstate_buf[ii++] = gs.p.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.p.im; }
                    outgoing_flowstate_buf[ii++] = gs.T.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.T.im; }
                    outgoing_flowstate_buf[ii++] = gs.u.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.u.im; }
                    outgoing_flowstate_buf[ii++] = gs.p_e.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.p_e.im; }
                    outgoing_flowstate_buf[ii++] = gs.a.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.a.im; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { outgoing_flowstate_buf[ii++] = gs.u_modes[j].re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.u_modes[j].im; } }
                        foreach (j; 0 .. nmodes) { outgoing_flowstate_buf[ii++] = gs.T_modes[j].re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.T_modes[j].im; } }
                    }
                    outgoing_flowstate_buf[ii++] = gs.mu.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.mu.im; }
                    outgoing_flowstate_buf[ii++] = gs.k.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.k.im; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { outgoing_flowstate_buf[ii++] = gs.k_modes[j].re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.k_modes[j].im; } }
                    }
                    outgoing_flowstate_buf[ii++] = gs.sigma.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.sigma.im; }
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) { outgoing_flowstate_buf[ii++] = gs.massf[j].re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.massf[j].im; } }
                        foreach (j; 0 .. nspecies) { outgoing_flowstate_buf[ii++] = gs.rho_s[j].re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.rho_s[j].im; } }
                    } else {
                        outgoing_flowstate_buf[ii++] = 1.0;  version(complex_numbers) { outgoing_flowstate_buf[ii++] = 0.0; } // single-species mass fraction
                    }
                    outgoing_flowstate_buf[ii++] = gs.quality.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = gs.quality.im; }
                    outgoing_flowstate_buf[ii++] = fs.vel.x.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.vel.x.im; }
                    outgoing_flowstate_buf[ii++] = fs.vel.y.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.vel.y.im; }
                    outgoing_flowstate_buf[ii++] = fs.vel.z.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.vel.z.im; }
                    version(MHD) {
                        outgoing_flowstate_buf[ii++] = fs.B.x.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.B.x.im; }
                        outgoing_flowstate_buf[ii++] = fs.B.y.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.B.y.im; }
                        outgoing_flowstate_buf[ii++] = fs.B.z.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.B.z.im; }
                        outgoing_flowstate_buf[ii++] = fs.psi.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.psi.im; }
                        outgoing_flowstate_buf[ii++] = fs.divB.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.divB.im; }
                    }
                    version(turbulence) {
                        foreach (j; 0 .. nturb) { outgoing_flowstate_buf[ii++] = fs.turb[j].re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.turb[j].im; } }
                    }
                    outgoing_flowstate_buf[ii++] = fs.mu_t.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.mu_t.im; }
                    outgoing_flowstate_buf[ii++] = fs.k_t.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.k_t.im; }
                    outgoing_flowstate_buf[ii++] = fs.S.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = fs.S.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_flowstate_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_flowstate_tag, MPI_COMM_WORLD);
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
        // Done with copying from source cells.
    } // end exchange_flowstate_phase1()

    // not @nogc
    void exchange_flowstate_phase2(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The source block is in another MPI process, go fetch the data via messages.
                //
                size_t nspecies = this_blk.myConfig.n_species;
                size_t nmodes = this_blk.myConfig.n_modes;
                size_t nturb = this_blk.myConfig.turb_model.nturb;
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_flowstate_request, &incoming_flowstate_status);
                } else {
                    MPI_Wait(&incoming_flowstate_request, &incoming_flowstate_status);
                }
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    FlowState* fs = &(c.fs);
                    GasState* gs = &(fs.gas);
                    gs.rho.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.rho.im = incoming_flowstate_buf[ii++]; }
                    gs.p.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.p.im = incoming_flowstate_buf[ii++]; }
                    gs.T.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.T.im = incoming_flowstate_buf[ii++]; }
                    gs.u.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.u.im = incoming_flowstate_buf[ii++]; }
                    gs.p_e.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.p_e.im = incoming_flowstate_buf[ii++]; }
                    gs.a.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.a.im = incoming_flowstate_buf[ii++]; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { gs.u_modes[j].re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.u_modes[j].im = incoming_flowstate_buf[ii++]; } }
                        foreach (j; 0 .. nmodes) { gs.T_modes[j].re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.T_modes[j].im = incoming_flowstate_buf[ii++]; } }
                    }
                    gs.mu.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.mu.im = incoming_flowstate_buf[ii++]; }
                    gs.k.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.k.im = incoming_flowstate_buf[ii++]; }
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) { gs.k_modes[j].re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.k_modes[j].im = incoming_flowstate_buf[ii++]; } }
                    }
                    gs.sigma.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.sigma.im = incoming_flowstate_buf[ii++]; }
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) { gs.massf[j].re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.massf[j].im = incoming_flowstate_buf[ii++]; } }
                        foreach (j; 0 .. nspecies) { gs.rho_s[j].re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.rho_s[j].im = incoming_flowstate_buf[ii++]; } }
                    } else {
                        double junk = incoming_flowstate_buf[ii++]; version(complex_numbers) { double more_junk = incoming_flowstate_buf[ii++]; }
                    }
                    gs.quality.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { gs.quality.im = incoming_flowstate_buf[ii++]; }
                    fs.vel.x.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.vel.x.im = incoming_flowstate_buf[ii++]; }
                    fs.vel.y.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.vel.y.im = incoming_flowstate_buf[ii++]; }
                    fs.vel.z.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.vel.z.im = incoming_flowstate_buf[ii++]; }
                    version(MHD) {
                        fs.B.x.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.B.x.im = incoming_flowstate_buf[ii++]; }
                        fs.B.y.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.B.y.im = incoming_flowstate_buf[ii++]; }
                        fs.B.z.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.B.z.im = incoming_flowstate_buf[ii++]; }
                        fs.psi.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.psi.im = incoming_flowstate_buf[ii++]; }
                        fs.divB.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.divB.im = incoming_flowstate_buf[ii++]; }
                    }
                    version(turbulence) {
                        foreach(j; 0 .. nturb) { fs.turb[j].re = incoming_flowstate_buf[ii++];  version(complex_numbers) { fs.turb[j].im = incoming_flowstate_buf[ii++]; } }
                    }
                    fs.mu_t.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.mu_t.im = incoming_flowstate_buf[ii++]; }
                    fs.k_t.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.k_t.im = incoming_flowstate_buf[ii++]; }
                    fs.S.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { fs.S.im = incoming_flowstate_buf[ii++]; }
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


    // not @nogc
    void exchange_shock_phase0(double t, int gtl, int ftl)
    // look in exchange_flowstate_phase0 for description of activities
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                size_t fs_size = 1; // we only want to transfer a single value
                size_t ne = ghost_cells.length * fs_size;
                version(complex_numbers) { ne *= 2; }
                if (incoming_flowstate_buf.length < ne) { incoming_flowstate_buf.length = ne; }
                incoming_flowstate_tag = make_mpi_tag(other_blk.id, other_face, 0);
                MPI_Irecv(incoming_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_flowstate_tag, MPI_COMM_WORLD, &incoming_flowstate_request);
            } else {
            }
        } else {
        }
    } // end exchange_shock_phase0()

    // not @nogc
    void exchange_shock_phase1(double t, int gtl, int ftl)
    // look in exchange_flowstate_phase1 for description of activities
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                //
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                const size_t fs_size = 1;
                size_t ne = ghost_cells.length * fs_size;
                version(complex_numbers) { ne *= 2; }
                if (outgoing_flowstate_buf.length < ne) { outgoing_flowstate_buf.length = ne; }
                outgoing_flowstate_tag = make_mpi_tag(blk.id, which_boundary, 0);
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    outgoing_flowstate_buf[ii++] = c.fs.S.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = c.fs.S.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_flowstate_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_flowstate_tag, MPI_COMM_WORLD);
                }
            } else {
            }
        } else {
        }
    } // end exchange_shock_phase1()

    // not @nogc
    void exchange_shock_phase2(double t, int gtl, int ftl)
    // look in exchange_flowstate_phase2 for description of activities
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_flowstate_request, &incoming_flowstate_status);
                } else {
                    MPI_Wait(&incoming_flowstate_request, &incoming_flowstate_status);
                }
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    c.fs.S.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { c.fs.S.im = incoming_flowstate_buf[ii++]; }
                }
            } else {
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].fs.S = mapped_cells[i].fs.S;
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].fs.S = mapped_cells[i].fs.S;
            }
        }
        // Done with copying from source cells.
    } // end exchange_flowstate_phase2()

    void exchange_turbulent_transprops_phase0()
    /*
        Exchange only mu_t and k_t between blocks, using MPI if required.

        @author: Nick Gibbons (20/08/21)
    */
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                size_t fs_size = 2;
                size_t ne = ghost_cells.length * fs_size;
                version(complex_numbers) { ne *= 2; }

                // We'll reuse the flowstate tags and buffers. This *should* be okay.
                if (incoming_flowstate_buf.length < ne) { incoming_flowstate_buf.length = ne; }
                incoming_flowstate_tag = make_mpi_tag(other_blk.id, other_face, 0);
                MPI_Irecv(incoming_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_flowstate_tag, MPI_COMM_WORLD, &incoming_flowstate_request);
            } else {
            }
        } else {
        }
    }

    void exchange_turbulent_transprops_phase1()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");

                const size_t fs_size = 2;
                size_t ne = ghost_cells.length * fs_size;
                version(complex_numbers) { ne *= 2; }
                if (outgoing_flowstate_buf.length < ne) { outgoing_flowstate_buf.length = ne; }
                outgoing_flowstate_tag = make_mpi_tag(blk.id, which_boundary, 0);
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    outgoing_flowstate_buf[ii++] = c.fs.mu_t.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = c.fs.mu_t.im; }
                    outgoing_flowstate_buf[ii++] = c.fs.k_t.re; version(complex_numbers) { outgoing_flowstate_buf[ii++] = c.fs.k_t.im; }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_flowstate_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(outgoing_flowstate_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_flowstate_tag, MPI_COMM_WORLD);
                }
            } else {
            }
        } else {
        }
    }

    void exchange_turbulent_transprops_phase2()
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_flowstate_request, &incoming_flowstate_status);
                } else {
                    MPI_Wait(&incoming_flowstate_request, &incoming_flowstate_status);
                }
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    c.fs.mu_t.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { c.fs.mu_t.im = incoming_flowstate_buf[ii++]; }
                    c.fs.k_t.re = incoming_flowstate_buf[ii++]; version(complex_numbers) { c.fs.k_t.im = incoming_flowstate_buf[ii++]; }
                }
            } else {
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].fs.mu_t = mapped_cells[i].fs.mu_t;
                    ghost_cells[i].fs.k_t = mapped_cells[i].fs.k_t;
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].fs.mu_t = mapped_cells[i].fs.mu_t;
                ghost_cells[i].fs.k_t = mapped_cells[i].fs.k_t;
            }
        }
    }


    @nogc
    size_t viscous_gradient_buffer_entry_size(const LocalConfig myConfig){
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
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
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
                //size_t nspecies = this_blk.myConfig.n_species;
                //size_t nmodes = this_blk.myConfig.n_modes;
                //size_t nitems = 12;
                //version(turbulence) { nitems += 6; }
                size_t grad_size = viscous_gradient_buffer_entry_size(this_blk.myConfig);
                size_t ne = ghost_cells.length * grad_size;
                if (incoming_viscous_gradient_buf.length < ne) { incoming_viscous_gradient_buf.length = ne; }
                //
                // Post non-blocking receive for geometry data that we expect to receive later
                // from the other_blk MPI process.
                incoming_viscous_gradient_tag = make_mpi_tag(other_blk.id, other_face, 2);
                MPI_Irecv(incoming_viscous_gradient_buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                          incoming_viscous_gradient_tag, MPI_COMM_WORLD, &incoming_viscous_gradient_request);
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
    } // end exchange_viscous_gradient_phase0()

    // not @nogc
    void exchange_viscous_gradient_phase1(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
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
                size_t nspecies = this_blk.myConfig.n_species;
                size_t nmodes = this_blk.myConfig.n_modes;
                size_t nturb = this_blk.myConfig.turb_model.nturb;
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                // Blocking send of this block's flow data
                // to the corresponding non-blocking receive that was posted
                // in the other MPI process.
                //size_t nitems = 12;
                //version(turbulence) { nitems += 6; }
                size_t grad_size = viscous_gradient_buffer_entry_size(this_blk.myConfig);
                size_t ne = ghost_cells.length * grad_size;
                if (outgoing_viscous_gradient_buf.length < ne) { outgoing_viscous_gradient_buf.length = ne; }
                outgoing_viscous_gradient_tag = make_mpi_tag(blk.id, which_boundary, 2);
                auto buf = outgoing_viscous_gradient_buf;
                size_t ii = 0;
                foreach (c; outgoing_mapped_cells) {
                    // velocity
                    buf[ii++] = c.grad.vel[0][0].re; version(complex_numbers) { buf[ii++] = c.grad.vel[0][0].im; }
                    buf[ii++] = c.grad.vel[0][1].re; version(complex_numbers) { buf[ii++] = c.grad.vel[0][1].im; }
                    buf[ii++] = c.grad.vel[0][2].re; version(complex_numbers) { buf[ii++] = c.grad.vel[0][2].im; }
                    buf[ii++] = c.grad.vel[1][0].re; version(complex_numbers) { buf[ii++] = c.grad.vel[1][0].im; }
                    buf[ii++] = c.grad.vel[1][1].re; version(complex_numbers) { buf[ii++] = c.grad.vel[1][1].im; }
                    buf[ii++] = c.grad.vel[1][2].re; version(complex_numbers) { buf[ii++] = c.grad.vel[1][2].im; }
                    buf[ii++] = c.grad.vel[2][0].re; version(complex_numbers) { buf[ii++] = c.grad.vel[2][0].im; }
                    buf[ii++] = c.grad.vel[2][1].re; version(complex_numbers) { buf[ii++] = c.grad.vel[2][1].im; }
                    buf[ii++] = c.grad.vel[2][2].re; version(complex_numbers) { buf[ii++] = c.grad.vel[2][2].im; }
                    // rho, p, T, u
                    buf[ii++] = c.grad.T[0].re; version(complex_numbers) { buf[ii++] = c.grad.T[0].im; }
                    buf[ii++] = c.grad.T[1].re; version(complex_numbers) { buf[ii++] = c.grad.T[1].im; }
                    buf[ii++] = c.grad.T[2].re; version(complex_numbers) { buf[ii++] = c.grad.T[2].im; }
                    // tke, omega
                    version(turbulence) {
                        foreach (j; 0 .. nturb) {
                            buf[ii++] = c.grad.turb[j][0].re; version(complex_numbers) { buf[ii++] = c.grad.turb[j][0].im; }
                            buf[ii++] = c.grad.turb[j][1].re; version(complex_numbers) { buf[ii++] = c.grad.turb[j][1].im; }
                            buf[ii++] = c.grad.turb[j][2].re; version(complex_numbers) { buf[ii++] = c.grad.turb[j][2].im; }
                        }
                    }
                    // multi-species
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) {
                            buf[ii++] = c.grad.massf[j][0].re; version(complex_numbers) { buf[ii++] = c.grad.massf[j][0].im; }
                            buf[ii++] = c.grad.massf[j][1].re; version(complex_numbers) { buf[ii++] = c.grad.massf[j][1].im; }
                            buf[ii++] = c.grad.massf[j][2].re; version(complex_numbers) { buf[ii++] = c.grad.massf[j][2].im; }
                        }
                    }
                    // multi-T
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) {
                            buf[ii++] = c.grad.T_modes[j][0].re; version(complex_numbers) { buf[ii++] = c.grad.T_modes[j][0].im; }
                            buf[ii++] = c.grad.T_modes[j][1].re; version(complex_numbers) { buf[ii++] = c.grad.T_modes[j][1].im; }
                            buf[ii++] = c.grad.T_modes[j][2].re; version(complex_numbers) { buf[ii++] = c.grad.T_modes[j][2].im; }
                        }
                    }
                }
                version(mpi_timeouts) {
                    MPI_Request send_request;
                    MPI_Isend(buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                              outgoing_viscous_gradient_tag, MPI_COMM_WORLD, &send_request);
                    MPI_Status send_status;
                    MPI_Wait_a_while(&send_request, &send_status);
                } else {
                    MPI_Send(buf.ptr, to!int(ne), MPI_DOUBLE, other_blk_rank,
                             outgoing_viscous_gradient_tag, MPI_COMM_WORLD);
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
        // Done with copying from source cells.
    } // end exchange_viscous_gradient_phase1()

    // not @nogc
    void exchange_viscous_gradient_phase2(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            if (find(GlobalConfig.localFluidBlockIds, other_blk.id).empty) {
                // The source block is in another MPI process, go fetch the data via messages.
                //
                size_t nspecies = this_blk.myConfig.n_species;
                size_t nmodes = this_blk.myConfig.n_modes;
                size_t nturb = this_blk.myConfig.turb_model.nturb;
                assert(outgoing_mapped_cells.length == ghost_cells.length,
                       "oops, mismatch in outgoing_mapped_cells and ghost_cells.");
                //
                // Wait for non-blocking receive to complete.
                // Once complete, copy the data back into the local context.
                version(mpi_timeouts) {
                    MPI_Wait_a_while(&incoming_viscous_gradient_request, &incoming_viscous_gradient_status);
                } else {
                    MPI_Wait(&incoming_viscous_gradient_request, &incoming_viscous_gradient_status);
                }
                auto buf = incoming_viscous_gradient_buf;
                size_t ii = 0;
                foreach (c; ghost_cells) {
                    // velocity
                    c.grad.vel[0][0].re = buf[ii++]; version(complex_numbers) { c.grad.vel[0][0].im = buf[ii++]; }
                    c.grad.vel[0][1].re = buf[ii++]; version(complex_numbers) { c.grad.vel[0][1].im = buf[ii++]; }
                    c.grad.vel[0][2].re = buf[ii++]; version(complex_numbers) { c.grad.vel[0][2].im = buf[ii++]; }
                    c.grad.vel[1][0].re = buf[ii++]; version(complex_numbers) { c.grad.vel[1][0].im = buf[ii++]; }
                    c.grad.vel[1][1].re = buf[ii++]; version(complex_numbers) { c.grad.vel[1][1].im = buf[ii++]; }
                    c.grad.vel[1][2].re = buf[ii++]; version(complex_numbers) { c.grad.vel[1][2].im = buf[ii++]; }
                    c.grad.vel[2][0].re = buf[ii++]; version(complex_numbers) { c.grad.vel[2][0].im = buf[ii++]; }
                    c.grad.vel[2][1].re = buf[ii++]; version(complex_numbers) { c.grad.vel[2][1].im = buf[ii++]; }
                    c.grad.vel[2][2].re = buf[ii++]; version(complex_numbers) { c.grad.vel[2][2].im = buf[ii++]; }
                    // T
                    c.grad.T[0].re = buf[ii++]; version(complex_numbers) { c.grad.T[0].im = buf[ii++]; }
                    c.grad.T[1].re = buf[ii++]; version(complex_numbers) { c.grad.T[1].im = buf[ii++]; }
                    c.grad.T[2].re = buf[ii++]; version(complex_numbers) { c.grad.T[2].im = buf[ii++]; }
                    // formerly tke, omega
                    version(turbulence) {
                        foreach (j; 0 .. nturb) {
                            c.grad.turb[j][0].re = buf[ii++]; version(complex_numbers) { c.grad.turb[j][0].im = buf[ii++]; }
                            c.grad.turb[j][1].re = buf[ii++]; version(complex_numbers) { c.grad.turb[j][1].im = buf[ii++]; }
                            c.grad.turb[j][2].re = buf[ii++]; version(complex_numbers) { c.grad.turb[j][2].im = buf[ii++]; }
                        }
                    }
                    // multi-species
                    version(multi_species_gas) {
                        foreach (j; 0 .. nspecies) {
                            c.grad.massf[j][0].re = buf[ii++]; version(complex_numbers) { c.grad.massf[j][0].im = buf[ii++]; }
                            c.grad.massf[j][1].re = buf[ii++]; version(complex_numbers) { c.grad.massf[j][1].im = buf[ii++]; }
                            c.grad.massf[j][2].re = buf[ii++]; version(complex_numbers) { c.grad.massf[j][2].im = buf[ii++]; }
                        }
                    }
                    // multi-T
                    version(multi_T_gas) {
                        foreach (j; 0 .. nmodes) {
                            c.grad.T_modes[j][0].re = buf[ii++]; version(complex_numbers) { c.grad.T_modes[j][0].im = buf[ii++]; }
                            c.grad.T_modes[j][1].re = buf[ii++]; version(complex_numbers) { c.grad.T_modes[j][1].im = buf[ii++]; }
                            c.grad.T_modes[j][2].re = buf[ii++]; version(complex_numbers) { c.grad.T_modes[j][2].im = buf[ii++]; }
                        }
                    }
                }
            } else {
                // The other block happens to be in this MPI process so
                // we know that we can just access the cell data directly.
                foreach (i; 0 .. ghost_cells.length) {
                    ghost_cells[i].grad.copy_values_from(*(mapped_cells[i].grad));
                }
            }
        } else { // not mpi_parallel
            // For a single process,
            // we know that we can just access the data directly.
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].grad.copy_values_from(*(mapped_cells[i].grad));
            }
        }
        // Rotate the viscous gradients stored at the cells
        if (reorient_vector_quantities) {
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].rotate_gradients(Rmatrix);
            }
        }
        // Done with copying from source cells.
    } // end exchange_viscous_gradient_phase2()

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	throw new Error("GhostCellFullFaceCopy.apply_for_interface_unstructured_grid() not implemented");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("GhostCellFullFaceCopy.apply_unstructured_grid() not implemented");
    }

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	throw new Error("GhostCellFullFaceCopy.apply_for_interface_structured_grid() not implemented");
    }

    @nogc
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
