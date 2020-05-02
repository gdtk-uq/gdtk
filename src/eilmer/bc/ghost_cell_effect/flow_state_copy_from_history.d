// flow_state_copy_from_history.d

module bc.ghost_cell_effect.flow_state_copy_from_history;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;

import geom;
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


class GhostCellFlowStateCopyFromHistory : GhostCellEffect {
public:
    this(int id, int boundary, string fileName)
    {
        super(id, boundary, "flowStateCopyFromHistory");
        fhistory = new FlowHistory(fileName);
    }

    override string toString() const
    {
        return format("flowStateCopyFromHistory(filename=\"%s\")", fhistory.fileName);
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	throw new Error("GhostCellFlowStateCopyFromHistory.apply_for_interface_unstructured_grid() not yet implemented");
    }

    // not @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        FVCell ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            if (bc.outsigns[i] == 1) {
                ghost0 = f.right_cell;
            } else {
                ghost0 = f.left_cell;
            }
            ghost0.fs.copy_values_from(fhistory.get_flowstate(t));
        } // end foreach face
    } // end apply_unstructured_grid()

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        // Fill ghost cells with data from just inside the boundary.
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVInterface dest_face;
        FlowState fstate;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        bool nghost3 = (blk.n_ghost_cell_layers == 3);

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    dest_cell = blk.get_cell(i,j+1,k);
                    fstate = fhistory.get_flowstate(t);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.fs.copy_values_from(fstate);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j+3,k);
                        dest_cell.fs.copy_values_from(fstate);
                    }
                } // end i loop
            } // for k
            break;
        case Face.east:
            i = blk.imax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i+1,j,k);
                    fstate = fhistory.get_flowstate(t);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.fs.copy_values_from(fstate);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i+3,j,k);
                        dest_cell.fs.copy_values_from(fstate);
                    }
                } // end j loop
            } // for k
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    dest_cell = blk.get_cell(i,j-1,k);
                    fstate = fhistory.get_flowstate(t);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.fs.copy_values_from(fstate);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j-3,k);
                        dest_cell.fs.copy_values_from(fstate);
                    }
                } // end i loop
            } // for k
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i-1,j,k);
                    fstate = fhistory.get_flowstate(t);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.fs.copy_values_from(fstate);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i-3,j,k);
                        dest_cell.fs.copy_values_from(fstate);
                    }
                } // end j loop
            } // for k
            break;
        case Face.top:
            k = blk.kmax;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i,j,k+1);
                    fstate = fhistory.get_flowstate(t);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.fs.copy_values_from(fstate);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j,k+3);
                        dest_cell.fs.copy_values_from(fstate);
                    }
                } // end j loop
            } // for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i,j,k-1);
                    fstate = fhistory.get_flowstate(t);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.fs.copy_values_from(fstate);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j,k-3);
                        dest_cell.fs.copy_values_from(fstate);
                    }
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()

private:
    FlowHistory fhistory;

} // end class GhostCellFlowStateCopyFromHistory
