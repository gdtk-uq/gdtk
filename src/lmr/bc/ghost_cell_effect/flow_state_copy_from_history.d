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
import fvinterface;
import lmr.fluidfvcell;
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
        my_fs = FlowState(GlobalConfig.gmodel_master, GlobalConfig.turb_model.nturb);
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
        FluidFVCell ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        fhistory.set_flowstate(my_fs, t, gmodel);
        foreach (i, f; bc.faces) {
            ghost0 = (bc.outsigns[i] == 1) ? f.right_cell : f.left_cell;
            ghost0.fs.copy_values_from(my_fs);
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	throw new Error("GhostCellFlowStateCopyFromHistory.apply_for_interface_structured_grid() not yet implemented");
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        fhistory.set_flowstate(my_fs, t, gmodel);
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
                ghost.fs.copy_values_from(my_fs);
            }
        }
    } // end apply_structured_grid()

private:
    FlowHistory fhistory;
    FlowState my_fs;
} // end class GhostCellFlowStateCopyFromHistory
