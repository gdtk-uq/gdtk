// synthesise_flow_state.d
// 2021-07-20 PJ
// This is Lachlan's boundary condition for a applying an inflow with synthetic disturbances.
// There is a companion class BIE_SynthesiseFlowState in the module boundary_interface_effect.d
// and the systhesis function is over in the SyntheticFlowState class in the flowstate.d module.


module bc.ghost_cell_effect.synthesise_flow_state;

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
import fluidblock;
import sfluidblock;
import gas;
import bc;


class GhostCellSynthesiseFlowState : GhostCellEffect {
public:
    this(int id, int boundary, string fileName)
    {
        super(id, boundary, "synthesiseFlowState");
        sfs = new SyntheticFlowState(fileName);
    }

    override string toString() const
    {
        return format("synthesiseFlowState(filename=\"%s\")", sfs.fileName);
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	throw new Error("GhostCellSynthesiseFlowState.apply_for_interface_unstructured_grid() not yet implemented");
    }

    // not @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        foreach (i, f; bc.faces) {
            auto ghost = (bc.outsigns[i] == 1) ? f.right_cell : f.left_cell;
            sfs.set_flowstate(*(ghost.fs), t, ghost.pos[0], gmodel);
        }
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
	throw new Error("GhostCellSynthesiseFlowState.apply_for_interface_structured_grid() not yet implemented");
    }

    // not @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
                sfs.set_flowstate(*(ghost.fs), t, ghost.pos[0], gmodel);
            }
        }
    } // end apply_structured_grid()

private:
    SyntheticFlowState sfs;
} // end class GhostCellSynthesiseFlowState
