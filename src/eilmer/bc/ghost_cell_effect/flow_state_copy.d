// flow_state_copy.d

module bc.ghost_cell_effect.flow_state_copy;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;

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


class GhostCellFlowStateCopy : GhostCellEffect {
public:
    this(int id, int boundary, in FlowState _fstate)
    {
        super(id, boundary, "flowStateCopy");
        fstate = new FlowState(_fstate);
    }

    override string toString() const
    {
        return "flowStateCopy(fstate=" ~ to!string(fstate) ~ ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        auto ghost = (bc.outsigns[f.i_bndry] == 1) ? f.right_cell : f.left_cell;
	ghost.fs.copy_values_from(fstate);
        if (blk.omegaz != 0.0) {
            into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
        }
    } // end apply_for_interface_unstructured_grid()

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            auto ghost = (bc.outsigns[i] == 1) ? f.right_cell : f.left_cell;
            ghost.fs.copy_values_from(fstate);
            if (blk.omegaz != 0.0) {
                into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
            }
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            auto ghost = (bc.outsigns[f.i_bndry] == 1) ? f.right_cells[n] : f.left_cells[n];
            ghost.fs.copy_values_from(fstate);
            if (blk.omegaz != 0.0) {
                into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
            }
        }
    } // end apply_for_interface_structured_grid()

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                auto ghost = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
                ghost.fs.copy_values_from(fstate);
                if (blk.omegaz != 0.0) {
                    into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz);
                }
            }
        }
    } // end apply_structured_grid()

private:
    FlowState fstate;

} // end class GhostCellFlowStateCopy
