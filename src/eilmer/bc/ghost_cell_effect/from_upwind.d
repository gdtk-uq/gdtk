// from_upwind.d
// We copy the ghost cell data from either the ambient FlowState or
// from the interior cell, whichever is upwind.
// PJ 2021-04-09
//
// Added the DualState class for Flynn Hack, PJ 2022-06-26.

module bc.ghost_cell_effect.from_upwind;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;

import geom;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvcell;
import fluidblock;
import sfluidblock;
import gas;
import bc;


class GhostCellFromUpwindCopy : GhostCellEffect {
public:
    this(int id, int boundary, in FlowState _fstate)
    {
        super(id, boundary, "fromUpwindCopy");
        fstate = new FlowState(_fstate);
    }

    override string toString() const
    {
        return "fromUpwindCopy(fstate=" ~ to!string(fstate) ~ ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        FVCell inside, ghost;
        bool upwind_is_inside;
        if (bc.outsigns[f.i_bndry] == 1) {
            ghost = f.right_cell;
            inside = f.left_cell;
            upwind_is_inside = dot(inside.fs.vel, f.n) > 0.0;
        } else {
            ghost = f.left_cell;
            inside = f.right_cell;
            upwind_is_inside = dot(inside.fs.vel, f.n) < 0.0;
        }
        if (upwind_is_inside) {
            ghost.fs.copy_values_from(inside.fs);
        } else {
            ghost.fs.copy_values_from(fstate);
            if (blk.omegaz != 0.0) { into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz); }
        }
    } // end apply_for_interface_unstructured_grid()

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            apply_for_interface_unstructured_grid(t, gtl, ftl, f);
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        FVCell inside, ghost;
        bool upwind_is_inside;
        if (bc.outsigns[f.i_bndry] == 1) {
            ghost = f.right_cells[0];
            inside = f.left_cells[0];
            upwind_is_inside = dot(inside.fs.vel, f.n) > 0.0;
        } else {
            ghost = f.left_cells[0];
            inside = f.right_cells[0];
            upwind_is_inside = dot(inside.fs.vel, f.n) < 0.0;
        }
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            ghost = (bc.outsigns[f.i_bndry] == 1) ? f.right_cells[n] : f.left_cells[n];
            if (upwind_is_inside) {
                ghost.fs.copy_values_from(inside.fs);
            } else {
                ghost.fs.copy_values_from(fstate);
                if (blk.omegaz != 0.0) { into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz); }
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
            apply_for_interface_structured_grid(t, gtl, ftl, f);
        }
    } // end apply_structured_grid()

private:
    FlowState fstate;

} // end class GhostCellFromUpwindCopy


class GhostCellFromUpwindCopyDualState : GhostCellEffect {
public:
    this(int id, int boundary, in FlowState fstate1, in FlowState fstate2, Vector3 p, Vector3 n)
    {
        super(id, boundary, "fromUpwindCopyDualState");
        this.fstate1 = new FlowState(fstate1);
        this.fstate2 = new FlowState(fstate2);
        this.p = Vector3(p);
        this.n = Vector3(n);
    }

    override string toString() const
    {
        return "fromUpwindCopyDualState(" ~
            "fstate1=" ~ to!string(fstate1) ~
            ", fstate2=" ~ to!string(fstate2) ~
            ", p=" ~ to!string(p) ~
            ", n=" ~ to!string(n) ~
            ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        FVCell inside, ghost;
        bool upwind_is_inside;
        if (bc.outsigns[f.i_bndry] == 1) {
            ghost = f.right_cell;
            inside = f.left_cell;
            upwind_is_inside = dot(inside.fs.vel, f.n) > 0.0;
        } else {
            ghost = f.left_cell;
            inside = f.right_cell;
            upwind_is_inside = dot(inside.fs.vel, f.n) < 0.0;
        }
        if (upwind_is_inside) {
            ghost.fs.copy_values_from(inside.fs);
        } else {
            compute_blended_flow_state(f);
            ghost.fs.copy_values_from(blended_fstate);
            if (blk.omegaz != 0.0) { into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz); }
        }
    } // end apply_for_interface_unstructured_grid()

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            apply_for_interface_unstructured_grid(t, gtl, ftl, f);
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        FVCell inside, ghost;
        bool upwind_is_inside;
        if (bc.outsigns[f.i_bndry] == 1) {
            ghost = f.right_cells[0];
            inside = f.left_cells[0];
            upwind_is_inside = dot(inside.fs.vel, f.n) > 0.0;
        } else {
            ghost = f.left_cells[0];
            inside = f.right_cells[0];
            upwind_is_inside = dot(inside.fs.vel, f.n) < 0.0;
        }
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            ghost = (bc.outsigns[f.i_bndry] == 1) ? f.right_cells[n] : f.left_cells[n];
            if (upwind_is_inside) {
                ghost.fs.copy_values_from(inside.fs);
            } else {
                compute_blended_flow_state(f);
                ghost.fs.copy_values_from(blended_fstate);
                if (blk.omegaz != 0.0) { into_rotating_frame(ghost.fs.vel, ghost.pos[gtl], blk.omegaz); }
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
            apply_for_interface_structured_grid(t, gtl, ftl, f);
        }
    } // end apply_structured_grid()

private:
    FlowState fstate1, fstate2, blended_fstate;
    Vector3 p, n;

    @nogc
    void compute_blended_flow_state(FVInterface f)
    {
        blended_fstate.copy_values_from(fstate1);  // FIX-ME for the details of the DualState flavour
    }

} // end class GhostCellFromUpwindCopyDualState
