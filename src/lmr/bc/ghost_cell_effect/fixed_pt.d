// fixed_pt.d

module bc.ghost_cell_effect.fixed_pt;

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
import lmr.fluidfvcell;
import fluidblock;
import sfluidblock;
import gas;
import bc;


class GhostCellFixedPT : GhostCellEffect {
public:
    double p_outside;
    double T_outside;

    this(int id, int boundary, double p_outside, double T_outside)
    {
        super(id, boundary, "FixedPT");
        this.p_outside = p_outside;
        this.T_outside = T_outside;
    }

    override string toString() const
    {
        return "FixedPT(p_outside=" ~ to!string(p_outside) ~ ", T_outside=" ~ to!string(T_outside) ~")";
    }

    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        FluidFVCell src_cell, ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        if (bc.outsigns[f.i_bndry] == 1) {
            src_cell = f.left_cell;
            ghost0 = f.right_cell;
        } else {
            src_cell = f.right_cell;
            ghost0 = f.left_cell;
        }
        ghost0.fs.copy_values_from(src_cell.fs);
        ghost0.fs.gas.p = p_outside;
        ghost0.fs.gas.T = T_outside;
        foreach(ref elem; ghost0.fs.gas.T_modes) { elem = T_outside; }
        gmodel.update_thermo_from_pT(ghost0.fs.gas);
    } // end apply_for_interface_unstructured_grid()

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        FluidFVCell src_cell, ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        auto gmodel = blk.myConfig.gmodel;
        foreach (i, f; bc.faces) {
            if (bc.outsigns[i] == 1) {
                src_cell = f.left_cell;
                ghost0 = f.right_cell;
            } else {
                src_cell = f.right_cell;
                ghost0 = f.left_cell;
            }
            ghost0.fs.copy_values_from(src_cell.fs);
            ghost0.fs.gas.p = p_outside;
            ghost0.fs.gas.T = T_outside;
            foreach(ref elem; ghost0.fs.gas.T_modes) { elem = T_outside; }
            gmodel.update_thermo_from_pT(ghost0.fs.gas);
        } // end foreach face
    } // end apply_unstructured_grid()

    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        FluidFVCell src_cell, dest_cell;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (n; 0 .. blk.n_ghost_cell_layers) {
            if (bc.outsigns[f.i_bndry] == 1) {
                src_cell = f.left_cells[n];
                dest_cell = f.right_cells[n];
            } else {
                src_cell = f.right_cells[n];
                dest_cell = f.left_cells[n];
            }
            dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
            dest_cell.fs.gas.p = p_outside;
            dest_cell.fs.gas.T = T_outside;
            foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; }
            gmodel.update_thermo_from_pT(dest_cell.fs.gas);
        }
    } // end apply_for_interface_structured_grid()

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        FluidFVCell src_cell, dest_cell;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                if (bc.outsigns[i] == 1) {
                    src_cell = f.left_cells[n];
                    dest_cell = f.right_cells[n];
                } else {
                    src_cell = f.right_cells[n];
                    dest_cell = f.left_cells[n];
                }
                dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                dest_cell.fs.gas.p = p_outside;
                dest_cell.fs.gas.T = T_outside;
                foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; }
                gmodel.update_thermo_from_pT(dest_cell.fs.gas);
            }
        }
    } // end apply_structured_grid()
} // end class GhostCellFixedPT
