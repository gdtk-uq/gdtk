// fixed_p.d

module bc.ghost_cell_effect.fixed_p;

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


class GhostCellFixedP : GhostCellEffect {
public:
    double p_outside;

    this(int id, int boundary, double p_outside)
    {
        super(id, boundary, "FixedP");
        this.p_outside = p_outside;
    }

    override string toString() const
    {
        return "FixedP(p_outside=" ~ to!string(p_outside) ~ ")";
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
       gmodel.update_thermo_from_pT(ghost0.fs.gas);
    }

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
            gmodel.update_thermo_from_pT(ghost0.fs.gas);
        } // end foreach face
    } // end apply_unstructured_grid()

    @nogc
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
                gmodel.update_thermo_from_pT(dest_cell.fs.gas);
            }
        }
    } // end apply_structured_grid()
} // end class GhostCellFixedP
