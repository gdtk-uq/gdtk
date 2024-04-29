// extrapolate_copy.d

module bc.ghost_cell_effect.extrapolate_copy;

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


class GhostCellExtrapolateCopy : GhostCellEffect {
public:
    int xOrder;

    this(int id, int boundary, int x_order)
    {
        super(id, boundary, "ExtrapolateCopy");
        xOrder = x_order;
    }

    override string toString() const
    {
        return "ExtrapolateCopy(x_order=" ~ to!string(xOrder) ~ ")";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        FluidFVCell src_cell, ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        if (bc.outsigns[f.i_bndry] == 1) {
            src_cell = f.left_cell;
            ghost0 = f.right_cell;
        } else {
        src_cell = f.right_cell;
        ghost0 = f.left_cell;
        }
        if (xOrder == 1) {
            throw new Error("Linear extrapolation not implemented.");
        } else {
        // Zero-order extrapolation.
        ghost0.fs.copy_values_from(src_cell.fs);
        }
    } // end apply_unstructured_grid()

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        FluidFVCell src_cell, ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            if (bc.outsigns[i] == 1) {
                src_cell = f.left_cell;
                ghost0 = f.right_cell;
            } else {
                src_cell = f.right_cell;
                ghost0 = f.left_cell;
            }
            if (xOrder == 1) {
                throw new Error("Linear extrapolation not implemented.");
            } else {
                // Zero-order extrapolation.
                ghost0.fs.copy_values_from(src_cell.fs);
            }
        } // end foreach face
    } // end apply_unstructured_grid()

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        if (xOrder == 0) {
            // Fill ghost cells with data from just inside the boundary
            // using zero-order extrapolation (i.e. just copy the data).
            // We assume that this boundary is an outflow boundary.
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                FluidFVCell src_cell, dest_cell;
                if (bc.outsigns[f.i_bndry] == 1) {
                    src_cell = f.left_cells[0];
                    dest_cell = f.right_cells[n];
                } else {
                    src_cell = f.right_cells[0];
                    dest_cell = f.left_cells[n];
                }
                dest_cell.fs.copy_values_from(src_cell.fs);
            }
        } else {
            // Extrapolate FlowState (presumably) from cells 0 and 1 into a destination state in cell 2.
            //    |---c0---|---c1---|---c2---|
            // This extrapolation assumes that cell-spacing is uniform.
            if (blk.n_ghost_cell_layers > 0)  {
                FluidFVCell c0, c1, c2;
                if (bc.outsigns[f.i_bndry] == 1) {
                    c0 = f.left_cells[1]; c1 = f.left_cells[0]; c2 = f.right_cells[0];
                } else {
                    c0 = f.right_cells[1]; c1 = f.right_cells[0]; c2 = f.left_cells[0];
                }
                linearly_extrapolate_flowstate(c0.fs, c1.fs, c2.fs);
                foreach (n; 1 .. blk.n_ghost_cell_layers) {
                    // Shuffle along and do next cell.
                    c0 = c1; c1 = c2;
                    c2 = (bc.outsigns[f.i_bndry] == 1) ? f.right_cells[n] : f.left_cells[n];
                    linearly_extrapolate_flowstate(c0.fs, c1.fs, c2.fs);
                }
            }
        } // end else
    } // end apply_for_interface_structured_grid

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        BoundaryCondition bc = blk.bc[which_boundary];
        if (xOrder == 0) {
            // Fill ghost cells with data from just inside the boundary
            // using zero-order extrapolation (i.e. just copy the data).
            // We assume that this boundary is an outflow boundary.
            foreach (i, f; bc.faces) {
                foreach (n; 0 .. blk.n_ghost_cell_layers) {
                    FluidFVCell src_cell, dest_cell;
                    if (bc.outsigns[i] == 1) {
                        src_cell = f.left_cells[0];
                        dest_cell = f.right_cells[n];
                    } else {
                        src_cell = f.right_cells[0];
                        dest_cell = f.left_cells[n];
                    }
                    dest_cell.fs.copy_values_from(src_cell.fs);
                }
            }
        } else {
            // Extrapolate FlowState (presumably) from cells 0 and 1 into a destination state in cell 2.
            //    |---c0---|---c1---|---c2---|
            // This extrapolation assumes that cell-spacing is uniform.
            foreach (i, f; bc.faces) {
                if (blk.n_ghost_cell_layers == 0) continue;
                FluidFVCell c0, c1, c2;
                if (bc.outsigns[i] == 1) {
                    c0 = f.left_cells[1]; c1 = f.left_cells[0]; c2 = f.right_cells[0];
                } else {
                    c0 = f.right_cells[1]; c1 = f.right_cells[0]; c2 = f.left_cells[0];
                }
                linearly_extrapolate_flowstate(c0.fs, c1.fs, c2.fs);
                foreach (n; 1 .. blk.n_ghost_cell_layers) {
                    // Shuffle along and do next cell.
                    c0 = c1; c1 = c2;
                    c2 = (bc.outsigns[i] == 1) ? f.right_cells[n] : f.left_cells[n];
                    linearly_extrapolate_flowstate(c0.fs, c1.fs, c2.fs);
                }
            } // end foreach face
        } // end else
    } // end apply_structured_grid()

private:
    @nogc
    void linearly_extrapolate_flowstate(FlowState* fs0, FlowState* fs1, FlowState* fs2)
    {
        auto gmodel = blk.myConfig.gmodel;
        // Extrapolate on primitive variables
        fs2.gas.rho = 2.0*fs1.gas.rho - fs0.gas.rho;
        fs2.gas.u = 2.0*fs1.gas.u - fs0.gas.u;
        version(multi_T_gas) {
            size_t nmodes = blk.myConfig.n_modes;
            foreach (i; 0 .. nmodes) { fs2.gas.u_modes[i] = 2.0*fs1.gas.u_modes[i] - fs0.gas.u_modes[i]; }
        }
        version(multi_species_gas) {
            size_t nsp = blk.myConfig.n_species;
            if (nsp > 1) {
                foreach (i; 0 .. nsp) { fs2.gas.massf[i] = 2.0*fs1.gas.massf[i] - fs0.gas.massf[i]; }
                scale_mass_fractions(fs2.gas.massf);
            } else {
                fs2.gas.massf[0] = 1.0;
            }
        }
        gmodel.update_thermo_from_rhou(fs2.gas);
        fs2.vel.set(2.0*fs1.vel.x-fs0.vel.x, 2.0*fs1.vel.y-fs0.vel.y, 2.0*fs1.vel.z-fs0.vel.z) ;
        version(MHD) {
            fs2.B.set(2.0*fs1.B.x-fs0.B.x, 2.0*fs1.B.y-fs0.B.y, 2.0*fs1.B.z-fs0.B.z);
        }
        version(turbulence) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            foreach(it; 0 .. nturb) { fs2.turb[it] = 2.0*fs1.turb[it] - fs0.turb[it]; }
        }
        fs2.mu_t = 2.0*fs1.mu_t - fs0.mu_t;
        fs2.k_t = 2.0*fs1.k_t - fs0.k_t;
    } // end linearly_extrapolate_flow()

} // end class GhostCellExtrapolateCopy
