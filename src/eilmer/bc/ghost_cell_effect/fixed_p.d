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
import fvcore;
import fvinterface;
import fvcell;
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
	throw new Error("GhostCellFixedP.apply_for_interface_unstructured_grid() not yet implemented");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        FVCell src_cell, ghost0;
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
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        bool nghost3 = (blk.n_ghost_cell_layers == 3);

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i,j+1,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j+3,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell.fs.gas.p = p_outside;
                        gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    }
                } // end i loop
            } // for k
            break;
        case Face.east:
            i = blk.imax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i+1,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i+3,j,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell.fs.gas.p = p_outside;
                        gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    }
                } // end j loop
            } // for k
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i,j-1,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j-3,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell.fs.gas.p = p_outside;
                        gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    }
                } // end i loop
            } // for k
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i-1,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i-3,j,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell.fs.gas.p = p_outside;
                        gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    }
                } // end j loop
            } // for k
            break;
        case Face.top:
            k = blk.kmax;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i,j,k+1);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j,k+3);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell.fs.gas.p = p_outside;
                        gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    }
                } // end j loop
            } // for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i,j,k-1);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    if (nghost3) {
                        dest_cell = blk.get_cell(i,j,k-3);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell.fs.gas.p = p_outside;
                        gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    }
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()
} // end class GhostCellFixedP
