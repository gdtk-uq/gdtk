// internal_copy_then_reflect.d

module bc.ghost_cell_effect.internal_copy_then_reflect;

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


class GhostCellInternalCopyThenReflect : GhostCellEffect {
public:
    
    this(int id, int boundary)
    {
        super(id, boundary, "InternalCopyThenReflect");
    }

    override string toString() const 
    {
        return "InternalCopyThenReflect()";
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        // In contrast with the structured-grid implementation, below,
        // we use a single source cell because we're not sure that
        // the next cell in the interior list (if there is one) will
        // be positioned as a mirror image of the second ghost cell.
        // If the grid is clustered toward the boundary, the error
        // introduced by this zero-order reconstruction will be mitigated.
        // PJ 2016-04-12
        FVCell src_cell, ghost0;
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            if (bc.outsigns[i] == 1) {
                src_cell = f.left_cell;
                ghost0 = f.right_cell;
            } else {
                src_cell = f.right_cell;
                ghost0 = f.left_cell;
            }
            ghost0.fs.copy_values_from(src_cell.fs);
            reflect_normal_velocity(ghost0.fs, f);
            
            if (blk.myConfig.MHD) {
                reflect_normal_magnetic_field(ghost0.fs, f);
            }
        } // end foreach face
    } // end apply_unstructured_grid()

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVInterface IFace;
        auto copy_opt = CopyDataOption.minimal_flow;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    // ghost cell 1.
                    src_cell = blk.get_cell(i,j,k);
                    IFace = src_cell.iface[Face.north];
                    dest_cell = blk.get_cell(i,j+1,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                    // ghost cell 2.
                    src_cell = blk.get_cell(i,j-1,k);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                } // end i loop
            } // for k
            break;
        case Face.east:
            i = blk.imax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    // ghost cell 1.
                    src_cell = blk.get_cell(i,j,k);
                    IFace = src_cell.iface[Face.east];
                    dest_cell = blk.get_cell(i+1,j,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                    // ghost cell 2.
                    src_cell = blk.get_cell(i-1,j,k);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                } // end j loop
            } // for k
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    // ghost cell 1.
                    src_cell = blk.get_cell(i,j,k);
                    IFace = src_cell.iface[Face.south];
                    dest_cell = blk.get_cell(i,j-1,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                    // ghost cell 2.
                    src_cell = blk.get_cell(i,j+1,k);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                } // end i loop
            } // for k
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    // ghost cell 1.
                    src_cell = blk.get_cell(i,j,k);
                    IFace = src_cell.iface[Face.west];
                    dest_cell = blk.get_cell(i-1,j,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                    // ghost cell 2.
                    src_cell = blk.get_cell(i+1,j,k);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                } // end j loop
            } // for k
            break;
        case Face.top:
            k = blk.kmax;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    // ghost cell 1.
                    src_cell = blk.get_cell(i,j,k);
                    IFace = src_cell.iface[Face.top];
                    dest_cell = blk.get_cell(i,j,k+1);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                    // ghost cell 2.
                    src_cell = blk.get_cell(i,j,k-1);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                } // end j loop
            } // for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    // ghost cell 1.
                    src_cell = blk.get_cell(i,j,k);
                    IFace = src_cell.iface[Face.bottom];
                    dest_cell = blk.get_cell(i,j,k-1);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                    // ghost cell 2.
                    src_cell = blk.get_cell(i,j,k+1);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.copy_values_from(src_cell, copy_opt);
                    reflect_normal_velocity(dest_cell.fs, IFace);
                    if (blk.myConfig.MHD) {
                        reflect_normal_magnetic_field(dest_cell.fs, IFace);
                    }
                } // end j loop
            } // for i
            break;
        } // end switch which_boundary
    } // end apply_structured_grid()
} // end class GhostCellInternalCopyThenReflect
