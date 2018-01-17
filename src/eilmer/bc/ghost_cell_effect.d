// ghost_cell_effect.d
//
// RG & PJ 2015-12-03 : first hack

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.algorithm;

import geom;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
import fvcore;
import fvinterface;
import fvcell;
import fluidblock;
import sfluidblock;
import gas;
import user_defined_effects;
import bc;

@nogc
void reflect_normal_velocity(ref FlowState fs, in FVInterface IFace)
// Reflects normal velocity with respect to the supplied interface.
//
// The process is to rotate the velocity vector into the local frame of
// the interface, negate the normal (local x-component) velocity and
// rotate back to the global frame.
{
    fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.vel.refx = -(fs.vel.x);
    fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

@nogc
void reverse_tangential_velocity(ref FlowState fs, in FVInterface IFace)
// Reverses the tangential velocity with respect to the supplied interface.
//
// The process is to rotate the velocity vector into the local frame of
// the interface, negate the tangential (local y-component and z-component) velocity and
// rotate back to the global frame.
{
    fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.vel.refy = -(fs.vel.y);
    fs.vel.refz = -(fs.vel.z);
    fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

@nogc
void reflect_normal_magnetic_field(ref FlowState fs, in FVInterface IFace)
{
    fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.B.refx = -(fs.B.x);
    fs.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        //Used for different boundary conditions in the divergence cleaning- not currently active
        /*if (GlobalConfig.divergence_cleaning) { 
                fs.psi = fs.psi + GlobalConfig.c_h * (fs.B.refx - fs.B.x);
        }*/

}

GhostCellEffect make_GCE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string gceType = jsonData["type"].str;
    // At the point at which we call this function, we may be inside the block-constructor.
    // Don't attempt the use the block-owned gas model.
    auto gmodel = GlobalConfig.gmodel_master; 
    GhostCellEffect newGCE;
    switch (gceType) {
    case "internal_copy_then_reflect":
        newGCE = new GhostCellInternalCopyThenReflect(blk_id, boundary);
        break;
    case "flowstate_copy":
        auto flowstate = new FlowState(jsonData["flowstate"], gmodel);
        newGCE = new GhostCellFlowStateCopy(blk_id, boundary, flowstate);
        break;
    case "flowstate_copy_from_profile":
        string fname = getJSONstring(jsonData, "filename", "");
        string match = getJSONstring(jsonData, "match", "xyz");
        newGCE = new GhostCellFlowStateCopyFromProfile(blk_id, boundary, fname, match);
        break;
    case "extrapolate_copy":
        int xOrder = getJSONint(jsonData, "x_order", 0);
        newGCE = new GhostCellExtrapolateCopy(blk_id, boundary, xOrder);
        break;
    case "fixed_pressure":
        double p_outside = getJSONdouble(jsonData, "p_outside", 1.0e5);
        newGCE = new GhostCellFixedP(blk_id, boundary, p_outside);
        break;
    case "fixed_pressure_temperature":
        double p_outside = getJSONdouble(jsonData, "p_outside", 1.0e5);
        double T_outside = getJSONdouble(jsonData, "T_outside", 300.0);
        newGCE = new GhostCellFixedPT(blk_id, boundary, p_outside, T_outside);
        break;
    case "from_stagnation_condition":
        auto stagnation_condition = new FlowState(jsonData["stagnation_condition"], gmodel);
        string direction_type = getJSONstring(jsonData, "direction_type", "normal");
        double direction_x = getJSONdouble(jsonData, "direction_x", 1.0);
        double direction_y = getJSONdouble(jsonData, "direction_y", 0.0);
        double direction_z = getJSONdouble(jsonData, "direction_z", 0.0);
        double alpha = getJSONdouble(jsonData, "alpha", 0.0);
        double beta = getJSONdouble(jsonData, "beta", 0.0);
        double mass_flux = getJSONdouble(jsonData, "mass_flux", 0.0);
        double relax_factor = getJSONdouble(jsonData, "relax_factor", 0.1);
        newGCE = new GhostCellFromStagnation(blk_id, boundary, stagnation_condition,
                                             direction_type, 
                                             Vector3(direction_x, direction_y, direction_z),
                                             alpha, beta,
                                             mass_flux, relax_factor);
        break;
    case "full_face_copy":
    case "full_face_exchange_copy": // old name is also allowed
        int otherBlock = getJSONint(jsonData, "other_block", -1);
        string otherFaceName = getJSONstring(jsonData, "other_face", "none");
        int neighbourOrientation = getJSONint(jsonData, "orientation", 0);
        bool rvq = getJSONbool(jsonData, "reorient_vector_quantities", false);
        double[] Rmatrix = getJSONdoublearray(jsonData, "Rmatrix",
                                              [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]);
        newGCE = new GhostCellFullFaceCopy(blk_id, boundary,
                                           otherBlock, face_index(otherFaceName),
                                           neighbourOrientation,
                                           rvq, Rmatrix);
        break;
    case "mapped_cell_copy":
    case "mapped_cell_exchange_copy": // old name is also allowed
        bool cmff = getJSONbool(jsonData, "cell_mapping_from_file", false);
        string fname = getJSONstring(jsonData, "filename", "none");
        bool transform_pos = getJSONbool(jsonData, "transform_position", false);
        Vector3 c0 = getJSONVector3(jsonData, "c0", Vector3(0.0,0.0,0.0));
        Vector3 n = getJSONVector3(jsonData, "n", Vector3(0.0,0.0,1.0));
        double alpha = getJSONdouble(jsonData, "alpha", 0.0);
        Vector3 delta = getJSONVector3(jsonData, "delta", Vector3(0.0,0.0,0.0));
        bool lmc = getJSONbool(jsonData, "list_mapped_cells", false);
        bool rvq = getJSONbool(jsonData, "reorient_vector_quantities", false);
        double[] Rmatrix = getJSONdoublearray(jsonData, "Rmatrix",
                                              [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0]);
        newGCE = new GhostCellMappedCellCopy(blk_id, boundary,
                                             cmff, fname,
                                             transform_pos, c0, n, alpha, delta, lmc,
                                             rvq, Rmatrix);
        break;
    case "user_defined":
        string fname = getJSONstring(jsonData, "filename", "none");
        newGCE = new UserDefinedGhostCell(blk_id, boundary, fname);
        break;
    default:
        string errMsg = format("ERROR: The GhostCellEffect type: '%s' is unknown.", gceType);
        throw new Exception(errMsg);
    }
    return newGCE;
}

class GhostCellEffect {
public:
    FluidBlock blk;
    int which_boundary;
    string type;

    this(int id, int boundary, string _type)
    {
        blk = globalFluidBlocks[id];
        which_boundary = boundary;
        type = _type;
    }
    // Most ghost cell effects will not need to do anything
    // special after construction.
    // However, the user-defined ghost cells bc need some
    // extra work done to set-up the Lua_state after all
    // of the blocks and bcs have been constructed.
    void post_bc_construction() {}
    override string toString() const
    {
        return "GhostCellEffect()";
    }
    void apply(double t, int gtl, int ftl)
    {
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid: 
            apply_unstructured_grid(t, gtl, ftl);
            break;
        case Grid_t.structured_grid:
            apply_structured_grid(t, gtl, ftl);
        }
    }
    abstract void apply_unstructured_grid(double t, int gtl, int ftl);
    abstract void apply_structured_grid(double t, int gtl, int ftl);
    abstract ref FVCell get_mapped_cell(size_t i);
} // end class GhostCellEffect

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

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

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

            // TO_DO: remove this ugly hack after we have removed ghost_cell_data dependence for
            //        wall-type boundary conditions. K.D. 17/06/2016
            if (bc.type == "wall_no_slip_fixed_t" || bc.type == "wall_no_slip_adiabatic") {
                reverse_tangential_velocity(ghost0.fs, f);
            }
            
            if (blk.myConfig.MHD) {
                reflect_normal_magnetic_field(ghost0.fs, f);
            }
        } // end foreach face
    } // end apply_unstructured_grid()

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

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

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
            ghost0.fs.copy_values_from(fstate);
        } // end foreach face
    } // end apply_unstructured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        // Fill ghost cells with data from just inside the boundary
        // using zero-order extrapolation (i.e. just copy the data).
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVInterface dest_face;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    dest_cell = blk.get_cell(i,j+1,k);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.fs.copy_values_from(fstate);
                } // end i loop
            } // for k
            break;
        case Face.east:
            i = blk.imax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i+1,j,k);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.fs.copy_values_from(fstate);
                } // end j loop
            } // for k
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    dest_cell = blk.get_cell(i,j-1,k);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.fs.copy_values_from(fstate);
                } // end i loop
            } // for k
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i-1,j,k);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.fs.copy_values_from(fstate);
                } // end j loop
            } // for k
            break;
        case Face.top:
            k = blk.kmax;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i,j,k+1);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.fs.copy_values_from(fstate);
                } // end j loop
            } // for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i,j,k-1);
                    dest_cell.fs.copy_values_from(fstate);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.fs.copy_values_from(fstate);
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()

private:
    FlowState fstate;

} // end class GhostCellFlowStateCopy

class GhostCellFlowStateCopyFromProfile : GhostCellEffect {
public:
    this(int id, int boundary, string fileName, string match)
    {
        super(id, boundary, "flowStateCopyFromProfile");
        fprofile = new FlowProfile(fileName, match);
    }
    
    override string toString() const
    {
        return format("flowStateCopyFromProfile(filename=\"%s\", match=\"%s\")",
                      fprofile.fileName, fprofile.posMatch);
    }

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

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
            ghost0.fs.copy_values_from(fprofile.get_flowstate(ghost0.id, ghost0.pos[0]));
            fprofile.adjust_velocity(ghost0.fs, ghost0.pos[0]);
        } // end foreach face
    } // end apply_unstructured_grid()

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        // Fill ghost cells with data from just inside the boundary
        // using zero-order extrapolation (i.e. just copy the data).
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVInterface dest_face;
        FlowState fstate;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    dest_cell = blk.get_cell(i,j+1,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                    dest_cell = blk.get_cell(i,j+2,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                } // end i loop
            } // for k
            break;
        case Face.east:
            i = blk.imax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i+1,j,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                    dest_cell = blk.get_cell(i+2,j,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                } // end j loop
            } // for k
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    dest_cell = blk.get_cell(i,j-1,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                    dest_cell = blk.get_cell(i,j-2,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                } // end i loop
            } // for k
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i-1,j,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                    dest_cell = blk.get_cell(i-2,j,k);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                } // end j loop
            } // for k
            break;
        case Face.top:
            k = blk.kmax;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i,j,k+1);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                    dest_cell = blk.get_cell(i,j,k+2);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                } // end j loop
            } // for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    dest_cell = blk.get_cell(i,j,k-1);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                    dest_cell = blk.get_cell(i,j,k-2);
                    fstate = fprofile.get_flowstate(dest_cell.id, dest_cell.pos[0]);
                    dest_cell.fs.copy_values_from(fstate);
                    fprofile.adjust_velocity(dest_cell.fs, dest_cell.pos[0]);
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()

private:
    FlowProfile fprofile;

} // end class GhostCellFlowStateCopyFromProfile

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

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
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
            if (xOrder == 1) {
                throw new Error("First order extrapolation not implemented.");
            } else {
                // Zero-order extrapolation.
                ghost0.fs.copy_values_from(src_cell.fs);
            }
        } // end foreach face
    } // end apply_unstructured_grid()
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        // Fill ghost cells with data from just inside the boundary
        // using zero-order extrapolation (i.e. just copy the data).
        // We assume that this boundary is an outflow boundary.
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVCell cell_1, cell_2;
        auto gmodel = blk.myConfig.gmodel;
        size_t nsp = gmodel.n_species;
        size_t nmodes = gmodel.n_modes;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    if ( xOrder == 1 ) {
                        //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
                        //      (j-1)        (j)           (j+1)
                        //  dest: ghost cell 1
                        //  [1]: first interior cell
                        //  [2]: second interior cell
                        // This extrapolation assumes that cell-spacing between
                        // cells 1 and 2 continues on in the exterior
                        cell_1 = blk.get_cell(i,j,k);
                        cell_2 = blk.get_cell(i,j-1,k);
                        dest_cell = blk.get_cell(i,j+1,k);
                        // Extrapolate on primitive variables
                        // 1. First exterior point
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                        // 2. Second exterior point
                        //  |---[2]---|||---[1]---|---[dest]------
                        //      (j)        (j+1)       (j+2)
                        cell_2 = cell_1;
                        cell_1 = dest_cell;
                        dest_cell = blk.get_cell(i,j+2,k);
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                    }
                    else {
                        // Zero-order extrapolation
                        src_cell = blk.get_cell(i,j,k);
                        dest_cell = blk.get_cell(i,j+1,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell = blk.get_cell(i,j+2,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    } 
                } // end i loop
            } // for k
            break;
        case Face.east:
            i = blk.imax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    if ( xOrder == 1 ) {
                        //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
                        //      (i-1)        (i)           (i+1)
                        //  dest: ghost cell 1
                        //  [1]: first interior cell
                        //  [2]: second interior cell
                        // This extrapolation assumes that cell-spacing between
                        // cells 1 and 2 continues on in the exterior
                        cell_1 = blk.get_cell(i,j,k);
                        cell_2 = blk.get_cell(i-1,j,k);
                        dest_cell = blk.get_cell(i+1,j,k);
                        // Extrapolate on primitive variables
                        // 1. First exterior point
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                        // 2. Second exterior point
                        //  |---[2]---|||---[1]---|---[dest]------
                        //      (i)        (i+1)       (i+2)
                        cell_2 = cell_1;
                        cell_1 = dest_cell;
                        dest_cell = blk.get_cell(i+2,j,k);
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                    }
                    else {
                        src_cell = blk.get_cell(i,j,k);
                        dest_cell = blk.get_cell(i+1,j,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell = blk.get_cell(i+2,j,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    }
                } // end j loop
            } // for k
            break;
        case Face.south:
            j = blk.jmin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    if ( xOrder == 1 ) {
                        //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
                        //      (j+1)        (j)           (j-1)
                        //  dest: ghost cell 1
                        //  [1]: first interior cell
                        //  [2]: second interior cell
                        // This extrapolation assumes that cell-spacing between
                        // cells 1 and 2 continues on in the exterior
                        cell_1 = blk.get_cell(i,j,k);
                        cell_2 = blk.get_cell(i,j+1,k);
                        dest_cell = blk.get_cell(i,j-1,k);
                        // Extrapolate on primitive variables
                        // 1. First exterior point
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                        // 2. Second exterior point
                        //  |---[2]---|||---[1]---|---[dest]------
                        //      (j)        (j-1)       (j-2)
                        cell_2 = cell_1;
                        cell_1 = dest_cell;
                        dest_cell = blk.get_cell(i,j-2,k);
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                    } else {
                        src_cell = blk.get_cell(i,j,k);
                        dest_cell = blk.get_cell(i,j-1,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell = blk.get_cell(i,j-2,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    }
                } // end i loop
            } // for k
            break;
        case Face.west:
            i = blk.imin;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    if ( xOrder == 1 ) {
                        //  ---[ghost cell 2]---|--- [dest] ---|||--- [1] ---|---[2]----
                        //      (i-2)                 (i-1)           (i)       (i+1)
                        //  dest: ghost cell 1
                        //  [1]: first interior cell
                        //  [2]: second interior cell
                        // This extrapolation assumes that cell-spacing between
                        // cells 1 and 2 continues on in the exterior
                        cell_1 = blk.get_cell(i,j,k);
                        cell_2 = blk.get_cell(i+1,j,k);
                        dest_cell = blk.get_cell(i-1,j,k);
                        // Extrapolate on primitive variables
                        // 1. First exterior point
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                        // 2. Second exterior point
                        //  |---[dest]---|---[1]---|||---[2]---|------|
                        //       (i-2)       (i-1)       (i)
                        cell_2 = cell_1;
                        cell_1 = dest_cell;
                        dest_cell = blk.get_cell(i-2,j,k);
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                    } else {
                        // Zero-order extrapolation
                        src_cell = blk.get_cell(i,j,k);
                        dest_cell = blk.get_cell(i-1,j,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell = blk.get_cell(i-2,j,k);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    }
                } // end j loop
            } // for k
            break;
        case Face.top:
            k = blk.kmax;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    if ( xOrder == 1 ) {
                        //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
                        //      (k-1)        (k)           (k+1)
                        //  dest: ghost cell 1
                        //  [1]: first interior cell
                        //  [2]: second interior cell
                        // This extrapolation assumes that cell-spacing between
                        // cells 1 and 2 continues on in the exterior
                        cell_1 = blk.get_cell(i,j,k);
                        cell_2 = blk.get_cell(i,j,k-1);
                        dest_cell = blk.get_cell(i,j,k+1);
                        // Extrapolate on primitive variables
                        // 1. First exterior point
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                        // 2. Second exterior point
                        //  |---[2]---|||---[1]---|---[dest]------
                        //      (k)        (k+1)       (k+2)
                        cell_2 = cell_1;
                        cell_1 = dest_cell;
                        dest_cell = blk.get_cell(i,j,k+2);
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                    } else {
                        // Zero-order extrapolation
                        src_cell = blk.get_cell(i,j,k);
                        dest_cell = blk.get_cell(i,j,k+1);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell = blk.get_cell(i,j,k+2);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    }
                } // end j loop
            } // for i
            break;
        case Face.bottom:
            k = blk.kmin;
            for (i = blk.imin; i <= blk.imax; ++i) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    if ( xOrder == 1 ) {
                        //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
                        //      (k+1)        (k)           (k-1)
                        //  dest: ghost cell 1
                        //  [1]: first interior cell
                        //  [2]: second interior cell
                        // This extrapolation assumes that cell-spacing between
                        // cells 1 and 2 continues on in the exterior
                        cell_1 = blk.get_cell(i,j,k);
                        cell_2 = blk.get_cell(i,j,k+2);
                        dest_cell = blk.get_cell(i,j,k-1);
                        // Extrapolate on primitive variables
                        // 1. First exterior point
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                        // 2. Second exterior point
                        //  |---[2]---|||---[1]---|---[dest]------
                        //      (k)        (k-1)       (k-2)
                        cell_2 = cell_1;
                        cell_1 = dest_cell;
                        dest_cell = blk.get_cell(i,j,k-2);
                        dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
                        dest_cell.fs.gas.u = 2.0*cell_1.fs.gas.u - cell_2.fs.gas.u;
                        for ( size_t imode = 0; imode < nmodes; ++imode ) {
                            dest_cell.fs.gas.u_modes[imode] = 2.0*cell_1.fs.gas.u_modes[imode] - cell_2.fs.gas.u_modes[imode];
                        }
                        if ( nsp > 1 ) {
                            for ( size_t isp = 0; isp < nsp; ++isp ) {
                                dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
                            }
                            scale_mass_fractions(dest_cell.fs.gas.massf);
                        }
                        else {
                            dest_cell.fs.gas.massf[0] = 1.0;
                        }
                        gmodel.update_thermo_from_rhou(dest_cell.fs.gas);
                        dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
                        dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
                        dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
                        dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
                        dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
                        dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
                        dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
                    } else {
                        src_cell = blk.get_cell(i,j,k);
                        dest_cell = blk.get_cell(i,j,k-1);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                        dest_cell = blk.get_cell(i,j,k-2);
                        dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    }
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()
} // end class GhostCellExtrapolateCopy

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

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

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
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

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
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()
} // end class GhostCellFixedP

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

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

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
            ghost0.fs.gas.T = T_outside;
            foreach(ref elem; ghost0.fs.gas.T_modes) { elem = T_outside; } 
            gmodel.update_thermo_from_pT(ghost0.fs.gas);
        } // end foreach face
    } // end apply_unstructured_grid()
    
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    dest_cell = blk.get_cell(i,j+1,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
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
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
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
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; }
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
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
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
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
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
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
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; }
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
                    dest_cell.fs.gas.p = p_outside;
                    dest_cell.fs.gas.T = T_outside;
                    foreach(ref elem; dest_cell.fs.gas.T_modes) { elem = T_outside; } 
                    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
                } // end j loop
            } // for i
            break;
        } // end switch
    } // end apply_structured_grid()
} // end class GhostCellFixedPT

class GhostCellFromStagnation : GhostCellEffect {
public:
    FlowState stagnation_condition;
    string direction_type;
    Vector3 direction_vector;
    double alpha;
    double beta;
    double mass_flux;
    double relax_factor;
private:
    FlowState inflow_condition;
    double stagnation_entropy;
    double stagnation_enthalpy;
    double p0_min;
    double p0_max;

public:    
    this(int id, int boundary, in FlowState stagnation_condition,
         string direction_type, in Vector3 vec, double alpha, double beta,
         double mass_flux, double relax_factor)
    {
        super(id, boundary, "FromStagnation");
        this.stagnation_condition = new FlowState(stagnation_condition);
        this.direction_type = direction_type;
        this.direction_vector = vec;
        this.direction_vector.normalize();
        this.alpha = alpha;
        this.beta = beta;
        this.mass_flux = mass_flux;
        this.relax_factor = relax_factor;
        auto gmodel = blk.myConfig.gmodel;
        stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
        stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
        inflow_condition = new FlowState(stagnation_condition);
        // The following limits are used when adjusting the stagnation pressure
        // to achieve a fixed mass flow per unit area.
        // The initially-specified stagnation condition needs to be a reasonable guess
        // because T0 will be held fixed and p0 adjusted within the following bounds.
        p0_min = 0.1 * stagnation_condition.gas.p;
        p0_max = 10.0 * stagnation_condition.gas.p;
    }

    override string toString() const 
    {
        return "FromStagnation(stagnation_condition=" ~ to!string(stagnation_condition) ~ 
            ", direction_type=" ~ direction_type ~ 
            ", direction_vector=" ~ to!string(direction_vector) ~
            ", alpha=" ~ to!string(alpha) ~ ", beta=" ~ to!string(beta) ~ 
            ", mass_flux=" ~ to!string(mass_flux) ~
            ", relax_factor=" ~ to!string(relax_factor) ~ ")";
    }

    void set_velocity_components(ref Vector3 vel, double speed, ref FVInterface face)
    {
        switch (direction_type) {
        case "uniform":
            // Given a flow direction.
            vel.set(speed*direction_vector.x, speed*direction_vector.y, speed*direction_vector.z);
            break;
        case "axial":
            // Axial-flow through a presumably circular surface.
            // [TODO] 27-Feb-2014 through 16-Oct-2015 PJ:
            // check that this fall-through is OK.
        case "radial":
            // For turbo inflow, beta sets the axial flow.
            double vz = speed * sin(beta);
            // ...and alpha sets the angle of the flow in the plane of rotation.
            double vt = speed * cos(beta) * sin(alpha);
            double vr = -speed * cos(beta) * cos(alpha);
            double x = face.pos.x;
            double y = face.pos.y;
            double rxy = sqrt(x*x + y*y);
            vel.set(vr*x/rxy - vt*y/rxy, vt*x/rxy + vr*y/rxy, vz);
            break;
        case "normal":
        default:
            // The flow direction is into the block along the local face normal.
            final switch (which_boundary) {
            case Face.north:
            case Face.east:
            case Face.top:
                // Outward-facing normal.
                vel.set(-speed*face.n.x, -speed*face.n.y, -speed*face.n.z);
                break;
            case Face.west:
            case Face.south:
            case Face.bottom:
                // Inward-facing normal.
                vel.set(speed*face.n.x, speed*face.n.y, speed*face.n.z);
            }
        }
    } // end set_velocity_components()

    override ref FVCell get_mapped_cell(size_t i)
    {
        assert(0, "not implemented for this ghost_cell_effect");
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("GhostCellFromStagnation.apply_unstructured_grid() not yet implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        size_t i, j, k;
        FVCell src_cell, dest_cell;
        FVInterface face;
        auto gmodel = blk.myConfig.gmodel;
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");

        double p_stag = 0.0;
        double T_stag = 0.0; // temporary

        final switch (which_boundary) {
        case Face.north:
            j = blk.jmax;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.north];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // North faces have unit normals that nominally point outward from the domain, hence '-='
                    rhoUA -= local_rhoA * dot(cell.fs.vel, face.n); // mass flux
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end k loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.north];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j+1,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j+2,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end k loop
            break;
        case Face.east:
            i = blk.imax;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.east];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // East interfaces have normal that points outwards, hence '-='
                    rhoUA -= local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end j loop
            } // end k loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.east];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i+1,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i+2,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end j loop
            } // end k loop
            break;
        case Face.south:
            j = blk.jmin;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.south];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end k loop
            if ( mass_flux > 0.0 && ftl == 0 ) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.south];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j-1,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j-2,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end k loop
            break;
        case Face.west:
            i = blk.imin;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.west];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // West faces have unit normals that nominally point inward to the domain.
                    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end j loop
            } // end k loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (k = blk.kmin; k <= blk.kmax; ++k) {
                for (j = blk.jmin; j <= blk.jmax; ++j) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.west];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i-1,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i-2,j,k);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end j loop
            } // end k loop
            break;
        case Face.top:
            k = blk.kmax;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.top];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    // Top interfaces have normal that points outwards, hence '-='
                    rhoUA -= local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end j loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.top];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j,k+1);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j,k+2);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end j loop
            break;
        case Face.bottom:
            k = blk.kmin;
            // First, estimate the current bulk inflow condition.
            double area = 0.0;
            double rhoUA = 0.0; // current mass_flux through boundary
            double rhovxA = 0.0; // mass-weighted x-velocity
            double rhovyA = 0.0;
            double rhovzA = 0.0;
            double rhoA = 0.0;
            double pA = 0.0;
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    auto cell = blk.get_cell(i,j,k);
                    face = cell.iface[Face.bottom];
                    area += face.area[0];
                    double local_rhoA = cell.fs.gas.rho * face.area[0];
                    rhoA += local_rhoA;
                    rhoUA += local_rhoA * dot(cell.fs.vel, face.n);
                    rhovxA += local_rhoA * cell.fs.vel.x;
                    rhovyA += local_rhoA * cell.fs.vel.y;
                    rhovzA += local_rhoA * cell.fs.vel.z;
                    pA += cell.fs.gas.p * face.area[0];
                } // end i loop
            } // end j loop
            if (mass_flux > 0.0 && ftl == 0) {
                // Adjust the stagnation pressure to better achieve the specified mass flux.
                // Note that we only do this adjustment once, at the start of a
                // multi-level gas-dynamic update.
                double p = pA / area;
                double dp_over_p = relax_factor * 0.5 / (rhoA/area) * 
                    (mass_flux*mass_flux - rhoUA*fabs(rhoUA)/(area*area)) / p;
                double new_p0 = (1.0 + dp_over_p) * stagnation_condition.gas.p;
                new_p0 = fmin(fmax(new_p0, p0_min), p0_max);
                stagnation_condition.gas.p = new_p0;
                gmodel.update_thermo_from_pT(stagnation_condition.gas);
                stagnation_enthalpy = gmodel.enthalpy(stagnation_condition.gas);
                stagnation_entropy = gmodel.entropy(stagnation_condition.gas);
            }
            double bulk_speed = sqrt((rhovxA/rhoA)^^2 + (rhovyA/rhoA)^^2 + (rhovzA/rhoA)^^2);
            if (rhoUA < 0.0) { bulk_speed = 0.0; } // Block any outflow with stagnation condition.
            // Assume an isentropic process from a known total enthalpy.
            double enthalpy = stagnation_enthalpy - 0.5 * bulk_speed^^2;
            gmodel.update_thermo_from_hs(inflow_condition.gas, enthalpy, stagnation_entropy);
            // Now, apply the ghost-cell conditions
            for (j = blk.jmin; j <= blk.jmax; ++j) {
            for (i = blk.imin; i <= blk.imax; ++i) {
                    src_cell = blk.get_cell(i,j,k);
                    face = src_cell.iface[Face.bottom];
                    // Velocity components may vary with position on the block face.
                    set_velocity_components(inflow_condition.vel, bulk_speed, face);
                    dest_cell = blk.get_cell(i,j,k-1);
                    dest_cell.fs.copy_values_from(inflow_condition);
                    dest_cell = blk.get_cell(i,j,k-2);
                    dest_cell.fs.copy_values_from(inflow_condition);
                } // end i loop
            } // end j loop
            break;
        } // end switch
    } // end apply()
} // end class GhostCellFixedStagnationPT

class GhostCellFullFaceCopy : GhostCellEffect {
public:
    FluidBlock neighbourBlock;
    int neighbourFace;
    int neighbourOrientation;
    bool reorient_vector_quantities;
    double[] Rmatrix;
    // For each ghost cell associated with the boundary,
    // we will have a corresponding "mapped cell" from which we will copy
    // the flow conditions.
    FVCell[] ghost_cells;
    FVCell[] mapped_cells;
    size_t[] mapped_cell_ids;

    this(int id, int boundary,
         int otherBlock, int otherFace, int orient,
         bool reorient_vector_quantities,
         ref const(double[]) Rmatrix)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        super(id, boundary, "FullFaceCopy");
        neighbourBlock = globalFluidBlocks[otherBlock];
        neighbourFace = otherFace;
        neighbourOrientation = orient;
        this.reorient_vector_quantities = reorient_vector_quantities;
        this.Rmatrix = Rmatrix.dup();
    }

    override string toString() const
    { 
        string str = "FullFaceCopy(otherBlock=" ~ to!string(neighbourBlock.id) ~ 
            ", otherFace=" ~ to!string(neighbourFace) ~ 
            ", orient=" ~ to!string(neighbourOrientation) ~
            ", reorient_vector_quantities=" ~ to!string(reorient_vector_quantities) ~
            ", Rmatrix=[";
        foreach(i, v; Rmatrix) {
            str ~= to!string(v);
            str ~= (i < Rmatrix.length-1) ? ", " : "]";
        }
        str ~= ")";
        return str;
    }

    void set_up_cell_mapping()
    {
        auto dest_blk = cast(SFluidBlock) blk;
        assert(dest_blk, "Destination FlowBlock must be a structured-grid block.");
        int destination_face = which_boundary;
        auto src_blk = cast(SFluidBlock) neighbourBlock;
        assert(src_blk, "Source FlowBlock must be a structured-grid block.");
        int src_face = neighbourFace;
        int src_orientation = neighbourOrientation;
        //
        size_t i_dest, i_src, j_dest, j_src, k_dest, k_src;
        FVCell dest0, dest1;
        size_t src0id, src1id;
        //
        if (blk.myConfig.dimensions == 2) {
            // Handle the 2D case separately.
            switch (destination_face) {
            case Face.north:
                j_dest = dest_blk.jmax;  // index of the north-most plane of active cells
                for (size_t i = 0; i < dest_blk.nicell; ++i) {
                    i_dest = i + dest_blk.imin;
                    switch (src_face) {
                    case Face.north:
                        i_src = (src_blk.nicell - i - 1) + src_blk.imin;
                        j_src = src_blk.jmax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        j_src = i + src_blk.jmin;
                        i_src = src_blk.imax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        i_src = i + src_blk.imin;
                        j_src = src_blk.jmin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        j_src = (src_blk.njcell - i - 1) + src_blk.jmin;
                        i_src = src_blk.imin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                    dest0 = dest_blk.get_cell(i_dest,j_dest+1);
                    dest1 = dest_blk.get_cell(i_dest,j_dest+2);
                    ghost_cells ~= dest0;
                    ghost_cells ~= dest1;
                    mapped_cell_ids ~= src0id;
                    mapped_cell_ids ~= src1id;
                } // i loop
                break;
            case Face.east:
                i_dest = dest_blk.imax;  // index of the east-most plane of active cells
                for (size_t j = 0; j < dest_blk.njcell; ++j) {
                    j_dest = j + dest_blk.jmin;
                    switch (src_face) {
                    case Face.north:
                        i_src = j + src_blk.imin;
                        j_src = src_blk.jmax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        j_src = (src_blk.njcell - j - 1) + src_blk.jmin;
                        i_src = src_blk.imax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        i_src = (src_blk.nicell - j - 1) + src_blk.imin;
                        j_src = src_blk.jmin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        j_src = j + src_blk.jmin;
                        i_src = src_blk.imin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                    dest0 = dest_blk.get_cell(i_dest+1,j_dest);
                    dest1 = dest_blk.get_cell(i_dest+2,j_dest);
                    ghost_cells ~= dest0;
                    ghost_cells ~= dest1;
                    mapped_cell_ids ~= src0id;
                    mapped_cell_ids ~= src1id;
                } // j loop
                break;
            case Face.south:
                j_dest = dest_blk.jmin;  // index of the south-most plane of active cells
                for (size_t i = 0; i < dest_blk.nicell; ++i) {
                    i_dest = i + dest_blk.imin;
                    switch (src_face) {
                    case Face.north:
                        i_src = i + src_blk.imin;
                        j_src = src_blk.jmax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        j_src = (src_blk.njcell - i - 1) + src_blk.jmin;
                        i_src = src_blk.imax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        i_src = (src_blk.nicell - i - 1) + src_blk.imin;
                        j_src = src_blk.jmin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        j_src = i + src_blk.jmin;
                        i_src = src_blk.imin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                    dest0 = dest_blk.get_cell(i_dest,j_dest-1);
                    dest1 = dest_blk.get_cell(i_dest,j_dest-2);
                    ghost_cells ~= dest0;
                    ghost_cells ~= dest1;
                    mapped_cell_ids ~= src0id;
                    mapped_cell_ids ~= src1id;
                } // i loop
                break;
            case Face.west:
                i_dest = dest_blk.imin;  // index of the west-most plane of active cells
                for (size_t j = 0; j < dest_blk.njcell; ++j) {
                    j_dest = j + dest_blk.jmin;
                    switch (src_face) {
                    case Face.north:
                        i_src = (src_blk.nicell - j - 1) + src_blk.imin;
                        j_src = src_blk.jmax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1);
                        break;
                    case Face.east:
                        j_src = j + src_blk.jmin;
                        i_src = src_blk.imax; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src);
                        break;
                    case Face.south:
                        i_src = j + src_blk.imin;
                        j_src = src_blk.jmin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1);
                        break;
                    case Face.west:
                        j_src = (src_blk.njcell - j - 1) + src_blk.jmin;
                        i_src = src_blk.imin; 
                        src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src);
                        src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src);
                        break;
                    default:
                        assert(false, "Incorrect boundary connection, source face.");
                    } // end switch src_face
                    dest0 = dest_blk.get_cell(i_dest-1,j_dest);
                    dest1 = dest_blk.get_cell(i_dest-2,j_dest);
                    ghost_cells ~= dest0;
                    ghost_cells ~= dest1;
                    mapped_cell_ids ~= src0id;
                    mapped_cell_ids ~= src1id;
                } // j loop
                break;
            default:
                assert(false, "Incorrect boundary connection, destination face.");
            } // end switch destination_face
        } else {
            // presume dimensions == 3
            // Continue on with 3D work...
            final switch (destination_face) {
            case Face.north:
                j_dest = dest_blk.jmax;  // index of the north-most plane of active cells
                for (size_t i = 0; i < dest_blk.nicell; ++i) {
                    i_dest = i + dest_blk.imin;
                    for (size_t k = 0; k < dest_blk.nkcell; ++k) {
                        k_dest = k + dest_blk.kmin;
                        final switch (src_face) {
                        case Face.north:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = k; break;
                            case 1: i_src = k; k_src = i; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmax; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = i; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            i_src = src_blk.imax; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmin; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = k; break;
                            case 1: j_src = k; k_src = i; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            i_src = src_blk.imin; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = src_blk.njcell - i - 1;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmax; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = k; break;
                            case 1: i_src = k; j_src = i; break;
                            case 2: i_src = i; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - i - 1;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmin; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                        dest0 = dest_blk.get_cell(i_dest,j_dest+1,k_dest);
                        dest1 = dest_blk.get_cell(i_dest,j_dest+2,k_dest);
                        ghost_cells ~= dest0;
                        ghost_cells ~= dest1;
                        mapped_cell_ids ~= src0id;
                        mapped_cell_ids ~= src1id;
                    } // k loop
                } // i loop
                break;
            case Face.east:
                i_dest = dest_blk.imax;  // index of the east-most plane of active cells
                for (size_t j = 0; j < dest_blk.njcell; ++j) {
                    j_dest = j + dest_blk.jmin;
                    for (size_t k = 0; k < dest_blk.nkcell; ++k) {
                        k_dest = k + dest_blk.kmin;
                        final switch (src_face) {
                        case Face.north:
                            final switch (src_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = k; k_src = src_blk.nkcell - j - 1; break;
                            case 2: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = j;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmax; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - j - 1; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 2: j_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = j;
                            }
                            i_src = src_blk.imax; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 2: i_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = j;
                            }
                            j_src = src_blk.jmin; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            final switch (src_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = k; k_src = src_blk.nkcell - j - 1; break;
                            case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = j;
                            }
                            i_src = src_blk.imin; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - j - 1; break;
                            case 2: i_src = j; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = j;
                            }
                            k_src = src_blk.kmax; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            final switch (src_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = k; j_src = src_blk.njcell - j - 1; break;
                            case 2: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = j;
                            }
                            k_src = src_blk.kmin; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                        dest0 = dest_blk.get_cell(i_dest+1,j_dest,k_dest);
                        dest1 = dest_blk.get_cell(i_dest+2,j_dest,k_dest);
                        ghost_cells ~= dest0;
                        ghost_cells ~= dest1;
                        mapped_cell_ids ~= src0id;
                        mapped_cell_ids ~= src1id;
                    } // k loop
                } // j loop
                break;
            case Face.south:
                j_dest = dest_blk.jmin;  // index of the south-most plane of active cells
                for (size_t i = 0; i < dest_blk.nicell; ++i) {
                    i_dest = i + dest_blk.imin;
                    for (size_t k = 0; k < dest_blk.nkcell; ++k) {
                        k_dest = k + dest_blk.kmin;
                        final switch (src_face) {
                        case Face.north:
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = k; break;
                            case 1: i_src = k; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = i;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmax; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = i;
                            } // end switch (src_orientation)
                            i_src = src_blk.imax; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = i;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmin; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = k; break;
                            case 1: j_src = k; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = i;
                            } // end switch (src_orientation)
                            i_src = src_blk.imin; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = i; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = i;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmax; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = k; break;
                            case 1: i_src = k; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = i;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmin; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                        dest0 = dest_blk.get_cell(i_dest,j_dest-1,k_dest);
                        dest1 = dest_blk.get_cell(i_dest,j_dest-2,k_dest);
                        ghost_cells ~= dest0;
                        ghost_cells ~= dest1;
                        mapped_cell_ids ~= src0id;
                        mapped_cell_ids ~= src1id;
                    } // k loop
                } // i loop
                break;
            case Face.west:
                i_dest = dest_blk.imin;  // index of the west-most plane of active cells
                for (size_t j = 0; j < dest_blk.njcell; ++j) {
                    j_dest = j + dest_blk.jmin;
                    for (size_t k = 0; k < dest_blk.nkcell; ++k) {
                        k_dest = k + dest_blk.kmin;
                        final switch (src_face) {
                        case Face.north:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; k_src = k; break;
                            case 1: i_src = k; k_src = j; break;
                            case 2: i_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; k_src = src_blk.nkcell - j - 1;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmax; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            final switch (src_orientation) {
                            case 0: j_src = j; k_src = k; break;
                            case 1: j_src = src_blk.njcell - k - 1; k_src = j; break;
                            case 2: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = k; k_src = src_blk.nkcell - j - 1;
                            } // end switch (src_orientation)
                            i_src = src_blk.imax; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            final switch (src_orientation) {
                            case 0: i_src = j; k_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; k_src = j; break;
                            case 2: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - k - 1; break;
                            case 3: i_src = k; k_src = src_blk.nkcell - j - 1;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmin; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - j - 1; k_src = k; break;
                            case 1: j_src = k; k_src = j; break;
                            case 2: j_src = j; k_src = src_blk.nkcell - k - 1; break;
                            case 3: j_src = src_blk.njcell - k - 1; k_src = src_blk.nkcell - j - 1;
                            } // end switch (src_orientation)
                            i_src = src_blk.imin; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            final switch (src_orientation) {
                            case 0: i_src = j; j_src = k; break;
                            case 1: i_src = src_blk.nicell - k - 1; j_src = j; break;
                            case 2: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = k; j_src = src_blk.njcell - j - 1;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmax; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - j - 1; j_src = k; break;
                            case 1: i_src = k; j_src = j; break;
                            case 2: i_src = j; j_src = src_blk.njcell - k - 1; break;
                            case 3: i_src = src_blk.nicell - k - 1; j_src = src_blk.njcell - j - 1;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmin; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                        dest0 = dest_blk.get_cell(i_dest-1,j_dest,k_dest);
                        dest1 = dest_blk.get_cell(i_dest-2,j_dest,k_dest);
                        ghost_cells ~= dest0;
                        ghost_cells ~= dest1;
                        mapped_cell_ids ~= src0id;
                        mapped_cell_ids ~= src1id;
                    } // k loop
                } // j loop
                break;
            case Face.top:
                k_dest = dest_blk.kmax;  // index of the top-most plane of active cells
                for (size_t j = 0; j < dest_blk.njcell; ++j) {
                    j_dest = j + dest_blk.jmin;
                    for (size_t i = 0; i < dest_blk.nicell; ++i) {
                        i_dest = i + dest_blk.imin;
                        final switch (src_face) {
                        case Face.north:
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = j; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; k_src = i;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmax; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = j; break;
                            case 1: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = j; k_src = i;
                            } // end switch (src_orientation)
                            i_src = src_blk.imax; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - i - 1; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = j; k_src = i;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmin; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = j; k_src = src_blk.nkcell - i - 1; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = src_blk.njcell - j - 1; k_src = i;
                            } // end switch (src_orientation)
                            i_src = src_blk.imin; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = i; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = j; j_src = i;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmax; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = j; j_src = src_blk.njcell - i - 1; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; j_src = i;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmin; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch (src_face)
                        dest0 = dest_blk.get_cell(i_dest,j_dest,k_dest+1);
                        dest1 = dest_blk.get_cell(i_dest,j_dest,k_dest+2);
                        ghost_cells ~= dest0;
                        ghost_cells ~= dest1;
                        mapped_cell_ids ~= src0id;
                        mapped_cell_ids ~= src1id;
                    } // i loop
                } // j loop
                break;
            case Face.bottom:
                k_dest = dest_blk.kmin;  // index of the bottom-most plane of active cells
                for (size_t j = 0; j < dest_blk.njcell; ++j) {
                    j_dest = j + dest_blk.jmin;
                    for (size_t i = 0; i < dest_blk.nicell; ++i) {
                        i_dest = i + dest_blk.imin;
                        final switch (src_face) {
                        case Face.north:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; k_src = j; break;
                            case 1: i_src = j; k_src = i; break;
                            case 2: i_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmax; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src-1,k_src);
                            break;
                        case Face.east:
                            final switch (src_orientation) {
                            case 0: j_src = i; k_src = j; break;
                            case 1: j_src = src_blk.njcell - j - 1; k_src = i; break;
                            case 2: j_src = src_blk.njcell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = j; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            i_src = src_blk.imax; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src-1,j_src,k_src);
                            break;
                        case Face.south:
                            final switch (src_orientation) {
                            case 0: i_src = i; k_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; k_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; k_src = src_blk.nkcell - j - 1; break;
                            case 3: i_src = j; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            j_src = src_blk.jmin; 
                            i_src += src_blk.imin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src+1,k_src);
                            break;
                        case Face.west:
                            final switch (src_orientation) {
                            case 0: j_src = src_blk.njcell - i - 1; k_src = j; break;
                            case 1: j_src = j; k_src = i; break;
                            case 2: j_src = i; k_src = src_blk.nkcell - j - 1; break;
                            case 3: j_src = src_blk.njcell - j - 1; k_src = src_blk.nkcell - i - 1;
                            } // end switch (src_orientation)
                            i_src = src_blk.imin; 
                            j_src += src_blk.jmin; k_src += src_blk.kmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src+1,j_src,k_src);
                            break;
                        case Face.top:
                            final switch (src_orientation) {
                            case 0: i_src = i; j_src = j; break;
                            case 1: i_src = src_blk.nicell - j - 1; j_src = i; break;
                            case 2: i_src = src_blk.nicell - i - 1; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = j; j_src = src_blk.njcell - i - 1;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmax; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src-1);
                            break;
                        case Face.bottom:
                            final switch (src_orientation) {
                            case 0: i_src = src_blk.nicell - i - 1; j_src = j; break;
                            case 1: i_src = j; j_src = i; break;
                            case 2: i_src = i; j_src = src_blk.njcell - j - 1; break;
                            case 3: i_src = src_blk.nicell - j - 1; j_src = src_blk.njcell - i - 1;
                            } // end switch (src_orientation)
                            k_src = src_blk.kmin; 
                            i_src += src_blk.imin; j_src += src_blk.jmin;
                            src0id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src);
                            src1id = src_blk.ijk_indices_to_cell_id(i_src,j_src,k_src+1);
                        } // end switch src_face
                        dest0 = dest_blk.get_cell(i_dest,j_dest,k_dest-1);
                        dest1 = dest_blk.get_cell(i_dest,j_dest,k_dest-2);
                        ghost_cells ~= dest0;
                        ghost_cells ~= dest1;
                        mapped_cell_ids ~= src0id;
                        mapped_cell_ids ~= src1id;
                    } // i loop
                } // j loop
            } // end switch destination_face
        } // end if dimensions == ...
        //
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        } else {
            foreach (i; 0 .. ghost_cells.length) {
                mapped_cells ~= src_blk.cells[mapped_cell_ids[i]];
            }
            foreach (i; 0 .. ghost_cells.length) {
                ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            }
        }
    } // end set_up_cell_mapping()
    
    override ref FVCell get_mapped_cell(size_t i)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        return mapped_cells[i];
    }

    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("GhostCellFullFaceCopy.apply_unstructured_grid() not implemented");
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        SFluidBlock nbblk = cast(SFluidBlock) this.neighbourBlock;
        assert(nbblk !is null, "Oops, this should be an SFluidBlock object.");
        foreach (i; 0 .. ghost_cells.length) {
            ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
            ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
            // The following call to encode_conserved is needed for
            // the block-marching process.
            ghost_cells[i].encode_conserved(gtl, ftl, blk.omegaz);
        }
    }
} // end class GhostCellFullFaceCopy

class GhostCellMappedCellCopy : GhostCellEffect {
public:
    // For each ghost cell associated with the boundary,
    // we will have a corresponding "mapped cell" from which we will copy
    // the flow conditions.
    FVCell[] ghost_cells;
    FVCell[] mapped_cells;
    // Parameters for the calculation of the mapped-cell location.
    bool cell_mapping_from_file;
    string mapped_cells_filename;
    bool transform_position;
    Vector3 c0 = Vector3(0.0, 0.0, 0.0); // default origin
    Vector3 n = Vector3(0.0, 0.0, 1.0); // z-axis
    double alpha = 0.0; // rotation angle (radians) about specified axis vector
    Vector3 delta = Vector3(0.0, 0.0, 0.0); // default zero translation
    bool list_mapped_cells;
    // Parameters for the optional rotation of copied vector data.
    bool reorient_vector_quantities;
    double[] Rmatrix;

    this(int id, int boundary,
         bool cell_mapping_from_file,
         string mapped_cells_filename,
         bool transform_pos,
         ref const(Vector3) c0, ref const(Vector3) n, double alpha,
         ref const(Vector3) delta,
         bool list_mapped_cells,
         bool reorient_vector_quantities,
         ref const(double[]) Rmatrix)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        super(id, boundary, "MappedCellCopy");
        this.cell_mapping_from_file = cell_mapping_from_file;
        this.mapped_cells_filename = mapped_cells_filename;
        this.transform_position = transform_pos;
        this.c0 = c0;
        this.n = n; this.n.normalize();
        this.alpha = alpha;
        this.delta = delta;
        this.list_mapped_cells = list_mapped_cells;
        this.reorient_vector_quantities = reorient_vector_quantities;
        this.Rmatrix = Rmatrix.dup();
    }

    override string toString() const
    { 
        string str = "MappedCellCopy(" ~
            "cell_mapping_from_file=" ~ to!string(cell_mapping_from_file) ~
            ", mapped_cells_filename=" ~ to!string(mapped_cells_filename) ~
            ", transform_position=" ~ to!string(transform_position) ~
            ", c0=" ~ to!string(c0) ~ 
            ", n=" ~ to!string(n) ~ 
            ", alpha=" ~ to!string(alpha) ~
            ", delta=" ~ to!string(delta) ~
            ", list_mapped_cells=" ~ to!string(list_mapped_cells) ~
            ", reorient_vector_quantities=" ~ to!string(reorient_vector_quantities) ~
            ", Rmatrix=[";
        foreach(i, v; Rmatrix) {
            str ~= to!string(v);
            str ~= (i < Rmatrix.length-1) ? ", " : "]";
        }
        str ~= ")";
        return str;
    }

    void set_up_cell_mapping()
    {
        if (cell_mapping_from_file)
            set_up_cell_mapping_from_file();
        else
            set_up_cell_mapping_via_search();
    }

    void set_up_cell_mapping_from_file()
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        string makeFaceTag(const size_t[] node_id_list)
        {
            // We make a tag for this face out of the vertex id numbers
            // sorted in ascending order, to be used as a key in an
            // associative array of indices.  Any cycle of the same vertices
            // should define the same face, so we will use this tag as
            // a way to check if we have made a particular face before.
            // Of course, permutations of the same values that are not
            // correct cycles will produce the same tag, however,
            // we'll not worry about that for now because we should only
            // ever be providing correct cycles.
            size_t[] my_id_list = node_id_list.dup();
            sort(my_id_list);
            string tag = "";
            size_t n = my_id_list.length;
            foreach(i; 0 .. n) {
                if (i > 0) { tag ~= "-"; }
                tag ~= format("%d", my_id_list[i]);
            }
            return tag;
        }
        int[2][string] mapped_cells_list; // list of cells to be mapped to ghost cells,
                                          // referenced by the neighbour cells id. 
        // Stage 1 -- read mapped_cells file
        if (!exists(mapped_cells_filename)) {
            assert(0, "mapped_cells file does not exist.");
        } // else if the file exists
        auto f = File(mapped_cells_filename, "r");
        string getHeaderContent(string target)
        // Helper function to proceed through file, line-by-line,
        // looking for a particular header line.
        // Returns the content from the header line and leaves the file
        // at the next line to be read, presumably with expected data.
        {
            while (!f.eof) {
                auto line = f.readln().strip();
                if (canFind(line, target)) {
                    auto tokens = line.split("=");
                    return tokens[1].strip();
                }
            } // end while
            return ""; // didn't find the target
        }
        string blk_id_str = to!string(blk.id);
        string mapped_cells_tag = "NMappedCells in BLOCK[" ~ blk_id_str ~ "]";
        auto nfaces  = to!int(getHeaderContent(mapped_cells_tag));
        foreach(i; 0..nfaces) {
            auto lineContent = f.readln().strip();
            auto tokens = lineContent.split();
            int secondary_blk_id = to!int(tokens[1]);
            int secondary_cell_id = to!int(tokens[2]);
            auto faceTag = tokens[0];
            int[2] mapped_cell;
            mapped_cell[0] = secondary_blk_id;
            mapped_cell[1] = secondary_cell_id;
            mapped_cells_list[faceTag] = mapped_cell;
        }
        // Stage 2 -- map cells
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid: 
            BoundaryCondition bc = blk.bc[which_boundary];
            foreach (i, face; bc.faces) {
                size_t[] my_vtx_list;
                foreach(vtx; face.vtx)
                    my_vtx_list ~= vtx.id;
                string faceTag =  makeFaceTag(my_vtx_list);
                if (bc.outsigns[i] == 1) {
                    ghost_cells ~= face.right_cell;
                    int ghost_cell_blk_id = mapped_cells_list[faceTag][0];
                    int ghost_cell_id = mapped_cells_list[faceTag][1];
                    mapped_cells ~= globalFluidBlocks[ghost_cell_blk_id].cells[ghost_cell_id];
                } else {
                    ghost_cells ~= face.left_cell;
                    int ghost_cell_blk_id = mapped_cells_list[faceTag][0];
                    int ghost_cell_id = mapped_cells_list[faceTag][1];
                    mapped_cells ~= globalFluidBlocks[ghost_cell_blk_id].cells[ghost_cell_id];
                }
                if (list_mapped_cells) {
                    writeln("    ghost-cell-pos=", to!string(ghost_cells[$-1].pos[0]), 
                            " mapped-cell-pos=", to!string(mapped_cells[$-1].pos[0]));
                }
                ghost_cells[$-1].copy_values_from(mapped_cells[$-1], CopyDataOption.grid);
            } // end foreach face
            break;
        case Grid_t.structured_grid:
            throw new Error("mapped cells from file not implemented for structured grids");
        }
    }
    
    void set_up_cell_mapping_via_search()
    {
        // Stage-2 construction for this boundary condition,
        // for the situation when we haven't been told where to find our mapped cells.
        //
        // Needs to be called after the cell geometries have been computed,
        // because the search sifts through the cells in blocks
        // that happen to be in the local process.
        // The search does not extend to cells in blocks in other MPI tasks.
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid: 
            BoundaryCondition bc = blk.bc[which_boundary];
            foreach (i, f; bc.faces) {
                if (bc.outsigns[i] == 1) {
                    ghost_cells ~= f.right_cell;
                } else {
                    ghost_cells ~= f.left_cell;
                }
            } // end foreach face
            break;
        case Grid_t.structured_grid:
            size_t i, j, k;
            auto blk = cast(SFluidBlock) this.blk;
            assert(blk !is null, "Oops, this should be an SFluidBlock object.");
            final switch (which_boundary) {
            case Face.north:
                j = blk.jmax;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        ghost_cells ~= blk.get_cell(i,j+1,k);
                        ghost_cells ~= blk.get_cell(i,j+2,k);
                    } // end i loop
                } // for k
                break;
            case Face.east:
                i = blk.imax;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i+1,j,k);
                        ghost_cells ~= blk.get_cell(i+2,j,k);
                    } // end j loop
                } // for k
                break;
            case Face.south:
                j = blk.jmin;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (i = blk.imin; i <= blk.imax; ++i) {
                        ghost_cells ~= blk.get_cell(i,j-1,k);
                        ghost_cells ~= blk.get_cell(i,j-2,k);
                    } // end i loop
                } // for k
                break;
            case Face.west:
                i = blk.imin;
                for (k = blk.kmin; k <= blk.kmax; ++k) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i-1,j,k);
                        ghost_cells ~= blk.get_cell(i-2,j,k);
                    } // end j loop
                } // for k
                break;
            case Face.top:
                k = blk.kmax;
                for (i = blk.imin; i <= blk.imax; ++i) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i,j,k+1);
                        ghost_cells ~= blk.get_cell(i,j,k+2);
                    } // end j loop
                } // for i
                break;
            case Face.bottom:
                k = blk.kmin;
                for (i = blk.imin; i <= blk.imax; ++i) {
                    for (j = blk.jmin; j <= blk.jmax; ++j) {
                        ghost_cells ~= blk.get_cell(i,j,k-1);
                        ghost_cells ~= blk.get_cell(i,j,k-2);
                    } // end j loop
                } // for i
                break;
            } // end switch
        } // end switch blk.grid_type
        // Now that we have a collection of the local ghost cells,
        // locate the corresponding active cell so that we can later
        // copy that cell's flow state.
        if (list_mapped_cells) {
            writefln("Mapped cells for block[%d] boundary[%d]:", blk.id, which_boundary);
        }
        foreach (mygc; ghost_cells) {
            Vector3 ghostpos = mygc.pos[0];
            Vector3 mypos = ghostpos;
            if (transform_position) {
                Vector3 c1 = c0 + dot(n, (ghostpos - c0)) * n;
                Vector3 t1 = (ghostpos - c1);
                t1.normalize();
                Vector3 t2 = cross(n, t1);
                mypos = c1 + cos(alpha) * t1 + sin(alpha) * t2;
                mypos += delta;
            }
            // Because we need to access all of the gas blocks in the following search,
            // we have to run this set_up_cell_mapping function from a serial loop.
            // In parallel code, threads other than the main thread get uninitialized
            // versions of the localFluidBlocks array.
            //
            // First, attempt to find the enclosing cell at the specified position.
            bool found = false;
            foreach (ib, blk; localFluidBlocks) {
                found = false;
                size_t indx = 0;
                blk.find_enclosing_cell(mypos, indx, found);
                if (found) {
                    mapped_cells ~= blk.cells[indx];
                    break;
                }
            }
            if (!found) {
                // Fall back to nearest cell search.
                FVCell closest_cell = localFluidBlocks[0].cells[0];
                Vector3 cellpos = closest_cell.pos[0];
                Vector3 dp = cellpos - mypos;
                double min_distance = abs(dp);
                foreach (blk; localFluidBlocks) {
                    foreach (cell; blk.cells) {
                        dp = cell.pos[0] - mypos;
                        double distance = abs(dp);
                        if (distance < min_distance) {
                            closest_cell = cell;
                            min_distance = distance;
                        }
                    }
                }
                mapped_cells ~= closest_cell;
            }
            if (list_mapped_cells) {
                writeln("    ghost-cell-pos=", to!string(mygc.pos[0]), 
                        " mapped-cell-pos=", to!string(mapped_cells[$-1].pos[0]));
            }
            mygc.copy_values_from(mapped_cells[$-1], CopyDataOption.grid);
        } // end foreach mygc
    } // end set_up_cell_mapping()

    override ref FVCell get_mapped_cell(size_t i)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        return mapped_cells[i];
    }
    
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        foreach (i; 0 .. ghost_cells.length) {
            ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
            ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
            // [TODO] PJ 2018-01-14 If unstructured blocks ever get used in
            // the block-marching process, we will need a call to encode_conserved
            // at this point.  See the GhostCellFullFaceCopy class.
        }
    }

    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        version(mpi_parallel) {
            // [TODO] PJ 2018-01-17 communication...
        }
        foreach (i; 0 .. ghost_cells.length) {
            ghost_cells[i].fs.copy_values_from(mapped_cells[i].fs);
            ghost_cells[i].copy_values_from(mapped_cells[i], CopyDataOption.grid);
            if (reorient_vector_quantities) {
                ghost_cells[i].fs.reorient_vector_quantities(Rmatrix);
            }
        }
    }
} // end class GhostCellMappedCellCopy

