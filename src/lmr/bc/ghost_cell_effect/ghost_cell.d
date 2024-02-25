// ghost_cell.d
//
// RG & PJ 2015-12-03 : first hack
//         2018-01-20 : refactor into a package

module bc.ghost_cell_effect.ghost_cell;

import std.json;
import std.string;
import std.conv;
import std.stdio;
import std.math;
import std.file;
import std.algorithm;
version(mpi_parallel) {
    import mpi;
}

import geom;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fluidblock;
import sfluidblock;
import gas;
import bc;
import bc.ghost_cell_effect.gas_solid_full_face_copy;

@nogc
void reflect_normal_velocity(ref FlowState fs, in FVInterface IFace)
// Reflects normal velocity with respect to the supplied interface.
//
// The process is to rotate the velocity vector into the local frame of
// the interface, negate the normal (local x-component) velocity and
// rotate back to the global frame.
{
    fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.vel.x = -(fs.vel.x);
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
    fs.vel.y = -(fs.vel.y);
    fs.vel.z = -(fs.vel.z);
    fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

@nogc
void reflect_normal_magnetic_field(ref FlowState fs, in FVInterface IFace)
{
    version(MHD) {
        fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        fs.B.x = -(fs.B.x);
        fs.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        // Used for different boundary conditions in the divergence cleaning- not currently active
        /*if (GlobalConfig.divergence_cleaning) {
                fs.psi = fs.psi + GlobalConfig.c_h * (fs.B.x - fs.B.x);
        }*/
    }
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
        auto flowstate = FlowState(jsonData["flowstate"], gmodel);
        double x0 = getJSONdouble(jsonData, "x0", 0.0);
        double y0 = getJSONdouble(jsonData, "y0", 0.0);
        double z0 = getJSONdouble(jsonData, "z0", 0.0);
        double r = getJSONdouble(jsonData, "r", 0.0);
        newGCE = new GhostCellFlowStateCopy(blk_id, boundary, flowstate, x0, y0, z0, r);
        break;
    case "flowstate_copy_from_profile":
        string fname = getJSONstring(jsonData, "filename", "");
        string match = getJSONstring(jsonData, "match", "xyz");
        newGCE = new GhostCellFlowStateCopyFromProfile(blk_id, boundary, fname, match);
        break;
    case "flowstate_copy_from_history":
        string fname = getJSONstring(jsonData, "filename", "");
        newGCE = new GhostCellFlowStateCopyFromHistory(blk_id, boundary, fname);
        break;
    case "synthesise_flowstate":
        string fname = getJSONstring(jsonData, "filename", "");
        newGCE = new GhostCellSynthesiseFlowState(blk_id, boundary, fname);
        break;
    case "extrapolate_copy":
        int xOrder = getJSONint(jsonData, "x_order", 0);
        newGCE = new GhostCellExtrapolateCopy(blk_id, boundary, xOrder);
        break;
    case "from_upwind_copy":
        auto flowstate = FlowState(jsonData["flowstate"], gmodel);
        newGCE = new GhostCellFromUpwindCopy(blk_id, boundary, flowstate);
        break;
    case "from_upwind_copy_dual_state":
        auto flowstate1 = FlowState(jsonData["flowstate1"], gmodel);
        auto flowstate2 = FlowState(jsonData["flowstate2"], gmodel);
        Vector3 p = getJSONVector3(jsonData, "p", Vector3(0.0,0.0,0.0));
        Vector3 n = getJSONVector3(jsonData, "n", Vector3(0.0,1.0,0.0));
        newGCE = new GhostCellFromUpwindCopyDualState(blk_id, boundary, flowstate1, flowstate2, p, n);
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
        auto stagnation_condition = FlowState(jsonData["stagnation_condition"], gmodel);
        string fname = getJSONstring(jsonData, "filename", "");
        string direction_type = getJSONstring(jsonData, "direction_type", "normal");
        double direction_x = getJSONdouble(jsonData, "direction_x", 1.0);
        double direction_y = getJSONdouble(jsonData, "direction_y", 0.0);
        double direction_z = getJSONdouble(jsonData, "direction_z", 0.0);
        double alpha = getJSONdouble(jsonData, "alpha", 0.0);
        double beta = getJSONdouble(jsonData, "beta", 0.0);
        double mass_flux = getJSONdouble(jsonData, "mass_flux", 0.0);
        double relax_factor = getJSONdouble(jsonData, "relax_factor", 0.1);
        newGCE = new GhostCellFromStagnation(blk_id, boundary,
                                             stagnation_condition, fname,
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
        bool symmetric_mapping = getJSONbool(jsonData, "symmetric_mapping", false);
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
                                             cmff, fname, symmetric_mapping,
                                             transform_pos, c0, n, alpha, delta, lmc,
                                             rvq, Rmatrix);
        break;
    case "gas_solid_full_face_copy":
        int otherBlock = getJSONint(jsonData, "otherBlock", -1);
        string otherFaceName = getJSONstring(jsonData, "otherFace", "none");
        int neighbourOrientation = getJSONint(jsonData, "orientation", 0);
        newGCE = new GhostCellGasSolidFullFaceCopy(blk_id, boundary,
                                                   otherBlock, face_index(otherFaceName),neighbourOrientation);
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
    string desc;

    this(int id, int boundary, string description)
    {
        blk = cast(FluidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        which_boundary = boundary;
        desc = description;
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
    void apply_for_interface(double t, int gtl, int ftl, FVInterface f)
    {
        final switch (blk.grid_type) {
        case Grid_t.unstructured_grid:
            apply_for_interface_unstructured_grid(t, gtl, ftl, f);
            break;
        case Grid_t.structured_grid:
            apply_for_interface_structured_grid(t, gtl, ftl, f);
            break;
        }
    }
    abstract void apply_for_interface_unstructured_grid(double t, int gtl, int ftl, FVInterface f);
    abstract void apply_unstructured_grid(double t, int gtl, int ftl);
    abstract void apply_for_interface_structured_grid(double t, int gtl, int ftl, FVInterface f);
    abstract void apply_structured_grid(double t, int gtl, int ftl);
} // end class GhostCellEffect
