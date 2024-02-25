// bc/boundary_condition.d
// Base class for boundary condition objects, for use in Eilmer4
//
// Peter J. 2014-07-20 : first cut.
// RG & PJ  2015-12-03 : Decompose boundary conditions into lists of actions
//

module bc.boundary_condition;

import std.conv;
import std.json;
import std.stdio;
import std.string;

import util.lua;
import util.lua_service;
import gas;
import gas.luagas_model;
import nm.luabbla;
import json_helper;
import geom;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvvertex;
import lmr.fluidfvcell;
import fluidblock;
import sfluidblock;
import fluxcalc;
import ssolidblock;
import solidfvcell;
import solidfvinterface;
import bc.ghost_cell_effect;
import bc.boundary_interface_effect;
import bc.boundary_cell_effect;
import bc.boundary_flux_effect;
import bc.user_defined_effects;
import lua_helper;
import grid_motion;
import grid_motion_udf;
import mass_diffusion;

BoundaryCondition make_BC_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    auto newBC = new BoundaryCondition(blk_id, boundary);
    newBC.label = getJSONstring(jsonData, "label", "");
    newBC.type = getJSONstring(jsonData, "type", "");
    newBC.group = getJSONstring(jsonData, "group", "");
    newBC.is_wall_with_viscous_effects = getJSONbool(jsonData, "is_wall_with_viscous_effects", true);
    newBC.ghost_cell_data_available = getJSONbool(jsonData, "ghost_cell_data_available", true);
    newBC.convective_flux_computed_in_bc = getJSONbool(jsonData, "convective_flux_computed_in_bc", false);
    newBC.is_design_surface = getJSONbool(jsonData, "is_design_surface", false);
    newBC.num_cntrl_pts = getJSONint(jsonData, "num_cntrl_pts", 0);
    // Assemble list of preReconAction effects
    auto preReconActions = jsonData["pre_recon_action"].array;
    foreach ( jsonObj; preReconActions ) {
        newBC.preReconAction ~= make_GCE_from_json(jsonObj, blk_id, boundary);
    }
    auto postConvFluxActions = jsonData["post_conv_flux_action"].array;
    foreach ( jsonObj; postConvFluxActions ) {
        newBC.postConvFluxAction ~= make_BFE_from_json(jsonObj, blk_id, boundary);
    }
    auto preSpatialDerivActionsAtBndryFaces = jsonData["pre_spatial_deriv_action_at_bndry_faces"].array;
    foreach ( jsonObj; preSpatialDerivActionsAtBndryFaces ) {
        newBC.preSpatialDerivActionAtBndryFaces ~= make_BIE_from_json(jsonObj, blk_id, boundary);
    }
    auto preSpatialDerivActionsAtBndryCells = jsonData["pre_spatial_deriv_action_at_bndry_cells"].array;
    foreach ( jsonObj; preSpatialDerivActionsAtBndryCells ) {
        newBC.preSpatialDerivActionAtBndryCells ~= make_BCE_from_json(jsonObj, blk_id, boundary);
    }
    auto postDiffFluxActions = jsonData["post_diff_flux_action"].array;
    foreach ( jsonObj; postDiffFluxActions ) {
        newBC.postDiffFluxAction ~= make_BFE_from_json(jsonObj, blk_id, boundary);
    }

    // In case we're reading an old config file with no field_bc, just continue, silently.
    try {
        newBC.field_bc = jsonData["field_bc"].toString.parseJSON; // Deep copy the field_bc data by parsing it again.
    } catch (JSONException e) {}

    return newBC;
} // end make_BC_from_json()


class BoundaryCondition {
    // Boundary condition is built from composable pieces.
public:
    // Location of the boundary condition.
    FluidBlock blk; // the block to which this BC is applied
    int which_boundary; // identity/index of the relevant boundary
    lua_State* myL; // Lua context per BC for user-defined effects.
    // We may have a label for this specific boundary.
    string label;
    // We have a symbolic name for the type of boundary condition
    // when thinking about the flow problem conceptually.
    // Since BC's are composed of generic arrays of effects to be applied,
    // we need a concise way to distinguish one boundary condition from another.
    string type;
    // Sometimes it is convenient to think of individual boundaries
    // grouped together.
    string group;
    // Nature of the boundary condition that may be checked
    // by other parts of the CFD code.
    bool is_design_surface = false;
    int num_cntrl_pts = 0;
    bool is_wall_with_viscous_effects = true;
    bool ghost_cell_data_available = true;
    bool convective_flux_computed_in_bc = false;
    double emissivity = 0.0;
    FVInterface[] faces;
    FluidFVCell[] ghostcells;
    int[] outsigns;

    // list of vertices along boundary, used for the mesh perturbation stage of the
    // shape sensitivity calculator.
    FVVertex[] vertices;
    version(shape_sensitivity) {
        // Bezier curve parameterisation data objects
        Bezier bezier;
        double[] ts;
        Vector3[] surfacePoints;
        // When forming the block local Jacobian matrices for parallel execution of the shape sensitivity calculator,
        // we need to have copies of the neighbour block cells and interfaces that are effected
        // by perturbations in the parent block. We will reference the objects in these arrays by their global ids.
        FluidFVCell[size_t] neighbour_block_cells;
        FVInterface[] neighbour_block_faces;
    }

    // Working storage required when BC is a connection
    // to solid domain.
    // We need to keep these public because the sub-component BCs are going
    // to poke and prod the data here as needed.
    FluidFVCell[] gasCells;
    SolidFVCell[] solidCells;
    FVInterface[] ifaces;

    // We're going to store the JSONdata for the field boundaries, rather then the boundary objects themselves
    // TODO: In the future, rethink if this is a good idea. (NNG)
    JSONValue field_bc;
private:
    // Working storage for boundary flux derivatives
    FlowState* _Lft, _Rght;

public:
    // Action lists.
    // The BoundaryCondition is called at four stages in a global timestep.
    // Those stages are:
    // 1. pre reconstruction
    // 2. post convective flux evaluation
    // 3. pre spatial derivative estimate
    //    (a) apply a boundary interface effect
    //    (b) apply a boundary cell effect
    // 4. post diffusive flux evaluation
    // Note the object may be called more than 4 times depending
    // on the type of time-stepping used to advance the solution.
    // At each of these stages, a series of effects are applied in order
    // with the end goal to leave the boundary values in an appropriate
    // state. We will call this series of effects an action.
    GhostCellEffect[] preReconAction;
    BoundaryFluxEffect[] postConvFluxAction;
    BoundaryInterfaceEffect[] preSpatialDerivActionAtBndryFaces;
    BoundaryCellEffect[] preSpatialDerivActionAtBndryCells;
    BoundaryFluxEffect[] postDiffFluxAction;
    // In block-marching, we will want to temporarily change the down-stream
    // boundary from connected to ExtrapolateCopy.
    // The following saved-XXX-Action variables are places to save the original
    // boundary-condition actions while the temporary actions are in use.
    GhostCellEffect[] savedPreReconAction;
    BoundaryFluxEffect[] savedPostConvFluxAction;
    BoundaryInterfaceEffect[] savedPreSpatialDerivActionAtBndryFaces;
    BoundaryCellEffect[] savedPreSpatialDerivActionAtBndryCells;
    BoundaryFluxEffect[] savedPostDiffFluxAction;

    this(int id, int boundary, bool isWallWithViscousEffects=true,
         bool ghostCellDataAvailable=true, double _emissivity=0.0)
    {
        blk = cast(FluidBlock) globalBlocks[id];  // pick the relevant block out of the collection
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        which_boundary = boundary;
        type = "";
        group = "";
        is_wall_with_viscous_effects = isWallWithViscousEffects;
        ghost_cell_data_available = ghostCellDataAvailable;
        emissivity = _emissivity;
        auto gm = GlobalConfig.gmodel_master;
        auto nturb = GlobalConfig.turb_model.nturb;
        _Lft = new FlowState(gm, nturb);
        _Rght = new FlowState(gm, nturb);
    }

    void finalize()
    {
        if (myL) {
            lua_close(myL);
            myL = null;
        }
    }

    void post_bc_construction()
    {
        foreach (gce; preReconAction) gce.post_bc_construction();
        foreach (bfe; postConvFluxAction) bfe.post_bc_construction();
        foreach (bie; preSpatialDerivActionAtBndryFaces) bie.post_bc_construction();
        foreach (bce; preSpatialDerivActionAtBndryCells) bce.post_bc_construction();
        foreach (bfe; postDiffFluxAction) bfe.post_bc_construction();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "BoundaryCondition(";
        repr ~= "label= \"" ~ label ~ "\", type= \"" ~ type ~ "\", group= \"" ~ group;
        repr ~= "\", is_wall_with_viscous_effects= " ~ to!string(is_wall_with_viscous_effects);
        repr ~= ", ghost_cell_data_available= " ~ to!string(ghost_cell_data_available);
        repr ~= ", convective_flux_computed_in_bc= " ~ to!string(convective_flux_computed_in_bc);
        repr ~= ", is_design_surface= " ~ to!string(is_design_surface);
        repr ~= ", num_cntrl_pts= " ~ to!string(num_cntrl_pts);
        if ( preReconAction.length > 0 ) {
            repr ~= ", preReconAction=[" ~ to!string(preReconAction[0]);
            foreach (i; 1 .. preReconAction.length) {
                repr ~= ", " ~ to!string(preReconAction[i]);
            }
            repr ~= "]";
        }
        if ( postConvFluxAction.length > 0 ) {
            repr ~= ", postConvFluxAction=[" ~ to!string(postConvFluxAction[0]);
            foreach (i; 1 .. postConvFluxAction.length) {
                repr ~= ", " ~ to!string(postConvFluxAction[i]);
            }
            repr ~= "]";
        }
        if ( preSpatialDerivActionAtBndryFaces.length > 0 ) {
            repr ~= ", preSpatialDerivActionAtBndryFaces=[" ~ to!string(preSpatialDerivActionAtBndryFaces[0]);
            foreach (i; 1 .. preSpatialDerivActionAtBndryFaces.length) {
                repr ~= ", " ~ to!string(preSpatialDerivActionAtBndryFaces[i]);
            }
            repr ~= "]";
        }
        if ( preSpatialDerivActionAtBndryCells.length > 0 ) {
            repr ~= ", preSpatialDerivActionAtBndryCells=[" ~ to!string(preSpatialDerivActionAtBndryCells[0]);
            foreach (i; 1 .. preSpatialDerivActionAtBndryCells.length) {
                repr ~= ", " ~ to!string(preSpatialDerivActionAtBndryCells[i]);
            }
            repr ~= "]";
        }
        if ( postDiffFluxAction.length > 0 ) {
            repr ~= ", postDiffFluxAction=[" ~ to!string(postDiffFluxAction[0]);
            foreach (i; 1 .. postDiffFluxAction.length) {
                repr ~= ", " ~ to!string(postDiffFluxAction[i]);
            }
            repr ~= "]";
        }
        repr ~= ")";
        return to!string(repr);
    } // end toString()

    final void pushExtrapolateCopyAction()
    {
        // First, save original boundary-condition actions.
        savedPreReconAction = preReconAction.dup();
        preReconAction.length = 0;
        savedPostConvFluxAction = postConvFluxAction.dup();
        postConvFluxAction.length = 0;
        savedPreSpatialDerivActionAtBndryFaces = preSpatialDerivActionAtBndryFaces.dup();
        preSpatialDerivActionAtBndryFaces.length = 0;
        savedPreSpatialDerivActionAtBndryCells = preSpatialDerivActionAtBndryCells.dup();
        preSpatialDerivActionAtBndryCells.length = 0;
        savedPostDiffFluxAction = postDiffFluxAction.dup();
        postDiffFluxAction.length = 0;
        // Then, make the one action that we need to have a simple outflow.
        preReconAction ~= new GhostCellExtrapolateCopy(blk.id, which_boundary, 0);
    }

    final void restoreOriginalActions()
    {
        preReconAction = savedPreReconAction.dup();
        savedPreReconAction.length = 0;
        postConvFluxAction = savedPostConvFluxAction.dup();
        savedPostConvFluxAction.length = 0;
        preSpatialDerivActionAtBndryFaces = savedPreSpatialDerivActionAtBndryFaces.dup();
        savedPreSpatialDerivActionAtBndryFaces.length = 0;
        preSpatialDerivActionAtBndryCells = savedPreSpatialDerivActionAtBndryCells.dup();
        savedPreSpatialDerivActionAtBndryCells.length = 0;
        postDiffFluxAction = savedPostDiffFluxAction.dup();
        savedPostDiffFluxAction.length = 0;
    }

    final void applyPreReconAction(double t, int gtl, int ftl)
    {
        foreach ( gce; preReconAction ) gce.apply(t, gtl, ftl);
    }

    final void applyPreReconAction(double t, int gtl, int ftl, FVInterface f)
    {
        foreach ( gce; preReconAction ) gce.apply_for_interface(t, gtl, ftl, f);
    }

    final void applyPostConvFluxAction(double t, int gtl, int ftl)
    {
        foreach ( bfe; postConvFluxAction ) bfe.apply(t, gtl, ftl);
    }

    final void applyPostConvFluxAction(double t, int gtl, int ftl, FVInterface f)
    {
        foreach ( bfe; postConvFluxAction ) bfe.apply_for_interface(t, gtl, ftl, f);
    }

    final void applyPreSpatialDerivActionAtBndryFaces(double t, int gtl, int ftl)
    {
        foreach ( bie; preSpatialDerivActionAtBndryFaces ) bie.apply(t, gtl, ftl);
    }

    final void applyPreSpatialDerivActionAtBndryFaces(double t, int gtl, int ftl, FVInterface f)
    {
        foreach ( bie; preSpatialDerivActionAtBndryFaces ) bie.apply_for_interface(t, gtl, ftl, f);
    }

    final void applyPreSpatialDerivActionAtBndryCells(double t, int gtl, int ftl)
    {
        foreach ( bce; preSpatialDerivActionAtBndryCells ) bce.apply(t, gtl, ftl);
    }

    final void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
        foreach ( bfe; postDiffFluxAction ) bfe.apply(t, gtl, ftl);
    }

    final void applyPostDiffFluxAction(double t, int gtl, int ftl, FVInterface f)
    {
        foreach ( bfe; postDiffFluxAction ) bfe.apply_for_interface(t, gtl, ftl, f);
    }

    // The Lua interpreter for the user-defined boundary condition belongs to
    // the boundary condition.  User-defined GhostCellEffect or
    // BoundaryInterfaceEffect objects may need to initialize it.
    void init_lua_State(string luafname)
    {
        if (GlobalConfig.verbosity_level > 1) {
            writefln("Starting new Lua interpreter in BC: user-file=%s, blk.id=%d, bndry=%d",
                     luafname, blk.id, which_boundary);
        }
        if (myL) {
            writeln("Oops, pointer to interpreter is already non-null.");
        } else {
            myL = luaL_newstate();
        }
        if (!myL) { throw new Error("Could not allocate memory for Lua interpreter."); }
        luaL_openlibs(myL);
        // Top-level, generic data.
        lua_pushinteger(myL, blk.id); lua_setglobal(myL, "blkId");
        lua_pushglobaltable(myL);
        registerGasModel(myL);
        registerBBLA(myL);
        luaL_dostring(myL, "require 'lua_helper'");
        pushObj!(GasModel, GasModelMT)(myL, blk.myConfig.gmodel);
        lua_setglobal(myL, "gmodel");
        lua_pushinteger(myL, blk.myConfig.n_species);
        lua_setglobal(myL, "n_species");
        lua_pushinteger(myL, blk.myConfig.n_modes);
        lua_setglobal(myL, "n_modes");
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(myL, GlobalConfig.userPad, "userPad");
        }
        if (blk.grid_type == Grid_t.structured_grid) {
            // Structured-block-specific data
            auto sblk = cast(SFluidBlock) blk;
            assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
            lua_pushinteger(myL, sblk.nic); lua_setglobal(myL, "nicell");
            lua_pushinteger(myL, sblk.njc); lua_setglobal(myL, "njcell");
            lua_pushinteger(myL, sblk.nkc); lua_setglobal(myL, "nkcell");
            lua_pushinteger(myL, Face.north); lua_setglobal(myL, "north");
            lua_pushinteger(myL, Face.east); lua_setglobal(myL, "east");
            lua_pushinteger(myL, Face.south); lua_setglobal(myL, "south");
            lua_pushinteger(myL, Face.west); lua_setglobal(myL, "west");
            lua_pushinteger(myL, Face.top); lua_setglobal(myL, "top");
            lua_pushinteger(myL, Face.bottom); lua_setglobal(myL, "bottom");
        }
        // Boundary-specific data
        lua_pushinteger(myL, which_boundary); lua_setglobal(myL, "boundaryId");
        lua_pushstring(myL, label.toStringz); lua_setglobal(myL, "boundaryLabel");
        lua_pushstring(myL, type.toStringz); lua_setglobal(myL, "boundaryType");
        lua_pushstring(myL, group.toStringz); lua_setglobal(myL, "boundaryGroup");
        // Although we make the helper functions available within
        // the boundary-condition-specific Lua interpreter, we should use
        // those functions only in the context of the master thread.
        setSampleHelperFunctions(myL);
        setGridMotionHelperFunctions(myL);
        // FIXME: Is this capability depreciated?
        // Give access to diffusion coefficients calculation
        //lua_pushcfunction(myL, &luafn_computeBinaryDiffCoeffs);
        //lua_setglobal(myL, "computeBinaryDiffCoeffs");
        // Finally, do the actual user-supplied file.
        if ( luaL_dofile(myL, luafname.toStringz) != 0 ) {
            luaL_error(myL, "error while loading user-defined b.c. file '%s':\n %s\n",
                       luafname.toStringz, lua_tostring(myL, -1));
        }
    }
} // end class BoundaryCondition
