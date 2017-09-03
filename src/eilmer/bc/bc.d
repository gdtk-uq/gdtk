// bc/bc.d
// Base class for boundary condition objects, for use in Eilmer4
//
// Peter J. 2014-07-20 : first cut.
// RG & PJ  2015-12-03 : Decompose boundary conditions into lists of actions
//    

module bc;

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
import sgrid;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvcell;
import block;
import sblock;
import fluxcalc;
import ghost_cell_effect;
import boundary_interface_effect;
import boundary_cell_effect;
import boundary_flux_effect;
import user_defined_effects;
import lua_helper;
import grid_motion;

BoundaryCondition make_BC_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    auto newBC = new BoundaryCondition(blk_id, boundary);
    newBC.label = getJSONstring(jsonData, "label", "");
    newBC.type = getJSONstring(jsonData, "type", "");
    newBC.group = getJSONstring(jsonData, "group", "");
    newBC.is_wall = getJSONbool(jsonData, "is_wall", true);
    newBC.ghost_cell_data_available = getJSONbool(jsonData, "ghost_cell_data_available", true);
    newBC.convective_flux_computed_in_bc = getJSONbool(jsonData, "convective_flux_computed_in_bc", false);
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
    return newBC;
} // end make_BC_from_json()


class BoundaryCondition {
    // Boundary condition is built from composable pieces.
public:
    // Location of the boundary condition.
    Block blk; // the block to which this BC is applied
    int which_boundary; // identity/index of the relevant boundary
    lua_State* myL; // Lua context per BC for user-defined effects.
    // We may have a label for this specific boundary.
    string label;
    // We have a symbolic name for the type of boundary condition
    // when thinking about the flow problem conceptually. 
    string type;
    // Sometimes it is convenient to think of individual boundaries
    // grouped together.
    string group;
    // Nature of the boundary condition that may be checked 
    // by other parts of the CFD code.
    bool is_wall = true;
    bool ghost_cell_data_available = true;
    bool convective_flux_computed_in_bc = false;
    double emissivity = 0.0;
    FVInterface[] faces;
    FVCell[] ghostcells;
    int[] outsigns;

private:
    // Working storage for boundary flux derivatives
    FlowState _Lft, _Rght;

public:
    this(int id, int boundary, bool isWall=true, bool ghostCellDataAvailable=true, double _emissivity=0.0)
    {
	blk = gasBlocks[id];  // pick the relevant block out of the collection
	which_boundary = boundary;
	type = "";
	group = "";
	is_wall = isWall;
	ghost_cell_data_available = ghostCellDataAvailable;
	emissivity = _emissivity;
	auto gm = GlobalConfig.gmodel_master;
	_Lft = new FlowState(gm);
	_Rght = new FlowState(gm);
    }
    ~this()
    {
	if (myL != null) lua_close(myL);
    }
    void post_bc_construction()
    {
	foreach (gce; preReconAction) gce.post_bc_construction();
	foreach (bfe; postConvFluxAction) bfe.post_bc_construction();
	foreach (bie; preSpatialDerivActionAtBndryFaces) bie.post_bc_construction();
	foreach (bce; preSpatialDerivActionAtBndryCells) bce.post_bc_construction();
	foreach (bfe; postDiffFluxAction) bfe.post_bc_construction();
    }

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

    override string toString() const
    {
	char[] repr;
	repr ~= "BoundaryCondition(";
	repr ~= "label= \"" ~ label ~ "\", type= \"" ~ type ~ "\", group= \"" ~ group;
	repr ~= "\", is_wall= " ~ to!string(is_wall);
	repr ~= ", ghost_cell_data_available= " ~ to!string(ghost_cell_data_available);
	repr ~= ", convective_flux_computed_in_bc= " ~ to!string(convective_flux_computed_in_bc);
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
    }

    final void applyPreReconAction(double t, int gtl, int ftl)
    {
	foreach ( gce; preReconAction ) gce.apply(t, gtl, ftl);
    }

    final void applyPostConvFluxAction(double t, int gtl, int ftl)
    {
	foreach ( bfe; postConvFluxAction ) bfe.apply(t, gtl, ftl);
    }
    
    final void applyPreSpatialDerivActionAtBndryFaces(double t, int gtl, int ftl)
    {
	foreach ( bie; preSpatialDerivActionAtBndryFaces ) bie.apply(t, gtl, ftl);
    }

    final void applyPreSpatialDerivActionAtBndryCells(double t, int gtl, int ftl)
    {
	foreach ( bce; preSpatialDerivActionAtBndryCells ) bce.apply(t, gtl, ftl);
    }
    
    final void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
	foreach ( bfe; postDiffFluxAction ) bfe.apply(t, gtl, ftl);
    }

    // The Lua interpreter for the user-defined boundary condition belongs to
    // the boundary condition.  User-defined GhostCellEffect or 
    // BoundaryInterfaceEffect objects may need to initialize it. 
    void init_lua_State(string luafname)
    {
	myL = luaL_newstate();
	luaL_openlibs(myL);
	// Top-level, generic data.
	lua_pushinteger(myL, blk.id); lua_setglobal(myL, "blkId");
	registerGasModel(myL, LUA_GLOBALSINDEX);
	registerBBLA(myL);
	luaL_dostring(myL, "require 'lua_helper'");
	pushObj!(GasModel, GasModelMT)(myL, blk.myConfig.gmodel);
	lua_setglobal(myL, "gmodel");
	lua_pushinteger(myL, blk.myConfig.gmodel.n_species);
	lua_setglobal(myL, "n_species");
	lua_pushinteger(myL, blk.myConfig.gmodel.n_modes);
	lua_setglobal(myL, "n_modes");
	// Structured-block-specific data
	lua_pushinteger(myL, blk.nicell); lua_setglobal(myL, "nicell");
	lua_pushinteger(myL, blk.njcell); lua_setglobal(myL, "njcell");
	lua_pushinteger(myL, blk.nkcell); lua_setglobal(myL, "nkcell");
	lua_pushinteger(myL, Face.north); lua_setglobal(myL, "north");
	lua_pushinteger(myL, Face.east); lua_setglobal(myL, "east");
	lua_pushinteger(myL, Face.south); lua_setglobal(myL, "south");
	lua_pushinteger(myL, Face.west); lua_setglobal(myL, "west");
	lua_pushinteger(myL, Face.top); lua_setglobal(myL, "top");
	lua_pushinteger(myL, Face.bottom); lua_setglobal(myL, "bottom");
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
	// Finally, do the actual user-supplied file.
	luaL_dofile(myL, luafname.toStringz);
    }
} // end class BoundaryCondition


