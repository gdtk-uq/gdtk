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
import gas;
import json_helper;
import geom;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvcell;
import block;
import sblock;
import ghost_cell_effect;
import boundary_interface_effect;
import boundary_flux_effect;
import user_defined_effects;

BoundaryCondition make_BC_from_json(JSONValue jsonData, int blk_id, int boundary,
				    size_t nicell, size_t njcell, size_t nkcell)
{
    auto newBC = new BoundaryCondition(blk_id, boundary);
    // Assemble list of preReconAction effects
    auto preReconActionList = jsonData["pre_recon_action"].array;
    foreach ( jsonObj; preReconActionList ) {
	newBC.preReconAction ~= make_GCE_from_json(jsonObj, blk_id, boundary);
	newBC.ghost_cell_data_available = true;
    }
    auto preSpatialDerivActionList = jsonData["pre_spatial_deriv_action"].array;
    foreach ( jsonObj; preSpatialDerivActionList ) {
	newBC.preSpatialDerivAction ~= make_BIE_from_json(jsonObj, blk_id, boundary);
    }
    auto postDiffFluxActionList = jsonData["post_diff_flux_action"].array;
    foreach ( jsonObj; postDiffFluxActionList ) {
	newBC.postDiffFluxAction ~= make_BFE_from_json(jsonObj, blk_id, boundary);
    }
    // [TODO] Only need to the post convective flux option now.
    return newBC;
} // end make_BC_from_json()


class BoundaryCondition {
    // Boundary condition is built from composable pieces.
public:
    // Location of the boundary condition.
    SBlock blk; // the block to which this BC is applied
    int which_boundary; // identity/index of the relevant boundary

    // Nature of the boundary condition that may be checked 
    // by other parts of the CFD code.
    // BCCode type_code = BCCode.slip_wall;
    bool is_wall = true;
    bool ghost_cell_data_available = true;
    double emissivity = 0.0;

    this(int id, int boundary, bool isWall=true, bool ghostCellDataAvailable=true, double _emissivity=0.0)
    {
	blk = gasBlocks[id];  // pick the relevant block out of the collection
	which_boundary = boundary;
	is_wall = isWall;
	ghost_cell_data_available = ghostCellDataAvailable;
	emissivity = _emissivity;
    }

    // Action lists.
    // The BoundaryCondition is called at four stages in a global timestep.
    // Those stages are:
    // 1. pre reconstruction
    // 2. post convective flux evaluation
    // 3. pre spatial derivative estimate
    // 4. post diffusive flux evaluation
    // Note the object may be called more than 4 times depending
    // on the type of time-stepping used to advance the solution.
    // At each of these stages, a series of effects are applied in order
    // with the end goal to leave the boundary values in an appropriate
    // state. We will call this series of effects an action.
    GhostCellEffect[] preReconAction;
    //    BoundaryFluxEffect[] postConvFluxAction;
    BoundaryInterfaceEffect[] preSpatialDerivAction;
    BoundaryFluxEffect[] postDiffFluxAction;

    override string toString() const
    {
	char[] repr;
	repr ~= "BoundaryCondition(";
	if ( preReconAction.length > 0 ) {
	    repr ~= "preReconAction=[" ~ to!string(preReconAction[0]);
	    foreach (i; 1 .. preReconAction.length) {
		repr ~= ", " ~ to!string(preReconAction[i]);
	    }
	    repr ~= "]";
	}
	repr ~= ", ";
	if ( preSpatialDerivAction.length > 0 ) {
	    repr ~= "preSpatialDerivAction=[" ~ to!string(preSpatialDerivAction[0]);
	    foreach (i; 1 .. preSpatialDerivAction.length) {
		repr ~= ", " ~ to!string(preSpatialDerivAction[i]);
	    }
	    repr ~= "]";
	}
	repr ~= ", ";
	if ( postDiffFluxAction.length > 0 ) {
	    repr ~= "postDiffFluxAction=[" ~ to!string(postDiffFluxAction[0]);
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
    /*
    final void applyPostConvFluxAction(double t)
    {
	foreach ( bfe; postConvFluxAction ) bfe.apply(t);
    }
    */
    final void applyPreSpatialDerivAction(double t, int gtl, int ftl)
    {
	foreach ( bie; preSpatialDerivAction ) bie.apply(t, gtl, ftl);
    }
    
    final void applyPostDiffFluxAction(double t, int gtl, int ftl)
    {
	foreach ( bfe; postDiffFluxAction ) bfe.apply(t, gtl, ftl);
    }
} // end class BoundaryCondition

