/**
 * solidbc.d
 *
 * Author: Rowan G. and Peter J.
 */

module solidbc;

import std.json;
import std.conv;

import util.lua;
import geom;
import util.json_helper;
import solid_ghost_cell;
import solid_boundary_interface_effect;
import solid_boundary_flux_effect;
import ssolidblock;
import globaldata;
import fvcell;
import solidfvcell;
import solidfvinterface;

SolidBoundaryCondition makeSolidBCFromJson(JSONValue jsonData, int blk_id, int boundary,
                                           size_t nicell, size_t njcell, size_t nkcell)
{
    auto setsFluxDirectly = getJSONbool(jsonData, "sets_flux_directly", false);
    auto newBC = new SolidBoundaryCondition(blk_id, boundary, setsFluxDirectly);
    newBC.label = getJSONstring(jsonData, "label", "");
    newBC.type = getJSONstring(jsonData, "type", "");
    newBC.group = getJSONstring(jsonData, "group", "");
    // Assemble list of preSpatialDerivAction effects
    auto preSpatialDerivActionAtBndryCellsList = jsonData["pre_spatial_deriv_action_at_bndry_cells"].array;
    foreach (jsonObj; preSpatialDerivActionAtBndryCellsList) {
        newBC.preSpatialDerivActionAtBndryCells ~= makeSolidGCEfromJson(jsonObj, blk_id, boundary);
    }
    auto preSpatialDerivActionAtBndryFacesList = jsonData["pre_spatial_deriv_action_at_bndry_faces"].array;
    foreach (jsonObj; preSpatialDerivActionAtBndryFacesList) {
        newBC.preSpatialDerivActionAtBndryFaces ~= makeSolidBIEfromJson(jsonObj, blk_id, boundary);
    }
    // Assemble list of postFluxAction effects
    auto postFluxActionList = jsonData["post_flux_action"].array;
    foreach (jsonObj; postFluxActionList) {
        newBC.postFluxAction ~= makeSolidBFEfromJson(jsonObj, blk_id, boundary);
    }
    return newBC;
}

class SolidBoundaryCondition {
public:
    SSolidBlock blk;
    int whichBoundary;
    lua_State* myL; // Lua context per BC for user-defined effects
    bool setsFluxDirectly;
    // We may have a label for this specific boundary.
    string label;
    // We have a symbolic name for the type of boundary condition
    // when thinking about the flow problem conceptually. 
    string type;
    SolidFVInterface[] faces;
    int[] outsigns;
    // Sometimes it is convenient to think of individual boundaries
    // grouped together.
    string group;
    SolidGhostCellEffect[] preSpatialDerivActionAtBndryCells;
    SolidBoundaryInterfaceEffect[] preSpatialDerivActionAtBndryFaces;
    SolidBoundaryFluxEffect[] postFluxAction;
    // data structures used for the special coupled fluid-solid boundary condition
    FVCell[] gasCells;
    SolidFVCell[] solidCells;
    SolidFVInterface[] ifaces;
    
    this(int blkId, int boundary, bool _setsFluxDirectly)
    {
        blk = cast(SSolidBlock) globalBlocks[blkId];
        assert(blk !is null, "Oops, this should be a SSolidBlock object.");
        whichBoundary = boundary;
        setsFluxDirectly = _setsFluxDirectly;
    }
    ~this()
    {
        if (myL != null) lua_close(myL);
    }
    void postBCconstruction()
    {
        foreach (sgce; preSpatialDerivActionAtBndryCells) sgce.postBCconstruction();
        foreach (bie; preSpatialDerivActionAtBndryFaces) bie.postBCconstruction();
        foreach (bfe; postFluxAction) bfe.postBCconstruction();
    }

    final void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel)
    {
        foreach (sgce; preSpatialDerivActionAtBndryCells) sgce.apply(t, tLevel);
    }

    final void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel, SolidFVInterface f)
    {
        foreach (sgce; preSpatialDerivActionAtBndryCells) sgce.apply_for_interface(t, tLevel, f);
    }

    final void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel)
    {
        foreach (bie; preSpatialDerivActionAtBndryFaces) bie.apply(t, tLevel);
    }

    final void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel, SolidFVInterface f)
    {
        foreach (bie; preSpatialDerivActionAtBndryFaces) bie.apply_for_interface(t, tLevel, f);
    }

    final void applyPostFluxAction(double t, int tLevel)
    {
        foreach (bfe; postFluxAction) bfe.apply(t, tLevel);
    }

    final void applyPostFluxAction(double t, int tLevel, SolidFVInterface f)
    {
        foreach (bfe; postFluxAction) bfe.apply_for_interface(t, tLevel, f);
    }

}
