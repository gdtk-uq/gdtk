// solid_ghost_cell.d

module solid_ghost_cell;

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

import simcore;
import util.json_helper;
import geom;
import globaldata;
import globalconfig;
import solidfvinterface;
import ssolidblock;
import solidfvcell;
import solidbc;
import solid_full_face_copy;
import solid_gas_full_face_copy;

SolidGhostCellEffect makeSolidGCEfromJson(JSONValue jsonData, int blk_id, int boundary)
{
    string sgceType = jsonData["type"].str;
    // At the point at which we call this function, we may be inside the block-constructor.
    // Don't attempt the use the block-owned gas model.
    auto gmodel = GlobalConfig.gmodel_master; 
    SolidGhostCellEffect newSGCE;
    switch (sgceType) {
    case "solid_full_face_copy":
    case "solid_full_face_exchange_copy": // old name is also allowed
        int otherBlock = getJSONint(jsonData, "otherBlock", -1);
        string otherFaceName = getJSONstring(jsonData, "otherFace", "none");
        int neighbourOrientation = getJSONint(jsonData, "orientation", 0);
        newSGCE = new SolidGCE_SolidGhostCellFullFaceCopy(blk_id, boundary,
                                           otherBlock, face_index(otherFaceName),
                                           neighbourOrientation);
        break;
    case "solid_gas_full_face_copy": 
        int otherBlock = getJSONint(jsonData, "otherBlock", -1);
        string otherFaceName = getJSONstring(jsonData, "otherFace", "none");
        int neighbourOrientation = getJSONint(jsonData, "orientation", 0);
        newSGCE = new GhostCellSolidGasFullFaceCopy(blk_id, boundary,
                                                    otherBlock, face_index(otherFaceName),
                                                    neighbourOrientation);
        break;
    default:
        string errMsg = format("ERROR: The SolidGhostCellEffect type: '%s' is unknown.", sgceType);
        throw new Exception(errMsg);
    }
    return newSGCE;
}

class SolidGhostCellEffect {
public:
    SSolidBlock blk;
    int which_boundary;
    string desc;

    this(int id, int boundary, string description)
    {
        blk = cast(SSolidBlock) globalBlocks[id];
        assert(blk !is null, "Oops, this should be a SSolidBlock object.");
        which_boundary = boundary;
        desc = description;
    }
    // Most ghost cell effects will not need to do anything
    // special after construction.
    void postBCconstruction() {}
    override string toString() const
    {
        return "SolidGhostCellEffect()";
    }
    abstract void apply(double t, int tLevel);
    abstract void apply_for_interface(double t, int tLevel, SolidFVInterface f);
} // end class SolidGhostCellEffect
