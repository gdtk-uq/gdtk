/**
 * boundary_flux_effct.d
 *
 * Authors: RG and PJ
 * Date: 2015-05-07
 **/

module boundary_flux_effect;

import std.stdio;
import std.json;
import std.string;
import std.conv;

import globalconfig;
import globaldata;
import block;
import sblock;
import geom;
import json_helper;
import fvcore;
import fvcell;
import fvinterface;
import solidfvcell;
import solidfvinterface;
import gas_solid_interface;

BoundaryFluxEffect make_BFE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bfeType = jsonData["type"].str;
    BoundaryFluxEffect newBFE;
    
    switch ( bfeType ) {
    case "energy_flux_from_adjacent_solid":
	int otherBlock = getJSONint(jsonData, "other_block", -1);
	string otherFaceName = getJSONstring(jsonData, "other_face", "none");
	int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
	newBFE = new BFE_EnergyFluxFromAdjacentSolid(blk_id, boundary,
						     otherBlock, face_index(otherFaceName),
						     neighbourOrientation);
	break;
    default:
	string errMsg = format("ERROR: The BoundaryFluxEffect type: '%s' is unknown.", bfeType);
	throw new Error(errMsg);
    }
    
    return newBFE;
}

class BoundaryFluxEffect {
public:
    SBlock blk;
    int which_boundary;
    string type;
    
    this(int id, int boundary, string _type)
    {
	blk = gasBlocks[id];
	which_boundary = boundary;
	type = _type;
    }
    override string toString() const
    {
	return "BoundaryFluxEffect()";
    }
    abstract void apply(double t, int gtl, int ftl);
} // end class BoundaryFluxEffect()

// NOTE: This GAS DOMAIN boundary effect has a large
//       and important side-effect:
//       IT ALSO SETS THE FLUX IN THE ADJACENT SOLID DOMAIN
//       AT THE TIME IT IS CALLED.

class BFE_EnergyFluxFromAdjacentSolid : BoundaryFluxEffect {
public:
    int neighbourSolidBlk;
    int neighbourSolidFace;
    int neighbourOrientation;

    this(int id, int boundary,
	 int otherBlock, int otherFace, int orient)
    {
	super(id, boundary, "EnergyFluxFromAdjacentSolid");
	neighbourSolidBlk = otherBlock;
	neighbourSolidFace = otherFace;
	neighbourOrientation = orient;
    }

    override string toString() const 
    {
	return "BFE_EnergyFluxFromAdjacentSolid()";
    }
    
    override void apply(double t, int gtl, int ftl)
    {
	double kS = solidBlocks[neighbourSolidBlk].sp.k;
	computeFluxesAndTemperatures(ftl, kS,
				     _gasCells, _gasIFaces,
				     _solidCells, _solidIFaces);
    }

private:
    // Some private working arrays.
    // We'll pack data into these can pass out
    // to a routine that can compute the flux and
    // temperatures that balance at the interface.
    FVCell[] _gasCells;
    FVInterface[] _gasIFaces;
    SolidFVCell[] _solidCells;
    SolidFVInterface[] _solidIFaces;

public:
    void initSolidCellsAndIFaces()
    {
	size_t i, j, k;
	auto blk = solidBlocks[neighbourSolidBlk];
	switch ( neighbourSolidFace ) {
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    _solidCells ~= blk.getCell(i, j, k);
		    _solidIFaces ~= _solidCells[$-1].iface[Face.south];
		}
	    }
	    break;
	default:
	    throw new Error("initSolidCellsAndIFaces() only implemented for SOUTH face.");
	}
    }

    void initGasCellsAndIFaces()
    {
	size_t i, j, k;
	switch ( which_boundary ) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    _gasCells ~= blk.get_cell(i, j, k);
		    _gasIFaces ~= _gasCells[$-1].iface[Face.north];
		}
	    }
	    break;
	default:
	    throw new Error("initGasCellsAndIFaces() only implemented for NORTH gas face.");
	}
    }
}

