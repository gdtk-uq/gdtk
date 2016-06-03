module solid_boundary_flux_effect;

import std.stdio;
import std.json;
import std.string;
import std.format;
import util.lua;

import simcore;
import json_helper;
import geom;
import globaldata;
import solidfvinterface;
import ssolidblock;

SolidBoundaryFluxEffect makeSolidBFEfromJson(JSONValue jsonData, int blk_id, int boundary)
{
    string bfeType = jsonData["type"].str;
    SolidBoundaryFluxEffect newBFE;
    switch (bfeType) {
    case "zero_flux":
	newBFE = new SolidBFE_ZeroFlux(blk_id, boundary);
	break;
    default:
	string errMsg = format("ERROR: The SolidBoundaryFluxEffect type: '%s' is unknown.", bfeType);
	throw new Exception(errMsg);
    }
    return newBFE;
}

class SolidBoundaryFluxEffect {
public:
    SSolidBlock blk;
    int whichBoundary;
    string type;

    this(int id, int boundary, string _type)
    {
	blk = solidBlocks[id];
	whichBoundary = boundary;
	type = _type;
    }

    abstract void apply(double t, int tLevel);
}

class SolidBFE_ZeroFlux : SolidBoundaryFluxEffect {
public:
    this(int id, int boundary)
    {
	super(id, boundary, "ZeroFlux");
    }

    override void apply(double t, int tLevel)
    {
	size_t i, j, k;
	SolidFVInterface IFace;

	final switch (whichBoundary) {
	case Face.north:
	    j = blk.jmax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.flux = 0.0;
		}
	    }
	    break;
	case Face.east:
	    i = blk.imax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.flux = 0.0;
		}
	    }
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.flux = 0.0;
		}
	    }
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.flux = 0.0;
		}
	    }
	    break;
	case Face.top:
	    throw new Error("[TODO] ZeroFlux bc not implemented for TOP face.");
	case Face.bottom:
	    throw new Error("[TODO] ZeroFlux bc not implemented for BOTTOM face.");

	}
    }

}
