module solid_boundary_interface_effect;

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

SolidBoundaryInterfaceEffect makeSolidBIEfromJson(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    SolidBoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "fixed_temperature":
	double Twall = getJSONdouble(jsonData, "Twall", 300.0);
	newBIE = new SolidBIE_FixedT(blk_id, boundary, Twall);
	break;
    case "user_defined":
	string fname = getJSONstring(jsonData, "filename", "none");
	newBIE = new SolidBIE_UserDefined(blk_id, boundary, fname);
	break;
    default:
	string errMsg = format("ERROR: The SolidBoundaryInterfaceEffect type: '%s' is unknown.", bieType);
	throw new Exception(errMsg);
    }
    return newBIE;
}

class SolidBoundaryInterfaceEffect {
public:
    SSolidBlock blk;
    int whichBoundary;
    string type;
    
    this(int id, int boundary, string _type) {
	blk = solidBlocks[id];
	whichBoundary = boundary;
	type = _type;
    }

    abstract void apply(double t, int tLevel);
}

class SolidBIE_FixedT : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary, double Twall)
    {
	super(id, boundary, "FixedT");
	_Twall = Twall;
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
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.east:
	    i = blk.imax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.top:
	    throw new Error("[TODO] FixedT bc not implemented for TOP face.");
	case Face.bottom:
	    throw new Error("[TODO] FixedT bc not implemented for BOTTOM face.");

	}

    }

private:
    double _Twall;
}

class SolidBIE_UserDefined : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary, string fname)
    {
	super(id, boundary, "UserDefined");
	luaL_dofile(blk.myL, fname.toStringz);
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
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "north");
		}
	    }
	    break;
	case Face.east:
	    i = blk.imax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "east");
		}
	    }
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "south");
		}
	    }
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace, "west");
		}
	    }
	    break;
	case Face.top:
	    throw new Error("[TODO] FixedT bc not implemented for TOP face.");
	case Face.bottom:
	    throw new Error("[TODO] FixedT bc not implemented for BOTTOM face.");

	}
    }

    void callSolidIfaceUDF(double t, int tLevel, size_t i, size_t j, size_t k,
			   SolidFVInterface IFace, string boundaryName)
    {
	auto L = blk.myL;
	lua_getglobal(L, toStringz("solidInterface_"~boundaryName));
	// Set some userful values for the caller in table
	lua_newtable(L);
	lua_pushnumber(L, t); lua_setfield(L, -2, "t");
	lua_pushnumber(L, dt_global); lua_setfield(L, -2, "dt");
	lua_pushinteger(L, step); lua_setfield(L, -2, "timeStep");
	lua_pushinteger(L, tLevel); lua_setfield(L, -2, "timeLevel");
	lua_pushnumber(L, IFace.pos.x); lua_setfield(L, -2, "x");
	lua_pushnumber(L, IFace.pos.y); lua_setfield(L, -2, "y");
	lua_pushnumber(L, IFace.pos.z); lua_setfield(L, -2, "z");
	lua_pushnumber(L, IFace.n.x); lua_setfield(L, -2, "csX");
	lua_pushnumber(L, IFace.n.y); lua_setfield(L, -2, "csY");
	lua_pushnumber(L, IFace.n.z); lua_setfield(L, -2, "csZ");
	lua_pushnumber(L, IFace.t1.x); lua_setfield(L, -2, "csX1");
	lua_pushnumber(L, IFace.t1.y); lua_setfield(L, -2, "csY1");
	lua_pushnumber(L, IFace.t1.z); lua_setfield(L, -2, "csZ1");
	lua_pushnumber(L, IFace.t2.x); lua_setfield(L, -2, "csX2");
	lua_pushnumber(L, IFace.t2.y); lua_setfield(L, -2, "csY2");
	lua_pushnumber(L, IFace.t2.z); lua_setfield(L, -2, "csZ2");
	lua_pushinteger(L, i); lua_setfield(L, -2, "i");
	lua_pushinteger(L, j); lua_setfield(L, -2, "j");
	lua_pushinteger(L, k); lua_setfield(L, -2, "k");

	// Call function and expect back a temperature value.
	int number_args = 1;
	int number_results = 1;
	if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	    luaL_error(L, "error running user user-defined b.c. solidInterface_%s function: %s\n",
		       toStringz(boundaryName), lua_tostring(L, -1));
	}
	
	IFace.T = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }
}
