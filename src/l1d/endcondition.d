// endcondition.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module endcondition;


import std.conv;
import std.stdio;
import std.string;
import std.json;

import json_helper;
import geom;
import gas;
import gasflow;
import config;
import tube;
import gasslug;
import piston;
import simcore;

enum End { L, R };

class EndCondition {
public:
    size_t indx;
    GasSlug slugL = null;
    GasSlug slugR = null;
    End slugL_end, slugR_end;
    Piston pistonL, pistonR;
    End pistonL_face, pistonR_face;
    
    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct end-condition[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        //
        int slug_id = getJSONint(jsonData, "left-slug-id", -1);
        if (slug_id >= 0) { slugL = gasslugs[slug_id]; }
        string endStr = getJSONstring(jsonData, "left-slug-end-id", "L");
        slugL_end = (endStr == "L") ? End.L : End.R;
        //
        slug_id = getJSONint(jsonData, "right-slug-id", -1);
        if (slug_id >= 0) { slugR = gasslugs[slug_id]; }
        endStr = getJSONstring(jsonData, "right-slug-end-id", "L");
        slugR_end = (endStr == "L") ? End.L : End.R;
        //
        int piston_id = getJSONint(jsonData, "left-piston-id", -1);
        if (piston_id >= 0) { pistonL = pistons[piston_id]; }
        endStr = getJSONstring(jsonData, "left-piston-face-id", "L");
        pistonL_face = (endStr == "L") ? End.L : End.R;
        //
        piston_id = getJSONint(jsonData, "right-piston-id", -1);
        if (piston_id >= 0) { pistonR = pistons[piston_id]; }
        endStr = getJSONstring(jsonData, "right-piston-face-id", "L");
        pistonR_face = (endStr == "L") ? End.L : End.R;
        //
        if (L1dConfig.verbosity_level >= 1) {
            writeln("  connections:");
            writeln("    slugL_id= ", ((slugL) ? to!int(slugL.indx) : -1),
                    " end= ", slugL_end);
            writeln("    slugR_id= ", ((slugR) ? to!int(slugR.indx) : -1),
                    " end= ", slugR_end);
            writeln("    pistonL_id= ", ((pistonL) ? to!int(pistonL.indx) : -1),
                    " face= ", pistonL_face);
            writeln("    pistonR_id= ", ((pistonR) ? to!int(pistonR.indx) : -1),
                    " face= ", pistonR_face);
        }        
    } // end constructor
    
} // end class EndCondition


class Diaphragm : EndCondition {
public:
    double p_burst;
    double dt_hold;
    double dxL;
    double dxR;
    
    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Diaphragm[", indx, "]:");
        }
        auto connectionData = jsonData["connections"];
        super(indx, connectionData);
        p_burst = getJSONdouble(jsonData, "p_burst", 0.0);
        dt_hold = getJSONdouble(jsonData, "dt_hold", 0.0);
        dxL = getJSONdouble(jsonData, "dxL", 0.0);
        dxR = getJSONdouble(jsonData, "dxR", 0.0);
        if (L1dConfig.verbosity_level >= 1) {
            writeln("  p_burst= ", p_burst);
            writeln("  dt_hold= ", dt_hold);
            writeln("  dxL= ", dxL);
            writeln("  dxR= ", dxR);
        }
    } // end constructor
    
} // end class Diaphragm


class GasInterface : EndCondition {
public:

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("GasInterface[", indx, "]:");
        }
        auto connectionData = jsonData["connections"];
        super(indx, connectionData);
    } // end constructor
    
} // end class GasInterface


class FreeEnd : EndCondition {
public:

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("FreeEnd[", indx, "]:");
        }
        auto connectionData = jsonData["connections"];
        super(indx, connectionData);
    } // end constructor
    
} // end class FreeEnd


class VelocityEnd : EndCondition {
public:
    double vel;
    
    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("VelocityEnd[", indx, "]:");
        }
        auto connectionData = jsonData["connections"];
        super(indx, connectionData);
        vel = getJSONdouble(jsonData, "vel", 0.0);
        if (L1dConfig.verbosity_level >= 1) {
            writeln("  vel= ", vel);
        }
    } // end constructor
    
} // end class VelocityEnd


class PistonFace : EndCondition {
public:

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("PistonFace[", indx, "]:");
        }
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct piston-face[", indx, "] from json=", jsonData);
        }
        auto connectionData = jsonData["connections"];
        super(indx, connectionData);
    } // end constructor
    
} // end class PistonFace
