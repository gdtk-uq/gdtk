// endcondition.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module endcondition;


import std.conv;
import std.stdio;
import std.file;
import std.string;
import std.json;
import std.format;
import std.algorithm;
import std.math;

import util.json_helper;
import geom;
import gas;
import gasdyn.gasflow;
import config;
import tube;
import gasslug;
import lcell;
import piston;
import simcore;
import misc;

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


enum DiaphragmState {closed=0, triggered=1, open=2};

class Diaphragm : EndCondition {
public:
    double p_burst;
    double dt_hold;
    double dxL;
    double dxR;
    DiaphragmState state;
    double t_open;

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

    void read_data(File fp, int tindx)
    {
        string text = fp.readln().chomp();
        while (text.canFind("#")) { text = fp.readln().chomp(); }
        string[] items = text.split();
        int myTindx = to!int(items[0]);
        while (myTindx < tindx) {
            text = fp.readln().chomp();
            items = text.split();
            myTindx = to!int(items[0]);
        }
        // We should be at the line that contains the requested tindx.
        state = to!DiaphragmState(to!int(items[1]));
    } // end read_data()

    void write_data(File fp, int tindx, bool write_header)
    {
        if (write_header) { fp.writeln("# tindx  is_burst"); }
        fp.writeln(format("%d %d", tindx, to!int(state)));
    } // end write_data()

    void update_state(double t)
    {
        final switch (state) {
        case DiaphragmState.closed:
            // When sampling the pressure either side of the diaphragm,
            // we average the cell pressures over distances dxL, dxR.
            double pL = 0.0;
            int ncellL = 0;
            if (slugL) {
                if (slugL_end == End.L) {
                    double x0 = slugL.faces[0].x;
                    foreach (i; 0 .. slugL.ncells) {
                        LCell cL = slugL.cells[i];
                        pL += cL.gas.p;
                        ncellL++;
                        // Doing the test after guarantees at least one cell is sampled,
                        // even for dxL zero.
                        if (fabs(cL.xmid - x0) > dxL) { break; }
                    }
                } else {
                    double x0 = slugL.faces[$-1].x;
                    foreach (i; 0 .. slugL.ncells) {
                        LCell cL = slugL.cells[$-1-i];
                        pL += cL.gas.p;
                        ncellL++;
                        if (fabs(cL.xmid - x0) > dxL) { break; }
                    }
                }
            }
            pL /= ncellL;
            double pR = 0.0;
            int ncellR = 0;
            if (slugR) {
                if (slugR_end == End.L) {
                    double x0 = slugR.faces[0].x;
                    foreach (i; 0 .. slugR.ncells) {
                        LCell cR = slugR.cells[0];
                        pR += cR.gas.p;
                        ncellR++;
                        if (fabs(cR.xmid - x0) > dxR) { break; }
                    }
                } else {
                    double x0 = slugR.faces[$-1].x;
                    foreach (i; 0 .. slugR.ncells) {
                        LCell cR = slugR.cells[$-1-i];
                        pR += cR.gas.p;
                        ncellR++;
                        if (fabs(cR.xmid - x0) > dxR) { break; }
                    }
                }
            }
            pR /= ncellR;
            if (fabs(pL-pR) > p_burst) {
                t_open = t + dt_hold;
                state = DiaphragmState.triggered;
                string msg = format("t=%e Diaphragm at ec_indx=%d triggered\n", t, indx);
                write(msg);
                append(L1dConfig.job_name~"/events.txt", msg);
            }
            break;
        case DiaphragmState.triggered:
            if (t > t_open) {
                state = DiaphragmState.open;
                string msg = format("t=%e Diaphragm at ec_indx=%d opened\n", t, indx);
                write(msg);
                append(L1dConfig.job_name~"/events.txt", msg);
            }
            break;
        case DiaphragmState.open:
            // do nothing
            break;
        }
    } // end update_state()
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
