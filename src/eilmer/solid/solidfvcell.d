/**
 * solidfvcell.d
 *
 * A solid finite-volume cell, to be held by SolidBlock objects.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module solidfvcell;

import std.conv;
import std.string;
import std.array;
import std.format;
import geom;
import fvcore;
import solidfvinterface;
import solidfvvertex;
import solidprops;
import std.stdio;
import globalconfig; //Anand added this to access GlobalConfig

class SolidFVCell {
public:
    size_t id;
    // Cell properties
    double volume;
    double areaxy;
    Vector3 pos;
    // Cell state
    double T;
    double[] e;
    double[] dedt;
    // Cell source term
    double Q;
    // Connections
    SolidFVInterface[] iface;
    SolidFVVertex[] vtx;

    this(size_t n_flow_time_levels)
    {
	e.length = n_flow_time_levels;
	dedt.length = n_flow_time_levels;
    }

    void scanValuesFromString(string buffer)
    {
	auto items = split(buffer);
	pos.refx = to!double(items.front); items.popFront();
	pos.refy = to!double(items.front); items.popFront();
	pos.refz = to!double(items.front); items.popFront();
	volume = to!double(items.front); items.popFront();
	e[0] = to!double(items.front); items.popFront();
	T = to!double(items.front); items.popFront();
    }

    string writeValuesToString() const
    {
	auto writer = appender!string();
	formattedWrite(writer, "%.16e %.16e %.16e %.16e %.16e %.16e",
		       pos.x, pos.y, pos.z, volume, e[0], T);
	return writer.data;
    }


    void timeDerivatives(int ftl, int dimensions)
    {
	SolidFVInterface IFn = iface[Face.north];
	SolidFVInterface IFe = iface[Face.east];
	SolidFVInterface IFs = iface[Face.south];
	SolidFVInterface IFw = iface[Face.west];
	SolidFVInterface IFt, IFb;
	if (dimensions == 3) {
	    IFt = iface[Face.top];
	    IFb = iface[Face.bottom];
	}
	// Cell volume (inverted).
	double volInv = 1.0 / volume;
	double integral;
	
	// Sum up fluxes (of form q.n)
	integral = -IFe.flux * IFe.area - IFn.flux * IFn.area
	    + IFw.flux * IFw.area + IFs.flux * IFs.area;
	dedt[ftl] = volInv * integral + Q;
    }
    void stage1Update(double dt)
    {
	double gamma1 = 1.0; // Assume Euler
//	if (!force_euler) {
	    final switch (GlobalConfig.gasdynamic_update_scheme) {
	    case GasdynamicUpdate.euler:
	    case GasdynamicUpdate.moving_grid_1_stage:
	    case GasdynamicUpdate.moving_grid_2_stage:
	    case GasdynamicUpdate.pc: gamma1 = 1.0; break;
	    case GasdynamicUpdate.midpoint: gamma1 = 0.5; break;
	    case GasdynamicUpdate.classic_rk3: gamma1 = 0.5; break;
	    case GasdynamicUpdate.tvd_rk3: gamma1 = 1.0; break;
	    case GasdynamicUpdate.denman_rk3: gamma1 = 8.0/15.0; break;
	   }	
//    }
   	e[1] = e[0] + dt*gamma1*dedt[0];
   }
    void stage2Update(double dt)
    {
	// Assuming predictor-corrector
	double gamma1 = 0.5;
	double gamma2 = 0.5;
	final switch (GlobalConfig.gasdynamic_update_scheme) {
	case GasdynamicUpdate.euler:
	case GasdynamicUpdate.moving_grid_1_stage: assert(false, "invalid for 1-stage update.");
	case GasdynamicUpdate.moving_grid_2_stage:
	case GasdynamicUpdate.pc: gamma1 = 0.5, gamma2 = 0.5; break;
	case GasdynamicUpdate.midpoint: gamma1 = 0.0; gamma2 = 1.0; break;
	case GasdynamicUpdate.classic_rk3: gamma1 = -1.0; gamma2 = 2.0; break;
	case GasdynamicUpdate.tvd_rk3: gamma1 = 0.25; gamma2 = 0.25; break;
	case GasdynamicUpdate.denman_rk3: gamma1 = -17.0/60.0; gamma2 = 5.0/12.0; break;
	}
	e[2] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1]);
    }
	void stage3Update(double dt)
    {
	// Assuming TVD_RK3 scheme as done in flow update
	double gamma1 = 1.0/6.0; // presume TVD_RK3 scheme.
	double gamma2 = 1.0/6.0;
	double gamma3 = 4.0/6.0;
	final switch (GlobalConfig.gasdynamic_update_scheme) {
	case GasdynamicUpdate.euler:
	case GasdynamicUpdate.moving_grid_1_stage:
	case GasdynamicUpdate.moving_grid_2_stage:
	case GasdynamicUpdate.pc:
	case GasdynamicUpdate.midpoint:
	    assert(false, "invalid for 2-stage update.");
	case GasdynamicUpdate.classic_rk3: gamma1 = 1.0/6.0; gamma2 = 4.0/6.0; gamma3 = 1.0/6.0; break;
	case GasdynamicUpdate.tvd_rk3: gamma1 = 1.0/6.0; gamma2 = 1.0/6.0; gamma3 = 4.0/6.0; break;
	    // FIX-ME: Check that we have Andrew Denman's scheme ported correctly.
	case GasdynamicUpdate.denman_rk3: gamma1 = 0.0; gamma2 = -5.0/12.0; gamma3 = 3.0/4.0; break;
	}
	e[3] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1] + gamma3*dedt[2]);
    }
}

string[] varListForSolidCell()
{
    string[] list;
    list ~= ["pos.x", "pos.y", "pos.z", "volume"];
    list ~= ["e", "T"];
    return list;
}
