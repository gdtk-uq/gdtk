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

class SolidFVCell {
public:
    size_t id;
    // Cell properties
    double volume;
    double areaxy;
    Vector3 pos;
    // Cell state
    double[] T;
    double[] e;
    double[] dedt;
    // Cell source term
    double Q;
    // Connections
    SolidFVInterface[] iface;
    SolidFVVertex[] vtx;

    this()
    {
	T.length = n_time_levels;
	e.length = n_time_levels;
	dedt.length = n_time_levels;
    }

    void scanValuesFromString(string buffer)
    {
	auto items = split(buffer);
	pos.refx = to!double(items.front); items.popFront();
	pos.refy = to!double(items.front); items.popFront();
	pos.refz = to!double(items.front); items.popFront();
	volume = to!double(items.front); items.popFront();
	e[0] = to!double(items.front); items.popFront();
	T[0] = to!double(items.front); items.popFront();
    }

    string writeValuesToString() const
    {
	auto writer = appender!string();
	formattedWrite(writer, "%.16e %.16e %.16e %.16e %.16e %.16e",
		       pos.x, pos.y, pos.z, volume, e[0], T[0]);
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
	double gamma1 = 1.0;
	e[1] = e[0] + dt*gamma1*dedt[0];
    }
    void stage2Update(double dt)
    {
	// Assuming predictor-corrector
	double gamma1 = 0.5;
	double gamma2 = 0.5;
	e[2] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1]);
    }

}

string[] varListForSolidCell()
{
    string[] list;
    list ~= ["pos.x", "pos.y", "pos.z", "volume"];
    list ~= ["e", "T"];
    return list;
}
