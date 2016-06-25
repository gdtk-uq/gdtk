/**
 * solidprops.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module solidprops;

import std.stdio;
import std.array;
import std.format;
import std.json;

import json_helper;
import geom;
import sgrid;
import globalconfig;
import solidfvcell;

class SolidProps {
public:
    double rho;
    double k;
    double Cp;
    double k11;
    double k12;
    double k22;

    this(double rho_, double k_, double Cp_,
	 double k11_=0.0, double k12_=0.0, double k22_=0.0)
    {
	rho = rho_;
	k = k_;
	Cp = Cp_;
	k11 = k11_;
	k12 = k12_;
	k22 = k22_;
    }
}

SolidProps makeSolidPropsFromJson(JSONValue jsonData)
{
    auto rho = getJSONdouble(jsonData, "rho", 8960.0);
    auto k = getJSONdouble(jsonData, "k", 401.0);
    auto Cp = getJSONdouble(jsonData, "Cp", 386.0);
    return new SolidProps(rho, k, Cp);
}

double updateEnergy(SolidProps sp, double T)
{
    return sp.rho*sp.Cp*T;
}

double updateTemperature(SolidProps sp, double e)
{
    return e/(sp.rho*sp.Cp);
}

void writeInitialSolidFile(string fileName, ref StructuredGrid grid,
			   double initTemperature, ref SolidProps sp, double t0)
{
    // Numbers of cells derived from numbers of vertices.
    auto nic = grid.niv - 1;
    auto njc = grid.njv - 1;
    auto nkc = grid.nkv - 1;
    if (GlobalConfig.dimensions == 2) nkc = 1;

    double T = initTemperature;
    double e = updateEnergy(sp, T);

    string cellDataToString(size_t i, size_t j, size_t k)
    {
	Vector3 p000 = *grid[i,j,k];
	Vector3 p100 = *grid[i+1,j,k];
	Vector3 p110 = *grid[i+1,j+1,k];
	Vector3 p010 = *grid[i,j+1,k];
	// [TODO] provide better calculation using geom module.
	// For the moment, it doesn't matter greatly because the solver 
	// will compute it's own approximations
	auto pos = 0.25*(p000 + p100 + p110 + p010);
	auto volume = 0.0;
	if (GlobalConfig.dimensions == 3) {
	    Vector3 p001 = *grid[i,j,k+1];
	    Vector3 p101 = *grid[i+1,j,k+1];
	    Vector3 p111 = *grid[i+1,j+1,k+1];
	    Vector3 p011 = *grid[i,j+1,k+1];
	    pos = 0.5*pos + 0.125*(p001 + p101 + p111 + p011);
	}
	// Should match SolidFVCell.writeValuesToString()
	auto writer = appender!string();
	formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e",
		       pos.x, pos.y, pos.z, volume, e, T);
	return writer.data;
    }

    // Write the data for the whole structured block.
    auto f = File(fileName, "w");
    f.writefln("%.18e", t0);
    // Variable list for cell on one line.
    auto writer = appender!string();
    foreach(varname; varListForSolidCell()) {
	formattedWrite(writer, " \"%s\"", varname);
    }
    f.writeln(writer.data);
    // Numbers of cells.
    f.writefln("%d %d %d", nic, njc, nkc);
    // The actual cell data.
    foreach (k; 0 .. nkc) {
	foreach (j; 0 .. njc) {
	    foreach (i; 0 .. nic) {
		f.writefln(" %s", cellDataToString(i,j,k));
	    }
	}
    }
    f.close();
    return;
} // end write_initial_flow_file()
