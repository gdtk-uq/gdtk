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
import ntypes.complex;
import nm.number;

import util.json_helper;
import geom;
import globalconfig;
import solidfvcell;

struct SolidProps {
public:
    double rho = 0.0;
    double k = 0.0;
    double Cp = 0.0;
    double k11 = 0.0;
    double k12 = 0.0;
    double k13 = 0.0;
    double k21 = 0.0;
    double k22 = 0.0;
    double k23 = 0.0;
    double k31 = 0.0;
    double k32 = 0.0;
    double k33 = 0.0;

    this(double rho_, double k_, double Cp_,
         double k11_=0.0, double k12_=0.0, double k13_=0.0,
         double k21_=0.0, double k22_=0.0, double k23_=0.0,
         double k31_=0.0, double k32_=0.0, double k33_=0.0)
    {
        rho = rho_;
        k = k_;
        Cp = Cp_;
        k11 = k11_; k12 = k12_; k13 = k13_;
        k21 = k21_; k22 = k22_; k23 = k23_;
        k31 = k31_; k32 = k32_; k33 = k33_;
    }

    // Create a new SolidProps object as a weighted average
    // of two other objects.
    this(SolidProps A, SolidProps B, double wA, double wB)
    {
        rho = wA*A.rho + wB*B.rho;
        k = wA*A.k + wB*B.k;
        Cp = wA*A.Cp + wB*B.Cp;
        k11 = wA*A.k11 + wB*B.k11;
        k12 = wA*A.k12 + wB*B.k12;
        k13 = wA*A.k13 + wB*B.k13;
        k21 = wA*A.k21 + wB*B.k21;
        k22 = wA*A.k22 + wB*B.k22;
        k23 = wA*A.k23 + wB*B.k23;
        k31 = wA*A.k31 + wB*B.k31;
        k32 = wA*A.k32 + wB*B.k32;
        k33 = wA*A.k33 + wB*B.k33;
    }

}

number updateEnergy(ref SolidProps sp, number T)
{
    return sp.rho*sp.Cp*T;
}

number updateTemperature(ref SolidProps sp, number e)
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

    number T = initTemperature;
    number e = updateEnergy(sp, T);

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
        formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
                       pos.x.re, pos.y.re, pos.z.re, volume.re, e.re, T.re,
                       sp.rho, sp.Cp, sp.k,
                       sp.k11, sp.k12, sp.k13,
                       sp.k21, sp.k22, sp.k23,
                       sp.k31, sp.k32, sp.k33);
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
