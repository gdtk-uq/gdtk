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
import util.lua;
import util.lua_service : getDouble, getDoubleWithDefault;

import json_helper;
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

    // Create a SolidProps object from a Lua table
    this(lua_State* L, int tblIdx)
    {
        rho = getDouble(L, tblIdx, "rho");
        k = getDouble(L, tblIdx, "k");
        Cp = getDouble(L, tblIdx, "Cp");
        k11 = getDoubleWithDefault(L, tblIdx, "k11", 0.0);
        k12 = getDoubleWithDefault(L, tblIdx, "k12", 0.0);
        k13 = getDoubleWithDefault(L, tblIdx, "k13", 0.0);
        k21 = getDoubleWithDefault(L, tblIdx, "k21", 0.0);
        k22 = getDoubleWithDefault(L, tblIdx, "k22", 0.0);
        k23 = getDoubleWithDefault(L, tblIdx, "k23", 0.0);
        k31 = getDoubleWithDefault(L, tblIdx, "k31", 0.0);
        k32 = getDoubleWithDefault(L, tblIdx, "k32", 0.0);
        k33 = getDoubleWithDefault(L, tblIdx, "k33", 0.0);
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
