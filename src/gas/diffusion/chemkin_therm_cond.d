/*
 * chemkin_therm_cond.d
 * Thermal conductivity computed using the Chemkin curve format.
 *
 * Author: Oliver S., adapted from cea_therm_cond.d by Rowan G. and Peter J.
 * Version: 2019-04-14
 */

module gas.diffusion.chemkin_therm_cond;

import std.math;
import std.stdio;
import std.string;
import std.conv;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas.gas_model;
import gas.gas_state;
import gas.diffusion.therm_cond;

struct ChemkinThermCondCurve
{
public:
    double T_lower;
    double T_upper;
    double A;
    double B;
    double C;
    double D;

    this(double[string] params)
    {
        T_lower = params["T_lower"];
        T_upper = params["T_upper"];
        if ( T_upper <= T_lower )
            throw new Exception("T_upper <= T_lower in ChemkinThermCondCurve:constructor.");
        A = params["A"];
        B = params["B"];
        C = params["C"];
        D = params["D"];
    }

    @nogc number eval(number T) const
    {
        number logT = log(T);
        return eval(T, logT);
    }

    @nogc number eval(number T, number logT) const
    {
        if ( T < T_lower )
            throw new Exception("temperature value lower than T_lower in ChemkinThermCondCurve:eval()");
        if ( T > T_upper )
            throw new Exception("temperature value greater than T_upper in ChemkinThermCondCurve:eval()");
        number log_k = A + B*logT + C*logT*logT + D*logT*logT*logT;

        /* Chemkin value is in W/(m.K), SI units. */
        number k = exp(log_k);
        return k;
    }
}


class ChemkinThermalConductivity : ThermalConductivity
{
public:
    this(in ChemkinThermalConductivity src) {
        _curves = src._curves.dup;
        _T_lowest = src._T_lowest;
        _T_highest = src._T_highest;
    }
    this(in ChemkinThermCondCurve[] curves)
    {
        _curves = curves.dup;
        /* We want the curves to be ordered based on their
           temperature ranges, and continuous in their
           temperature ranges. This is the caller's responsibility.
           If they aren't, we won't construct an object.
        */
        _T_lowest = _curves[0].T_lower;
        _T_highest = _curves[$-1].T_upper;
        for ( auto i=1; i < _curves.length; ++i ) {
            if ( _curves[i].T_lower != _curves[i-1].T_upper ) {
                throw new Exception("ChemkinThermalConductivity: curves are not continuous in temperature.");
            }
        }
    }
    override ChemkinThermalConductivity dup() const
    {
        return new ChemkinThermalConductivity(this);
    }
    override number eval(number T) const
    {
        number logT = log(T);
        return eval(T, logT);
    }
    override number eval(number T, number logT) const
    {
        // At the limits of the curve, extrapolate value as a constant.
        if ( T < _T_lowest ) {
            return _curves[0].eval(to!number(_T_lowest));
        }
        if ( T > _T_highest ) {
            return _curves[$-1].eval(to!number(_T_highest));
        }
        // Search for curve segment and evaluate
        foreach ( c; _curves ) {
            if ( T >= c.T_lower && T <= c.T_upper ) {
                return c.eval(T, logT);
            }
        }
        // We should not reach this point.
        throw new Exception("ChemkinThermalConductivity:eval() -- we should never reach this point.");
    }

private:
    ChemkinThermCondCurve[] _curves;
    double _T_lowest;
    double _T_highest;
}

ChemkinThermalConductivity createChemkinThermalConductivity(lua_State* L)
{
    string[6] pList = ["T_lower", "T_upper", "A", "B", "C", "D"];
    double[string] params;
    auto nseg = getInt(L, -1, "nsegments");
    ChemkinThermCondCurve[] curves;
    foreach ( i; 0..nseg ) {
        auto key = format("segment%d", i);
        try {
            getAssocArrayOfDoubles(L, key, pList, params);
        }
        catch ( Exception e ) {
            writeln("ERROR: There was a problem reading in a value");
            writeln("ERROR: when initialising a ChemkinThermalConductivity object");
            writeln("ERROR: in function 'createChemkinThermalConductivity'.");
            writeln("ERROR: Exception message is:\n");
            writeln(e.msg);
            throw new Exception(e.msg);
        }
        curves ~= ChemkinThermCondCurve(params);
    }
    return new ChemkinThermalConductivity(curves);
}

version(chemkin_therm_cond_test) {
    import util.msg_service;
    int main() {

        /// First, let's test the ChemkinThermCondCurve on its own.
        double[string] params = ["T_lower":200.0, "T_upper":3500.0,
                                 "A":-2.119099884033e+01, "B":5.186927559697e+00,
                                 "C":-4.741229077145e-01, "D":1.610702319175e-02];
        auto chemkinCurve = ChemkinThermCondCurve(params);
        assert(isClose(0.124382229603, chemkinCurve.eval(2000.0), 1.0e-6), failedUnitTest());

        /*
        /// Next, let's test the creation and functionality
        /// of a ChemkinThermalConductivity object.
        auto L = init_lua_State();
        doLuaFile(L, "sample-data/CO2-therm-cond.lua");
        lua_getglobal(L, "chemkin");
        auto co2Chemkin = createChemkinThermalConductivity(L);
        lua_close(L);
        auto Q = GasState(1, 1);
        // Q.T = 3500.0;
        assert(isClose(0.185719655303, co2Chemkin.eval(Q, 3500.0), 1.0e-6), failedUnitTest());
        */

        return 0;
    }
}
