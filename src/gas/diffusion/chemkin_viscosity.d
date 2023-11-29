/*
 * chemkin_viscosity.d
 * Viscosity computed using the Chemkin curves.
 *
 * Author: Oliver S., adapted from cea_viscosity.d by Rowan G. and Peter J.
 * Version: 2019-04-14
 */

module gas.diffusion.chemkin_viscosity;

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
import gas.diffusion.viscosity;

struct ChemkinViscCurve
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
            throw new Exception("T_upper <= T_lower in ChemkinViscCurve:constructor.");
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
            throw new Exception("temperature value lower than T_lower in ChemkinViscCurve:eval()");
        if ( T > T_upper )
            throw new Exception("temperature value greater than T_upper in ChemkinViscCurve:eval()");
        number log_mu = A + B*logT + C*logT*logT + D*logT*logT*logT;

        /* Chemkin value is in SI units.
        */
        number mu = exp(log_mu);
        return mu;
    }
}


class ChemkinViscosity : Viscosity
{
public:
    this(in ChemkinViscosity src) {
        _curves = src._curves.dup;
        _T_lowest = src._T_lowest;
        _T_highest = src._T_highest;
    }
    this(in ChemkinViscCurve[] curves)
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
                throw new Exception("ChemkinViscosity: curves are not continuous in temperature.");
            }
        }
    }
    override ChemkinViscosity dup() const
    {
        return new ChemkinViscosity(this);
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
                return c.eval(T);
            }
        }
        // We should not reach this point.
        throw new Exception("ChemkinViscosity:eval() -- we should never reach this point.");
    }
    override number eval(number T) const
    {
        number logT = log(T);
        return eval(T, logT);
    }

private:
    ChemkinViscCurve[] _curves;
    double _T_lowest;
    double _T_highest;
}

ChemkinViscosity createChemkinViscosity(lua_State* L)
{
    string[6] pList = ["T_lower", "T_upper", "A", "B", "C", "D"];
    double[string] params;
    auto nseg = getInt(L, -1, "nsegments");
    ChemkinViscCurve[] curves;
    foreach ( i; 0..nseg ) {
        auto key = format("segment%d", i);
        try {
            getAssocArrayOfDoubles(L, key, pList, params);
        }
        catch ( Exception e ) {
            writeln("ERROR: There was a problem reading in a value");
            writeln("ERROR: when initialising a ChemkinViscosity object");
            writeln("ERROR: in function 'createChemkinViscosity'.");
            writeln("ERROR: Exception message is:\n");
            writeln(e.msg);
            throw new Exception(e.msg);
        }
        curves ~= ChemkinViscCurve(params);
    }
    return new ChemkinViscosity(curves);
}

version(chemkin_viscosity_test) {

    import util.msg_service;
    int main() {
        /// First, let's test the ChemkinViscCurve on its own.
        double[string] params = ["T_lower":200.0, "T_upper":3500.0,
                                 "A":-2.000056657134e+01, "B":2.898260367046e+00,
                                 "C":-3.018511062902e-01, "D":1.346540184145e-02];
        auto chemkinCurve = ChemkinViscCurve(params);
        assert(isClose(4.47483851e-05, chemkinCurve.eval(900.0), 1.0e-6), failedUnitTest());

        /*
        /// Next, let's test the creation and functionality
        /// of a ChemkinViscosity object.
        auto L = init_lua_State();
        doLuaFile(L, "sample-data/O2-viscosity.lua");
        lua_getglobal(L, "chemkin");
        auto o2Chemkin = createChemkinViscosity(L);
        lua_close(L);
        auto Q = GasState(1, 1);
        Q.T = 1500.0;
        assert(isClose(6.238853e-05, o2chemkin.eval(Q), 1.0e-6), failedUnitTest());
        */

        return 0;
    }
}
