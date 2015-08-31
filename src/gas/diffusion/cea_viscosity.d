/*
 * cea_viscosity.d
 * Viscosity computed using the CEA curves.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-01-13
 */

module gas.diffusion.cea_viscosity;

import gas.gas_model;
import gas.diffusion.viscosity;
import std.math;
import util.lua;
import util.lua_service;
import std.c.stdlib : exit;
import std.stdio;
import std.string;

struct CEAViscCurve
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
	    throw new Exception("T_upper <= T_lower in CeaViscCurve:constructor.");
	A = params["A"];
	B = params["B"];
	C = params["C"];
	D = params["D"];
    }

    double eval(double T) const
    {
	if ( T < T_lower )
	    throw new Exception("temperature value lower than T_lower in CeaViscCurve:eval()");
	if ( T > T_upper )
	    throw new Exception("temperature value greater than T_upper in CeaViscCurve:eval()");
	double log_mu = A*log(T) + B/T + C/(T*T) + D;

	/* CEA value is in microPoise, so convert to SI units.
	   1 P = 0.1 kg/(m.s)
	   1 microP = 1.0e-6 * 0.1 = 1.0e-7 kg/(m.s)
	*/
	double mu = exp(log_mu)*1.0e-7;
	return mu;
    }
}


class CEAViscosity : Viscosity
{
public:
    this(in CEAViscosity src) {
	_curves = src._curves.dup;
	_T_lowest = src._T_lowest;
	_T_highest = src._T_highest;
    }
    this(in CEAViscCurve[] curves)
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
		throw new Exception("CEAViscosity: curves are not continuous in temperature.");
	    }
	}
    }
    override CEAViscosity dup() const
    {
	return new CEAViscosity(this);
    }
    override double eval(in GasState Q) const
    {
	double T = Q.T[0];
	// At the limits of the curve, extrapolate value as a constant.
	if ( T < _T_lowest ) {
	    return _curves[0].eval(_T_lowest);
	}
	if ( T > _T_highest ) {
	    return _curves[$-1].eval(_T_highest);
	}
	// Search for curve segment and evaluate
	foreach ( c; _curves ) {
	    if ( T >= c.T_lower && T <= c.T_upper ) {
		return c.eval(T);
	    }
	}
	// We should not reach this point.
	throw new Exception("CEAViscosity:eval() -- we should never reach this point.");
    }

private:
    CEAViscCurve[] _curves;
    double _T_lowest;
    double _T_highest;
}

CEAViscosity createCEAViscosity(lua_State* L)
{
    string[6] pList = ["T_lower", "T_upper", "A", "B", "C", "D"];
    double[string] params;
    auto nseg = getInt(L, -1, "nsegments");
    CEAViscCurve[] curves;
    foreach ( i; 0..nseg ) {
	auto key = format("segment%d", i);
	try {
	    getAssocArrayOfDoubles(L, key, pList, params);
	}
	catch ( Exception e ) {
	    writeln("ERROR: There was a problem reading in a value");
	    writeln("ERROR: when initialising a CEAViscosity object");
	    writeln("ERROR: in function 'createCEAViscosity'.");
	    writeln("ERROR: Exception message is:\n");
	    writeln(e.msg);
	    exit(1);
	}
	curves ~= CEAViscCurve(params);
    }
    return new CEAViscosity(curves);
}

unittest
{
    import util.msg_service;
    /// First, let's test the CEAViscCurve on its own.
    double[string] params = ["T_lower":200.0, "T_upper":1000.0,
			     "A":0.62526577, "B":-0.31779652e2,
			     "C":-0.1640798e4, "D":0.17454992e01];
    auto ceaCurve = CEAViscCurve(params);
    assert(approxEqual(3.8818e-5, ceaCurve.eval(900.0)), failedUnitTest());

    /// Next, let's test the creation and functionality
    /// of a CEAViscosity object.
    auto L = init_lua_State("sample-data/O2-viscosity.lua");
    lua_getglobal(L, "cea");
    auto o2CEA = createCEAViscosity(L);
    lua_close(L);
    auto Q = new GasState(1, 1);
    Q.T[0] = 1500.0;
    assert(approxEqual(6.407851e-05, o2CEA.eval(Q)), failedUnitTest());
}
