/*
 * cea_therm_cond.d
 * Thermal conductivity computed using the CEA curve format.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-01-19
 */

module gas.diffusion.cea_therm_cond;

import gas.gas_model;
import gas.diffusion.therm_cond;
import std.math;
import util.lua;
import util.lua_service;
import std.c.stdlib : exit;
import std.stdio;
import std.string;

struct CEAThermCondCurve
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
	    throw new Exception("T_upper <= T_lower in CeaThermCondCurve:constructor.");
	A = params["A"];
	B = params["B"];
	C = params["C"];
	D = params["D"];
    }

    double eval(double T) const
    {
	if ( T < T_lower )
	    throw new Exception("temperature value lower than T_lower in CEAThermCondCurve:eval()");
	if ( T > T_upper )
	    throw new Exception("temperature value greater than T_upper in CEAThermCondCurve:eval()");
	double log_k = A*log(T) + B/T + C/(T*T) + D;

	/* CEA value is in microWatts/(cm.K), so convert to SI units of W/(m.K). */
	double k = exp(log_k)*1.0e-4;
	return k;
    }
}


class CEAThermalConductivity : ThermalConductivity
{
public:
    this(in CEAThermalConductivity src) {
	_curves = src._curves.dup;
	_T_lowest = src._T_lowest;
	_T_highest = src._T_highest;
    }
    this(in CEAThermCondCurve[] curves)
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
		throw new Exception("CEAThermalConductivity: curves are not continuous in temperature.");
	    }
	}
    }
    override CEAThermalConductivity dup() const
    {
	return new CEAThermalConductivity(this);
    }
    override double eval(in GasState Q, int imode) const
    {
	double T = Q.T[imode];
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
	throw new Exception("CEAThermalConductivity:eval() -- we should never reach this point.");
    }

private:
    CEAThermCondCurve[] _curves;
    double _T_lowest;
    double _T_highest;
}

CEAThermalConductivity createCEAThermalConductivity(lua_State* L)
{
    string[6] pList = ["T_lower", "T_upper", "A", "B", "C", "D"];
    double[string] params;
    auto nseg = getInt(L, -1, "nsegments");
    CEAThermCondCurve[] curves;
    foreach ( i; 0..nseg ) {
	auto key = format("segment%d", i);
	try {
	    getAssocArrayOfDoubles(L, key, pList, params);
	}
	catch ( Exception e ) {
	    writeln("ERROR: There was a problem reading in a value");
	    writeln("ERROR: when initialising a CEAThermalConductivity object");
	    writeln("ERROR: in function 'createCEAThermalConductivity'.");
	    writeln("ERROR: Exception message is:\n");
	    writeln(e.msg);
	    exit(1);
	}
	curves ~= CEAThermCondCurve(params);
    }
    return new CEAThermalConductivity(curves);
}

unittest
{
    import util.msg_service;
    /// First, let's test the CEAThermCondCurve on its own.
    double[string] params = ["T_lower":500.0, "T_upper":15000.0,
			     "A":0.76269502, "B":0.62341752e3,
			     "C":-0.71899552e6, "D":0.56927918];
    auto ceaCurve = CEAThermCondCurve(params);
    assert(approxEqual(0.1662583, ceaCurve.eval(7200.0)), failedUnitTest(__LINE__, __FILE__));

    /// Next, let's test the creation and functionality
    /// of a CEAThermalConductivity object.
    auto L = init_lua_State("sample-data/CO2-therm-cond.lua");
    lua_getglobal(L, "cea");
    auto co2CEA = createCEAThermalConductivity(L);
    lua_close(L);
    auto Q = new GasState(1, 1);
    Q.T[0] = 3500.0;
    assert(approxEqual(1.859070e-01, co2CEA.eval(Q, 0)), failedUnitTest(__LINE__, __FILE__));

}
