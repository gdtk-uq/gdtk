/**
 * cea_thermo_curves.d
 * Implements evaluation of Cp, h and s based
 * on the curve form used by the CEA program.
 *
 * Author: Rowan G. and Peter J.
 */

module gas.thermo.cea_thermo_curves;

import std.math;
import std.stdio;
import std.c.stdlib : exit;
import std.string;
import util.lua;
import util.lua_service;

struct CEAThermoCurve
{
public:
    double T_lower;
    double T_upper;
    double R;
    double[9] a;

    this(this)
    {
	a = a.dup;
    }
    double eval_Cp(double T) const
    {
	
	if ( T < T_lower ) 
	    throw new Exception("temperature value lower than T_lower in CEAThermoCurve.eval_Cp()");
	if ( T > T_upper )
	    throw new Exception("temperature value greater than T_upper in CEAThermoCurve.eval_Cp()");
	
	double Cp_on_R = a[0]/(T*T) + a[1]/T + a[2] + a[3]*T;
	Cp_on_R += a[4]*T*T + a[5]*T*T*T + a[6]*T*T*T*T;
	return R*Cp_on_R;
    }
    double eval_h(double T) const
    {
	if ( T < T_lower ) 
	    throw new Exception("temperature value lower than T_lower in CEAThermoCurve.eval_h()");
	if ( T > T_upper )
	    throw new Exception("temperature value greater than T_upper in CEAThermoCurve.eval_h()");
	double h_on_RT = -a[0]/(T*T) + a[1]/T * log(T) + a[2] + a[3]*T/2.0;
	h_on_RT +=  a[4]*T*T/3.0 + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]/T;
	return R*T*h_on_RT;
    }
    double eval_s(double T) const
    {
	if ( T < T_lower ) 
	    throw new Exception("temperature value lower than T_lower in CEAThermoCurve.eval_s()");
	if ( T > T_upper )
	    throw new Exception("temperature value greater than T_upper in CEAThermoCurve.eval_s()");
	double s_on_R = -a[0]/(2.0*T*T) - a[1]/T + a[2]*log(T) + a[3]*T;
	s_on_R += a[4]*T*T/2.0 + a[5]*T*T*T/3.0 + a[6]*T*T*T*T/4.0 + a[8];
	return R*s_on_R;
    }
}

class CEAThermo
{
public:
    this(in CEAThermo src) {
	_curves = src._curves.dup;
	_R = src._R;
	_T_lowest = src._T_lowest;
	_T_highest = src._T_highest;
    }
    this(in CEAThermoCurve[] curves)
    {
	_curves = curves.dup;
	/* We want the curves to ordered based on
	   temperature ranges from lowest to highest,
	   and the range of temperatures needs to be
	   continuous. This is the callers responsibility
	   to pass in the curves this way. If they aren't
	   correct, we won't contruct the CEAThermo object.
	*/
	_T_lowest = _curves[0].T_lower;
	_T_highest = _curves[$-1].T_upper;
	for ( auto i=1; i < _curves.length; ++i ) {
	    if ( _curves[i].T_lower != _curves[i-1].T_upper )
		throw new Exception("CEAThermo: curves are not continuous in temperature.");
	}
    }
    CEAThermo dup() const
    {
	return new CEAThermo(this);
    }
    double eval_Cp(double T) const
    {
	/* We don't try to any fancy extrapolation
	   off the ends of the curves. Beyond the range
	   of validity we simply take the value at
	   the end of the range.
	*/
	if ( T < _T_lowest ) {
	    return _curves[0].eval_Cp(_T_lowest);
	}
	if ( T > _T_highest ) {
	    return _curves[$-1].eval_Cp(_T_highest);
	}
	// Search for appropriate curve segment and evaluate.
	foreach ( c; _curves ) {
	    if ( T >= c.T_lower && T <= c.T_upper )
		return c.eval_Cp(T);
	}
	// We should never reach this point.
	throw new Exception("CEAThermo.eval_Cp(): we should never have reached this point.");
    }
    double eval_h(double T) const
    {
	/* We don't try any fancy extrapolation beyond the limits
	   of the curves. We will simply take Cp as constant and
	   integrate using that assumption to get the portion of
	   enthalpy that is contributed beyond the temperature range.
	*/
	if ( T < _T_lowest ) {
	    double h_lowest = _curves[0].eval_h(_T_lowest);
	    double Cp_lowest = _curves[0].eval_Cp(_T_lowest);
	    double h = h_lowest - Cp_lowest*(_T_lowest - T);
	    return h;
	}
	if ( T > _T_highest ) {
	    double h_highest = _curves[$-1].eval_h(_T_highest);
	    double Cp_highest = _curves[$-1].eval_Cp(_T_highest);
	    double h = h_highest + Cp_highest*(T - _T_highest);
	    return h;
	}
	// Search for appropriate curve segment and evaluate.
	foreach ( ref c; _curves ) {
	    if ( T >= c.T_lower && T <= c.T_upper )
		return c.eval_h(T);
	}
	assert(0);
	/* 
	// We should never reach this point.
	throw new Exception("CEAThermo.eval_h(): we should never have reached this point.");
	*/
    }
    double eval_s(double T) const
    {
	/* We don't try any fancy extrapolation beyond the limits
	   of the curves. We will simply take Cp as constant and
	   integrate using that assumption to get the portion of
	   entropy that is contributed beyond the temperature range.
	*/
	if ( T < _T_lowest ) {
	    double s_lowest = _curves[0].eval_s(_T_lowest);
	    double Cp_lowest = _curves[0].eval_Cp(_T_lowest);
	    double s = s_lowest - Cp_lowest*log(_T_lowest/T);
	    return s;
	}
	if ( T > _T_highest ) {
	    double s_highest = _curves[$-1].eval_s(_T_highest);
	    double Cp_highest = _curves[$-1].eval_Cp(_T_highest);
	    double s = s_highest + Cp_highest*(T/_T_highest);
	    return s;
	}
	// Search for appropriate curve segment and evaluate.
	foreach ( c; _curves ) {
	    if ( T >= c.T_lower && T <= c.T_upper )
		return c.eval_s(T);
	}
	// We should never reach this point.
	throw new Exception("CEAThermo.eval_s(): we should never have reached this point.");
    }
private:
    double _R;
    double _T_lowest;
    double _T_highest;
    CEAThermoCurve[] _curves;
}

CEAThermo createCEAThermo(lua_State* L, double R)
{
    auto nseg = getInt(L, -1, "nsegments");
    CEAThermoCurve[] curves;
    foreach ( i; 0..nseg ) {
	auto key = format("segment%d", i);
	lua_getfield(L, -1, key.toStringz);
	auto c = CEAThermoCurve();
	try {
	    c.R = R;
	    c.T_lower = getDouble(L, -1, "T_lower");
	    c.T_upper = getDouble(L, -1, "T_upper");
	    double[] coeffs;
	    getArrayOfDoubles(L, -1, "coeffs", coeffs);
	    foreach ( ic; 0..c.a.length ) {
		c.a[ic] = coeffs[ic];
	    }
	}
	catch ( Exception e ){
	    writeln("ERROR: There was a problem reading in one of the");
	    writeln("ERROR: thermo curves while trying to construct");
	    writeln("ERROR: a CEAThermo object in function 'createCEAThermo'.");
	    writeln("ERROR: Exception message is:\n");
	    writeln(e.msg);
	    exit(1);
	}
	curves ~= c;
	lua_pop(L, 1);
    }
    return new CEAThermo(curves);
}

unittest
{
    import util.msg_service;
    // 1. Test a segment (CEAThermoCurve) on its own.
    auto curve = CEAThermoCurve();
    curve.T_lower = 6000.000;
    curve.T_upper = 20000.000;
    curve.R =  8.31451/14.0067000e-3;
    curve.a[0] = 5.475181050e+08; curve.a[1] = -3.107574980e+05; curve.a[2] = 6.916782740e+01;
    curve.a[3] = -6.847988130e-03; curve.a[4] = 3.827572400e-07; curve.a[5] = -1.098367709e-11;
    curve.a[6] = 1.277986024e-16; curve.a[7] = 2.550585618e+06; curve.a[8] = -5.848769753e+02;
    assert(approxEqual(2022.9958, curve.eval_Cp(7500.0)), failedUnitTest(__LINE__, __FILE__));
    // 2. Test full curve
    auto L = init_lua_State("sample-data/O-thermo.lua");
    lua_getglobal(L, "CEA_coeffs");
    double R = 8.31451/0.0159994;
    auto oThermo = createCEAThermo(L, R);
    lua_close(L);
    assert(approxEqual(1328.627, oThermo.eval_Cp(500.0)), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(20030794.683, oThermo.eval_h(3700.0)), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(14772.717, oThermo.eval_s(10000.0)), failedUnitTest(__LINE__, __FILE__));
}
