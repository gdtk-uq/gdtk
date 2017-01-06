/**
 * rate_constant.d
 * Implements the calculation of reaction rate coefficients.
 *
 * Author: Rowan G.
 */

module kinetics.rate_constant;

import std.math;
import std.algorithm;
import std.typecons;
import std.string;
import util.lua;
import util.lua_service;
import gas;

/++
  RateConstant is a class for computing the rate constant of
  a chemical reaction.
+/

interface RateConstant {
    RateConstant dup();
    double eval(in GasState Q);
}

/++
 ArrheniusRateConstant uses the Arrhenius model to compute
 a chemical rate constant.

 Strictly speaking, this class implements a modified or
 generalised version of the Arrhenius rate constant.
 Arrhenius expression contains no temperature-dependent
 pre-exponential term. That is, the Arrhenius expression for
 rate constant is:
   k = A exp(-C/T))
 The modified expression is:
   k = A T^n exp(-C/T))

 This class implements that latter expression. It seems overly
 pedantic to provide two separate classes: one for the original
 Arrhenius expression and one for the modified version. Particularly
 given that almost all modern reaction schemes for gas phase chemistry
 provide rate constant inputs in terms of the modified expression.
 Simply setting n=0 gives the original Arrhenius form of the expression.
+/
class ArrheniusRateConstant : RateConstant {
public:
    this(double A, double n, double C)
    {
	_A = A;
	_n = n;
	_C = C;
    }
    this(lua_State* L)
    {
	_A = getDouble(L, -1, "A");
	_n = getDouble(L, -1, "n");
	_C = getDouble(L, -1, "C");
    }
    ArrheniusRateConstant dup()
    {
	return new ArrheniusRateConstant(_A, _n, _C);
    }
    override double eval(in GasState Q)
    {
	double T = Q.Ttr;
	return _A*pow(T, _n)*exp(-_C/T);
    }
private:
    double _A, _n, _C;
}


double thirdBodyConcentration(in GasState Q, Tuple!(int, double)[] efficiencies, GasModel gmodel)
{
    double val = 0.0;
    foreach (e; efficiencies) {
	int isp = e[0];
	double eff = e[1];
	val += eff * (Q.massf[isp]*Q.rho/gmodel.mol_masses[isp]);
    }
    return val;
}

class LHRateConstant : RateConstant {
public:
    this(ArrheniusRateConstant kInf, ArrheniusRateConstant k0,
	 Tuple!(int, double)[] efficiencies, GasModel gmodel)
    {
	_kInf = kInf.dup();
	_k0 = k0.dup();
	_efficiencies = efficiencies.dup();
	_gmodel = gmodel;
    }
    this(lua_State* L, Tuple!(int, double)[] efficiencies, GasModel gmodel)
    {
	lua_getfield(L, -1, "kInf");
	_kInf = new ArrheniusRateConstant(L);
	lua_pop(L, 1);

	lua_getfield(L, -1, "k0");
	_k0 = new ArrheniusRateConstant(L);
	lua_pop(L, 1);
	
	_efficiencies = efficiencies.dup();
	_gmodel = gmodel;
    }
    LHRateConstant dup()
    {
	return new LHRateConstant(_kInf, _k0, _efficiencies, _gmodel);
    }
    override double eval(in GasState Q)
    {
	double M = thirdBodyConcentration(Q, _efficiencies, _gmodel);
	double kInf = _kInf.eval(Q);
	double k0 = _k0.eval(Q);
	return k0*kInf*M/(k0*M + kInf);
    }
private:
    ArrheniusRateConstant _kInf, _k0;
    Tuple!(int, double)[] _efficiencies;
    GasModel _gmodel;
}

class TroeRateConstant : RateConstant {
public:
    this(ArrheniusRateConstant kInf, ArrheniusRateConstant k0, double Fcent,
	 Tuple!(int, double)[] efficiencies, GasModel gmodel)
    {
	_kInf = kInf.dup();
	_k0 = k0.dup();
	_Fcent = Fcent;
	_Fcent_supplied = true;
	_efficiencies = efficiencies.dup();
	_gmodel = gmodel;
    }
    this(lua_State* L, Tuple!(int, double)[] efficiencies, GasModel gmodel)
    {
	lua_getfield(L, -1, "kInf");
	_kInf = new ArrheniusRateConstant(L);
	lua_pop(L, 1);

	lua_getfield(L, -1, "k0");
	_k0 = new ArrheniusRateConstant(L);
	lua_pop(L, 1);

	// Look for an F_cent value.
	lua_getfield(L, -1, "F_cent");
	if ( !lua_isnil(L, -1) ) {
	    _Fcent_supplied = true;
	    _Fcent = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	else {
	    // Failng that, look for a, T1, T3 and possibly T2
	    _Fcent_supplied = false;
	    lua_pop(L, 1);
	    lua_getfield(L, -1, "a"); _a = luaL_checknumber(L, -1); lua_pop(L, 1);
	    lua_getfield(L, -1, "T1"); _T1 = luaL_checknumber(L, -1); lua_pop(L, 1);
	    lua_getfield(L, -1, "T3"); _T3 = luaL_checknumber(L, -1); lua_pop(L, 1);
	    lua_getfield(L, -1, "T2");
	    if ( !lua_isnumber(L, -1) ) {
		_T2_supplied = false;
		_T2 = 0.0;
	    }
	    else {
		_T2_supplied = true;
		_T2 = luaL_checknumber(L, -1);
	    }
	    lua_pop(L, 1);
	}
	
	_efficiencies = efficiencies.dup();
	_gmodel = gmodel;
    }
    TroeRateConstant dup()
    {
	return new TroeRateConstant(_kInf, _k0, _Fcent, _efficiencies, _gmodel);
    }
    override double eval(in GasState Q)
    {
	immutable double essentially_zero = 1.0e-30;
	double M = thirdBodyConcentration(Q, _efficiencies, _gmodel);
	double kInf = _kInf.eval(Q);
	double k0 = _k0.eval(Q);
	double p_r = k0*M/kInf;
	double log_p_r = log10(max(p_r, essentially_zero));
	double T = Q.Ttr;

	if ( !_Fcent_supplied ) {
	    _Fcent = (1.0 - _a)*exp(-T/_T3) + _a*exp(-T/_T1);
	    if ( _T2_supplied ) {
		_Fcent += exp(-_T2/T);
	    }
	}

	double log_F_cent = log10(max(_Fcent, essentially_zero));
	double c = -0.4 - 0.67*log_F_cent;
	double n = 0.75 - 1.27*log_F_cent;
	double d = 0.14;

	double numer = log_p_r + c; 
	double denom = n - d*numer;
	double frac = numer/denom;
	double log_F = log_F_cent / (1.0 + frac*frac);
	double F = pow(10,log_F);

	return F*k0*kInf*M/(k0*M + kInf);
    }
private:
    ArrheniusRateConstant _kInf, _k0;
    double _Fcent, _a, _T1, _T2, _T3;
    bool _Fcent_supplied, _T2_supplied;
    Tuple!(int, double)[] _efficiencies;
    GasModel _gmodel;
}


/++
 + Create a RateConstant object based on information in a LuaTable.
 +
 + Expected table format is:
 + tab = {model='Arrhenius',
 +        -- then model-specific parameters follow.
 +        A=..., n=..., C=...}
 +/
RateConstant createRateConstant(lua_State* L, Tuple!(int, double)[] efficiencies, GasModel gmodel)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "Arrhenius":
	return new ArrheniusRateConstant(L);
    case "Lindemann-Hinshelwood":
	return new LHRateConstant(L, efficiencies, gmodel);
    case "Troe":
	return new TroeRateConstant(L, efficiencies, gmodel);
    case "fromEqConst":
	return null;
    default:
	string msg = format("The rate constant model: %s could not be created.", model);
	throw new Exception(msg);
    }
}


version(rate_constant_test) {
    import util.msg_service;
    int main() {
	// Test 1. Rate constant for H2 + I2 reaction.
	auto rc = new ArrheniusRateConstant(1.94e14*1e-6, 0.0, 20620.0);
	auto gd = new GasState(1, 1);
	gd.Ttr = 700.0;
	assert(approxEqual(3.10850956e-5, rc.eval(gd), 1.0e-6), failedUnitTest());
	// Test 2. Read rate constant parameters for nitrogen dissociation
	// from Lua input and compute rate constant at 4000.0 K
	auto L = init_lua_State("sample-input/N2-diss.lua");
	lua_getglobal(L, "rate");
	auto rc2 = new ArrheniusRateConstant(L);
	gd.Ttr = 4000.0;
	assert(approxEqual(0.00159439, rc2.eval(gd), 1.0e-6), failedUnitTest());

	// Test 3. Try the pressure-dependent rate of Lindemann-Hinshelwood type
	

	return 0;
    }
}
