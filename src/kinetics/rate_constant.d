/**
 * rate_constant.d
 * Implements the calculation of reaction rate coefficients.
 *
 * Author: Rowan G.
 */

module kinetics.rate_constant;

import std.math;
import std.string;
import util.lua;
import util.lua_service;
import gas;

/++
  RateConstant is a class for computing the rate constant of
  a chemical reaction.
+/

interface RateConstant {
    RateConstant dup() const;
    double eval(in GasState Q) const;
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
    ArrheniusRateConstant dup() const
    {
	return new ArrheniusRateConstant(_A, _n, _C);
    }
    override double eval(in GasState Q) const
    {
	double T = Q.T[0];
	return _A*pow(T, _n)*exp(-_C/T);
    }
private:
    double _A, _n, _C;
}

/++
 + Create a RateConstant object based on information in a LuaTable.
 +
 + Expected table format is:
 + tab = {model='Arrhenius',
 +        -- then model-specific parameters follow.
 +        A=..., n=..., C=...}
 +/
RateConstant createRateConstant(lua_State* L)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "Arrhenius":
	return new ArrheniusRateConstant(L);
    default:
	string msg = format("The rate constant model: %s could not be created.", model);
	throw new Exception(msg);
    }
}


unittest {
    import util.msg_service;
    // Test 1. Rate constant for H2 + I2 reaction.
    auto rc = new ArrheniusRateConstant(1.94e14*1e-6, 0.0, 20620.0);
    auto gd = new GasState(1, 1);
    gd.T[0] = 700.0;
    assert(approxEqual(3.10850956e-5, rc.eval(gd)), failedUnitTest());
    // Test 2. Read rate constant parameters for nitrogen dissociation
    // from Lua input and compute rate constant at 4000.0 K
    auto L = init_lua_State("sample-input/N2-diss.lua");
    lua_getglobal(L, "rate");
    auto rc2 = new ArrheniusRateConstant(L);
    gd.T[0] = 4000.0;
    assert(approxEqual(0.00159439, rc2.eval(gd)), failedUnitTest());

}
