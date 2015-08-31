/**
 * sutherland_visc.d
 * Implements the Sutherland "law" for
 * viscosity.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

module gas.diffusion.sutherland_viscosity;

import std.math;
import gas.gas_model;
import gas.diffusion.viscosity;
import util.lua;
import util.lua_service;
import util.msg_service;

/++
  Compute the viscosity using Sutherland's expression.
  
  Params:
     T = temperature of gas in K
     T_ref = reference temperature in K
     mu_ref = reference viscosity in Pa.s
     S = Sutherland constant in K

   Returns:
     The viscosity in Pa.s.
+/
pure double sutherland_viscosity(double T, double T_ref, double mu_ref, double S)
in {
    assert(T > 0.0, brokenPreCondition("temperature", __LINE__, __FILE__));
    assert(T_ref > 0.0, brokenPreCondition("T_ref", __LINE__, __FILE__));
    assert(mu_ref > 0.0, brokenPreCondition("mu_ref", __LINE__, __FILE__));
}
out(result) {
    assert(result > 0.0, brokenPostCondition("viscosity", __LINE__, __FILE__));
}
body{
    double mu = mu_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
    return mu;
}

/++
  SutherlandViscosity is a viscosity model.
+/
class SutherlandViscosity : Viscosity {
public:
    this(in SutherlandViscosity src) {
	_T_ref = src._T_ref;
	_mu_ref = src._mu_ref;
	_S = src._S;
    }
    this(double T_ref, double mu_ref, double S) {
	_T_ref = T_ref;
	_mu_ref = mu_ref;
	_S = S;
    }
    override SutherlandViscosity dup() const {
	return new SutherlandViscosity(this);
    }
    /++
      Compute the viscosity assuming that temperature is
      up-to-date in GasState Q.
    +/
    override double eval(in GasState Q) const {
	return sutherland_viscosity(Q.T[0], _T_ref, _mu_ref, _S);
    }

private:
    double _T_ref;
    double _mu_ref;
    double _S;
}

SutherlandViscosity createSutherlandViscosity(lua_State* L)
{
    auto T_ref = getDouble(L, -1, "T_ref");
    auto mu_ref = getDouble(L, -1, "mu_ref");
    auto S = getDouble(L, -1, "S");
    return new SutherlandViscosity(T_ref, mu_ref, S);
}

unittest {
    double T = 300.0;
    double T_ref = 273.0; 
    double mu_ref = 1.716e-5;
    double S = 111.0;
    assert(approxEqual(sutherland_viscosity(T, T_ref, mu_ref, S), 1.84691e-05), failedUnitTest(__LINE__, __FILE__));

    auto vm = new SutherlandViscosity(T_ref, mu_ref, S);
    auto gd = new GasState(1, 1);
    gd.T[0] = 300.0;
    vm.update_viscosity(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), failedUnitTest(__LINE__, __FILE__));

    lua_State* L = init_lua_State("sample-data/O2-viscosity.lua");
    lua_getglobal(L, "Sutherland");
    vm = createSutherlandViscosity(L);
    lua_close(L);
    vm.update_viscosity(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), failedUnitTest(__LINE__, __FILE__));
}
