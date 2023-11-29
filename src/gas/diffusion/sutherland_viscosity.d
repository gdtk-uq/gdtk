/**
 * sutherland_visc.d
 * Implements the Sutherland "law" for
 * viscosity.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 *          2018-06-02 complex number flavour
 */

module gas.diffusion.sutherland_viscosity;

import std.math;
import ntypes.complex;
import nm.number;
import gas.gas_model;
import gas.gas_state;
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
@nogc pure number sutherland_viscosity(number T, double T_ref, double mu_ref, double S)
in {
    debug {
        assert(T > 0.0, brokenPreCondition("temperature", __LINE__, __FILE__));
        assert(T_ref > 0.0, brokenPreCondition("T_ref", __LINE__, __FILE__));
        assert(mu_ref > 0.0, brokenPreCondition("mu_ref", __LINE__, __FILE__));
    }
}
out(result) {
    debug {
        assert(result > 0.0, brokenPostCondition("viscosity", __LINE__, __FILE__));
    }
}
do{
    number mu = mu_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
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
    @nogc override number eval(number T) const {
        return sutherland_viscosity(T, _T_ref, _mu_ref, _S);
    }
    @nogc override number eval(number T, number logT) const {
        return sutherland_viscosity(T, _T_ref, _mu_ref, _S);
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

version(sutherland_viscosity_test) {
    import std.conv;
    int main() {
        number T = 300.0;
        double T_ref = 273.0;
        double mu_ref = 1.716e-5;
        double S = 111.0;
        assert(isClose(sutherland_viscosity(T, T_ref, mu_ref, S), 1.84691e-05, 1.0e-3), failedUnitTest());

        auto vm = new SutherlandViscosity(T_ref, mu_ref, S);
        double mu = vm.eval(300.0);
        assert(approxEqualNumbers(mu, to!number(1.84691e-05), 1.0e-5), failedUnitTest());

        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/O2-viscosity.lua");
        lua_getglobal(L, "Sutherland");
        vm = createSutherlandViscosity(L);
        lua_close(L);
        mu = vm.eval(300.0);
        assert(approxEqualNumbers(mu, to!number(1.84691e-05), 1.0e-3), failedUnitTest());

        return 0;
    }
}
