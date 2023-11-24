/**
 * sutherland_therm_cond.d
 * Implements the Sutherland "law" for
 * thermal conductivity/
 *
 * Author: Rowan G. and Peter J.
 */

module gas.diffusion.sutherland_therm_cond;

import std.math;
import ntypes.complex;
import nm.number;
import gas.gas_model;
import gas.gas_state;
import gas.diffusion.therm_cond;
import util.lua;
import util.lua_service;
import util.msg_service;
/++
  Compute the thermal conductivity using Sutherland's expression.

  Params:
     T = temperature of gas in K
     T_ref = reference temperature in K
     k_ref = reference thermal conductivity in W/(m.K)
     S = Sutherland constant in K

   Returns:
     The thermal conductivity in W/(m.K)
+/
@nogc pure
number sutherland_thermal_conductivity(number T, double T_ref, double k_ref, double S)
in {
    debug {
        assert(T > 0.0, brokenPreCondition("temperature", __LINE__, __FILE__));
        assert(T_ref > 0.0, brokenPreCondition("T_ref", __LINE__, __FILE__));
        assert(k_ref > 0.0, brokenPreCondition("k_ref", __LINE__, __FILE__));
    }
}
out(result) {
    debug {
        assert(result > 0.0, brokenPostCondition("thermal conductivity", __LINE__, __FILE__));
    }
}
do{
    number k = k_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
    return k;
}

/++
  SutherlandThermCond is a thermal conductivity model.
+/
class SutherlandThermCond : ThermalConductivity {
public:
    this(in SutherlandThermCond src) {
        _T_ref = src._T_ref;
        _k_ref = src._k_ref;
        _S = src._S;
    }
    this(double T_ref, double k_ref, double S) {
        _T_ref = T_ref;
        _k_ref = k_ref;
        _S = S;
    }
    override SutherlandThermCond dup() const {
        return new SutherlandThermCond(this);
    }
    @nogc override number eval(number T) const {
        return sutherland_thermal_conductivity(T, _T_ref, _k_ref, _S);
    }
    @nogc override number eval(number T, number logT) const {
        return sutherland_thermal_conductivity(T, _T_ref, _k_ref, _S);
    }

private:
    double _T_ref;
    double _k_ref;
    double _S;
}

SutherlandThermCond createSutherlandThermalConductivity(lua_State* L)
{
    auto T_ref = getDouble(L, -1, "T_ref");
    auto k_ref = getDouble(L, -1, "k_ref");
    auto S = getDouble(L, -1, "S");
    return new SutherlandThermCond(T_ref, k_ref, S);
}

version(sutherland_therm_cond_test) {
    import std.conv;
    int main() {
        number T = 300.0;
        double T_ref = 273.0;
        double k_ref = 0.0241;
        double S = 194.0;
        assert(approxEqualNumbers(sutherland_thermal_conductivity(T, T_ref, k_ref, S), to!number(0.0262449),
                                  1.0e-6), failedUnitTest());

        auto tcm = new SutherlandThermCond(T_ref, k_ref, S);
        double k = tcm.eval(300.0);
        assert(approxEqualNumbers(k, to!number(0.0262449), 1.0e-6), failedUnitTest());

        return 0;
    }
}
