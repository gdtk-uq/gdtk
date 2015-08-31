/**
 * sutherland_therm_cond.d
 * Implements the Sutherland "law" for
 * thermal conductivity/
 *
 * Author: Rowan G. and Peter J.
 */

module gas.diffusion.sutherland_therm_cond;

import std.math;
import gas.gas_model;
import gas.diffusion.therm_cond;
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
pure double sutherland_thermal_conductivity(double T, double T_ref, double k_ref, double S)
in {
    assert(T > 0.0, brokenPreCondition("temperature", __LINE__, __FILE__));
    assert(T_ref > 0.0, brokenPreCondition("T_ref", __LINE__, __FILE__));
    assert(k_ref > 0.0, brokenPreCondition("k_ref", __LINE__, __FILE__));
}
out(result) {
    assert(result > 0.0, brokenPostCondition("thermal conductivity", __LINE__, __FILE__));
}
body{
    double k = k_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
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
    override double eval(in GasState Q, int imode) const {
	return sutherland_thermal_conductivity(Q.T[imode], _T_ref, _k_ref, _S);
    }

private:
    double _T_ref;
    double _k_ref;
    double _S;
}

unittest {
    double T = 300.0;
    double T_ref = 273.0; 
    double k_ref = 0.0241;
    double S = 194.0;
    assert(approxEqual(sutherland_thermal_conductivity(T, T_ref, k_ref, S), 0.0262449), failedUnitTest(__LINE__, __FILE__));

    auto tcm = new SutherlandThermCond(T_ref, k_ref, S);
    auto gd = new GasState(1, 1);
    gd.T[0] = 300.0;
    gd.k[0] = 0.0;
    tcm.update_thermal_conductivity(gd);
    assert(approxEqual(gd.k[0], 0.0262449), failedUnitTest(__LINE__, __FILE__));
}
