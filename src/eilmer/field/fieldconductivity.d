/**
 * fieldconductivity.d
 * Models for gas conductivity. Perhaps this should be a gas model
 *
 * Author: Nick Gibbons
 * Version: 2021-05-24
 */

module fieldconductivity;

import std.stdio;
import std.format;
import std.math;

import geom;
import flowstate;

interface ConductivityModel{
    const double opCall(const FlowState fs, const Vector3 pos);
}


class TestConductivity : ConductivityModel{
    this() {}
    final const double opCall(const FlowState fs, const Vector3 pos){
        return -1.0*exp(pos.x.re)*cos(pos.y.re);
    }
}

class ConstantConductivity : ConductivityModel{
/*
    Test with a just constant conductivity.
*/
    this() {}
    final const double opCall(const FlowState fs, const Vector3 pos){
        return 1.0;
    }
}

class RaizerConductivity : ConductivityModel{
/*
    Test with the formula from: Y. P. Razier, Gas Discharge Physics (Springer-Verlag, 1991)
     - Valid for air, nitrogen and argon when weakly ionised.
*/
    this() {}
    final const double opCall(const FlowState fs, const Vector3 pos){
        version(multi_T_gas) {
            double Tref = fs.gas.T_modes[0].re; // Hmmm. This will crash in single temp
        } else {
            double Tref = fs.gas.T.re;
        }
        double sigma = 8300.0*exp(-36000.0/Tref);
        return sigma;
    }
}


ConductivityModel create_conductivity_model(string name){
    ConductivityModel conductivity_model;
    switch (name) {
    case "test":
        conductivity_model = new TestConductivity();
        break;
    case "constant":
        conductivity_model = new ConstantConductivity();
        break;
    case "raizer":
        conductivity_model = new RaizerConductivity();
        break;
    case "none":
        throw new Error("User has asked for solve_electric_field but failed to specify a conductivity model.");
    default:
        string errMsg = format("The conductivity model '%s' is not available.", name);
        throw new Error(errMsg);
    }
    return conductivity_model;
}
