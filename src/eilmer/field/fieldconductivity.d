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
    Interestingly, as far as the electric potential cares its the gradient of the
    conductivity that matters, not the conductivity value itself. What this means
    is that if we want to solve a constant conductivity problem, we can just say
    that it is equal to one, and get the same answer as any other value.

    For the *current* though, things are different. At that stage we really do need
    the actual values; but for the moment calculating currents is a post processing
    thing only.
*/
    this() {}
    final const double opCall(const FlowState fs, const Vector3 pos){
        return 1.0;
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
    case "none":
        throw new Error("User has asked for solve_electric_field but failed to specify a conductivity model.");
    default:
        string errMsg = format("The conductivity model '%s' is not available.", name);
        throw new Error(errMsg);
    }
    return conductivity_model;
}
