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
        return -1.0*exp(pos.x)*cos(pos.y);
    }
}


ConductivityModel create_conductivity_model(string name){
    ConductivityModel conductivity_model;
    switch (name) {
    case "test":
        conductivity_model = new TestConductivity();
        break;
    default:
        string errMsg = format("The conductivity model '%s' is not available.", name);
        throw new Error(errMsg);
    }
    return conductivity_model;
}
