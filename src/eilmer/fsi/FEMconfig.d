module fsi.femconfig;

import std.json;
import std.stdio;

import util.json_helper;
import geom;
import fsi;

// We call the correct FEM model in the main code, the commented parts are from when
// I chose the model internally to the FSI operation. Keeping the code here for now,
// for posterity.
//enum PlateModel {
    //EulerBernoulli,
    //KirchhoffLove
//}

//PlateModel PlateModel_from_name(string name) {
    //switch (name) {
        //case "EulerBernoulli": return PlateModel.EulerBernoulli;
        //case "KirchhoffLove": return PlateModel.KirchhoffLove;
        //default: return PlateModel.EulerBernoulli;
    //}
//}

// Use the same enum template as in the main code, with the associated
// "from_name" methods.
enum ForcingType {
    Fluid, UserDefined
}

ForcingType ForcingType_from_name(string name) {
    switch (name) {
        case "Fluid": return ForcingType.Fluid;
        case "UserDefined": return ForcingType.UserDefined;
        default: return ForcingType.Fluid;
    }
}

enum DampingModel {
    None, Rayleigh2Parameter
}

DampingModel DampingModel_from_name(string name) {
    switch (name) {
        case "Rayleigh2Parameter": return DampingModel.Rayleigh2Parameter;
        default: return DampingModel.None;
    }
}

enum FEMTemporalScheme {
    euler, RK4
}

FEMTemporalScheme FEMTemporalScheme_from_name(string name) {
    switch (name) {
        case "euler": return FEMTemporalScheme.euler;
        case "RK4": return FEMTemporalScheme.RK4;
        // Include a catch for "rk4"?
        case "rk4": return FEMTemporalScheme.RK4;
        // Make RK4 the default because it's stable
        default: return FEMTemporalScheme.RK4;
    }
}

class FEMConfig {
public:
    //PlateModel plateModel;
    int Nx, Nz;
    ForcingType northForcing, southForcing;
    double length, width, density, thickness;
    double youngsModulus, poissonsRatio;
    string BCs;
    int couplingStep;
    double[3] plateNormal;
    bool quasi3D;
    bool writeMatrices;
    int[] movingBlks;
    FEMTemporalScheme temporalScheme;
    DampingModel dampingModel;
    double[] dampingRatios, naturalFrequencies;
    int[] historyNodes;

    this(string jobName) {
        string fileName = "config/" ~ jobName ~ ".fsi";
        JSONValue jsonData = readJSONfile(fileName);
        //this.plateModel = PlateModel_from_name(getJSONstring(jsonData, "model", "EulerBernoulli"));
        // The discretization
        this.Nx = getJSONint(jsonData, "Nx", 5);
        this.Nz = getJSONint(jsonData, "Nz", 0);
        // The forcing
        this.northForcing = ForcingType_from_name(getJSONstring(jsonData, "northForcing", "Fluid"));
        this.southForcing = ForcingType_from_name(getJSONstring(jsonData, "nouthForcing", "Fluid"));
        // Geometry
        this.length = getJSONdouble(jsonData, "length", 1.0);
        this.width = getJSONdouble(jsonData, "width", 1.0);
        // Material properties
        this.density = getJSONdouble(jsonData, "density", 8000.0);
        this.thickness = getJSONdouble(jsonData, "thickness", 1e-3);
        this.youngsModulus = getJSONdouble(jsonData, "youngsModulus", 190e9);
        this.poissonsRatio = getJSONdouble(jsonData, "poissonsRatio", 0.28);
        // Boundary conditions
        this.BCs = getJSONstring(jsonData, "BCs", "CF");
        // Orientation
        this.plateNormal = getJSONdoublearray(jsonData, "plateNormal", [0.0, 1.0, 0.0]);
        // Coupling to fluid
        this.couplingStep = getJSONint(jsonData, "couplingStep", 1);
        this.movingBlks = getJSONintarray(jsonData, "movingBlks", [-1]);
        // Whether to operate in quasi-3D mode
        this.quasi3D = getJSONbool(jsonData, "quasi3D", false);
        // Schemes for solid
        this.temporalScheme = FEMTemporalScheme_from_name(getJSONstring(jsonData, "temporalScheme", "RK4"));
        this.dampingModel = DampingModel_from_name(getJSONstring(jsonData, "dampingModel", "none"));
        // For Rayleigh2Parameter damping
        this.dampingRatios = getJSONdoublearray(jsonData, "dampingRatios", [0.001, 0.001]);
        this.naturalFrequencies = getJSONdoublearray(jsonData, "naturalFrequencies", [100.0, 100.0]);
        // IO
        this.writeMatrices = getJSONbool(jsonData, "writeMatrices", false);
        this.historyNodes = getJSONintarray(jsonData, "historyNodes", []);
    }
}
