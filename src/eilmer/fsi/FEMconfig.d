module fsi.femconfig;

import std.json;
import std.stdio;

import json_helper;
import geom;
import fsi;

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
    DampingModel dampingModel;
    double dampingRatio;
    int[] historyNodes;

    this(string jobName) {
        string fileName = "config/" ~ jobName ~ ".fsi";
        JSONValue jsonData = readJSONfile(fileName);
        //this.plateModel = PlateModel_from_name(getJSONstring(jsonData, "model", "EulerBernoulli"));
        this.Nx = getJSONint(jsonData, "Nx", 5);
        this.Nz = getJSONint(jsonData, "Nz", 0);
        this.northForcing = ForcingType_from_name(getJSONstring(jsonData, "northForcing", "Fluid"));
        this.southForcing = ForcingType_from_name(getJSONstring(jsonData, "nouthForcing", "Fluid"));
        this.length = getJSONdouble(jsonData, "length", 1.0);
        this.width = getJSONdouble(jsonData, "width", 1.0);
        this.density = getJSONdouble(jsonData, "density", 8000.0);
        this.thickness = getJSONdouble(jsonData, "thickness", 1e-3);
        this.youngsModulus = getJSONdouble(jsonData, "youngsModulus", 190e9);
        this.poissonsRatio = getJSONdouble(jsonData, "poissonsRatio", 0.28);
        this.BCs = getJSONstring(jsonData, "BCs", "CF");
        this.plateNormal = getJSONdoublearray(jsonData, "plateNormal", [0.0, 1.0, 0.0]);
        this.couplingStep = getJSONint(jsonData, "couplingStep", 1);
        this.quasi3D = getJSONbool(jsonData, "quasi3D", false);
        this.writeMatrices = getJSONbool(jsonData, "writeMatrices", false);
        this.movingBlks = getJSONintarray(jsonData, "movingBlks", [-1]);
        this.dampingModel = DampingModel_from_name(getJSONstring(jsonData, "dampingModel", "none"));
        this.dampingRatio = getJSONdouble(jsonData, "dampingRatio", 0.001);
        this.historyNodes = getJSONintarray(jsonData, "historyNodes", []);
    }
}
