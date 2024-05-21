/**
 * solidblock.d
 *
 * Base class for a block representing a solid.
 * Typically, we want to compute the heat transfer
 * through the solid and its effect on the adjoining
 * flow field.
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-22-04
 *
 * Now a derived class from the Block base class
 * Kyle A. Damm 2020-02-11
 */

module solidblock;

import std.json;
import std.conv;

import util.lua;
import geom;
import globaldata;
import globalconfig;
import solidfvcell;
import solidfvinterface;
import solidbc;
import block;
import jacobian;

import nm.number;
import ntypes.complex;

import nm.smla;
import nm.bbla;

import lmr.solid.solidthermalmodel;

class SolidBlock : Block {
public:
    SolidThermalModel stm;
    double energyResidual; // monitor this for steady state
    Vector3 energyResidualLoc; // location of worst case
    int hncell; // number of history cells

    SolidFVCell[] cells; // collection of references to active cells in the domain
    SolidFVInterface[] faces; // collection of references to active faces in the domain
    SolidBoundaryCondition[] bc; // collection of references to boundary conditions

    FlowJacobian jacobian; // storage space for a Jacobian matrix

    this(int id, string label)
    {
        super(id, label);
    }

    override string toString() const { return "SolidBlock(id=" ~ to!string(id) ~ ")"; }

    void initSolidThermalModel(JSONValue jsonData)
    {
        auto model = jsonData["solid_block_" ~ to!string(id)]["model"].str;
        stm = lmr.solid.solidthermalmodel.initSolidThermalModel(jsonData["solid_thermal_models"][model]);
    }

    abstract void assembleArrays();
    abstract void initLuaGlobals();
    abstract void initBoundaryConditions(JSONValue jsonData);
    abstract void bindFacesAndVerticesToCells();
    abstract void bindCellsToFaces();
    abstract void assignCellLocationsForDerivCalc();
    abstract void readGrid(string filename);
    abstract void writeGrid(string filename, double sim_time);
    abstract void readSolution(string filename);
    abstract void computePrimaryCellGeometricData();
    abstract double determine_time_step_size(double cfl_value);

    abstract void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel);
    abstract void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel, SolidFVInterface f);
    abstract void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel);
    abstract void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel, SolidFVInterface f);
    abstract void applyPostFluxAction(double t, int tLevel);
    abstract void applyPostFluxAction(double t, int tLevel, SolidFVInterface f);
    abstract void computeSpatialDerivatives(int ftl);
    abstract void averageTemperatures();
    abstract void averageProperties();
    abstract void averageTGradients();
    abstract void computeFluxes();
    abstract void clearSources();
}

