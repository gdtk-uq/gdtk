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
import solidbc;
import solidprops;
import block;

class SolidBlock : Block {
public:
    double energyResidual; // monitor this for steady state
    Vector3 energyResidualLoc; // location of worst case
    int hncell; // number of history cells

    SolidFVCell[] activeCells; // collection of references to active cells in the domain
    SolidBoundaryCondition[] bc; // collection of references to boundary conditions

    this(int id, string label)
    {
        super(id, label);
        myConfig = dedicatedConfig[id];
    }

    override string toString() const { return "SolidBlock(id=" ~ to!string(id) ~ ")"; }
    abstract void assembleArrays();
    abstract void initLuaGlobals();
    abstract void initBoundaryConditions(JSONValue jsonData);
    abstract void bindFacesAndVerticesToCells();
    abstract void assignCellLocationsForDerivCalc();
    abstract void readGrid(string filename);
    abstract void writeGrid(string filename, double sim_time);
    abstract void readSolution(string filename);
    abstract void writeSolution(string fileName, double simTime);
    abstract void computePrimaryCellGeometricData();

    abstract void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel);
    abstract void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel);
    abstract void applyPostFluxAction(double t, int tLevel);
    abstract void computeSpatialDerivatives(int ftl);
    abstract void computeFluxes();
    abstract void clearSources();
}

