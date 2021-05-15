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

import nm.number;
import nm.complex;

import nm.smla;
import nm.bbla;

class SolidBlock : Block {
public:
    double energyResidual; // monitor this for steady state
    Vector3 energyResidualLoc; // location of worst case
    int hncell; // number of history cells

    SolidFVCell[] activeCells; // collection of references to active cells in the domain
    SolidBoundaryCondition[] bc; // collection of references to boundary conditions

    version(nk_accelerator)
    {
        // Work-space for Newton-Krylov accelerator
        // These arrays and matrices are directly tied to using the
        // GMRES iterative solver.
        SMatrix!number JcT; // transposed Jacobian (w.r.t conserved variables)
        number maxRate, residuals;
        double normAcc, dotAcc;
        size_t nvars;
        Matrix!number Minv;
        number[] Fe, de, Dinv, r0, x0;
        number[] v, w, zed;
        number[] g0, g1;
        Matrix!number Q1;
        Matrix!number V;
    }

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
    abstract double determine_time_step_size(double cfl_value);

    abstract void applyPreSpatialDerivActionAtBndryFaces(double t, int tLevel);
    abstract void applyPreSpatialDerivActionAtBndryCells(double t, int tLevel);
    abstract void applyPostFluxAction(double t, int tLevel);
    abstract void computeSpatialDerivatives(int ftl);
    abstract void computeFluxes();
    abstract void clearSources();

    version(nk_accelerator) {
    void allocate_GMRES_workspace()
    {
        size_t mOuter = to!size_t(GlobalConfig.sssOptions.maxOuterIterations);
        size_t mInner = to!size_t(GlobalConfig.sssOptions.nInnerIterations);
        size_t n = 1*activeCells.length;
        nvars = n;
        // Now allocate arrays and matrices
        Fe.length = n;
        de.length = n; de[] = to!number(0.0);
        r0.length = n;
        x0.length = n;
        Dinv.length = n;
        v.length = n;
        w.length = n;
        zed.length = n;
        g0.length = mOuter+1;
        g1.length = mOuter+1;
        //h_outer.length = mOuter+1;
        //hR_outer.length = mOuter+1;
        V = new Matrix!number(n, mOuter+1);
        //H0_outer = new Matrix!number(mOuter+1, mOuter);
        //H1_outer = new Matrix!number(mOuter+1, mOuter);
        //Gamma_outer = new Matrix!number(mOuter+1, mOuter+1);
        //Q0_outer = new Matrix!number(mOuter+1, mOuter+1);
        Q1 = new Matrix!number(mOuter+1, mOuter+1);
    }
    }
}

