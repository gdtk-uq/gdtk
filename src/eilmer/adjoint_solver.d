/** adjoint_solver.d
 * Code to construct, and solve the adjoint equations.
 *
 * Author: Kyle D.
 * Date: 2017-09-18
 *
 *
 * History:
 *   2017-09-18 : started -- 2D, 1st order, single-block, structured grid, Euler solver.
**/

import core.stdc.stdlib : exit;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.array;
import std.math;
import std.datetime;

import nm.smla;
import nm.bbla;
import nm.rsla;

import fluxcalc;
import grid;
import sgrid;
import usgrid;
import block;
import sblock;
import ublock;
import fvcell;
import fvinterface;
import sblock;
import globaldata;
import globalconfig;
import simcore;
import fvcore;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess;
import loads;

import std.math;
import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.regex;
import std.algorithm;
import std.bitmanip;
import std.stdint;
import std.range;
import gzip;
import fvcore;
import fileutil;
import geom;
import sgrid;
import grid;
import gas;
import globalconfig;
import flowsolution;
import solidsolution;

// EPSILON parameter for numerical differentiation of flux jacobian
// Value used based on Vanden and Orkwis (1996), AIAA J. 34:6 pp. 1125-1129
immutable double EPSILON = 1.0e-2;
immutable double ESSENTIALLY_ZERO = 1.0e-15;

void main(string[] args) {

    writeln("Eilmer compressible-flow simulation code -- adjoint solver.");
    
    // -----------------------------------------------------
    // 1. Read in flow solution
    // -----------------------------------------------------

    string msg = "Usage:                              Comment:\n";
    msg       ~= "e4adjoint  [--job=<string>]     name of job\n";
    
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    try {
	getopt(args,
	       "job", &jobName,
	       );
    } catch (Exception e) {
	writeln("Problem parsing command-line options.");
	writeln("Arguments not processed:");
	args = args[1 .. $]; // Dispose of program in first arg
	foreach (arg; args) writeln("   arg: ", arg);
	write(msg);
	exit(1);
    }

    GlobalConfig.base_file_name = jobName;
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    auto last_tindx = tindx_list[$-1];

    int maxCPUs = 1;
    int maxWallClock = 5*24*3600; // 5 days default
    init_simulation(last_tindx, maxCPUs, maxWallClock);
    //auto soln = new FlowSolution(jobName, ".", last_tindx, GlobalConfig.nBlocks);
    
    writeln("simulation initialised");
    
    // -----------------------------------------------------
    // 2. store the stencil of effected cells for each cell
    // -----------------------------------------------------

    // NB: currently only for 1st order interpolation
    // TODO: high order interpolation
    
    FVCell[][] stencils;
    foreach (blk; gasBlocks) {
	foreach(i, cell; blk.cells) {
	    FVCell[] cell_refs;
	    cell_refs ~= cell; // add the parent cell as the first reference
	    foreach(f; cell.iface) {
		if (f.left_cell.id != cell.id && f.left_cell.id < 1000000) cell_refs ~= f.left_cell;
		if (f.right_cell.id != cell.id && f.right_cell.id < 1000000) cell_refs ~= f.right_cell;
	    }
	    stencils ~= cell_refs;
	}
    }

    // -----------------------------------------------------
    // 3. Compute and store perturbed flux
    // -----------------------------------------------------
    FVCell cellPp; FVCell cellPm; FVCell cellR; FVCell cellL;
    double h; double diff;
    FVInterface ifacePp;
    FVInterface ifacePm;

    double[][] Jac;
    size_t nc = 3; // number of primitive variables
    size_t ncells = gasBlocks[0].ncells;
    // currently stores the entire Jacobian -- this is quite wasteful
    // TODO: sparse matrix storage
    foreach (i; 0..nc*ncells) {
	double[] row;
	foreach (j; 0..nc*ncells) {
	    row ~= 0.0;
	}
	Jac ~= row;
    }
    // hard-coded for perturbing rho, u, v, P
    foreach (blk; gasBlocks) {
	foreach(ci, cell; blk.cells) {
	    cellPp = new FVCell(dedicatedConfig[blk.id]);
	    ifacePp = new FVInterface(dedicatedConfig[blk.id], false);
	    cellPm = new FVCell(dedicatedConfig[blk.id]);
	    ifacePm = new FVInterface(dedicatedConfig[blk.id], false);
	    
	    // -----------------------------------------------------
	    // perturb rho
	    // -----------------------------------------------------
	    h = cell.fs.gas.rho * EPSILON + EPSILON;
	    cellPm.copy_values_from(cell, CopyDataOption.all);
	    cellPm.fs.gas.rho -= h;
	    cellPp.copy_values_from(cell, CopyDataOption.all);
	    cellPp.fs.gas.rho += h;
	    
	    foreach(iface; cell.iface) {
		// - perturbation
		// update thermo
		blk.myConfig.gmodel.update_thermo_from_rhop(cellPm.fs.gas);
		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.rho -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.rho += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}

		if(iface.left_cell.id == cellPm.id) {
		    cellR = iface.right_cell;
		    cellL = cellPm;
		}
		else {
		    cellR = cellPm;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.rho -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.rho += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		ifacePm.copy_values_from(iface, CopyDataOption.all);
		
		// + perturbation
		// update thermo
		blk.myConfig.gmodel.update_thermo_from_rhop(cellPp.fs.gas);
		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.rho += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.rho -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		
		if(iface.left_cell.id == cellPp.id) {
		    cellR = iface.right_cell;
		    cellL = cellPp;
		}
		else {
		    cellR = cellPp;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.rho += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.rho -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		ifacePp.copy_values_from(iface, CopyDataOption.all);
		
		diff = ifacePp.F.mass - ifacePm.F.mass;
		iface.dFdU[0][0] = diff/(2.0*h);
	    
		diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;
		iface.dFdU[1][0] = diff/(2.0*h);

		//diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;
		//iface.dFdU[2][0] = diff/(2.0*h);

		diff = ifacePp.F.total_energy - ifacePm.F.total_energy;
		//writef("iface.id = %d", ifacePp.id);
		//writef(", diff/(2h) = %.16f \n", diff/(2.0*h));
		//writef(", perturbed flux = %.16f", ifacePp.F.total_energy);
		//writef(", unperturbed flux = %.16f \n", ifacePm.F.total_energy);
		iface.dFdU[2][0] = diff/(2.0*h);
	    }
	    // -----------------------------------------------------
	    // perturb u
	    // -----------------------------------------------------

	    h = cell.fs.vel.refx * EPSILON + EPSILON;
	    cellPm.copy_values_from(cell, CopyDataOption.all);
	    cellPm.fs.vel.refx -= h;
	    cellPp.copy_values_from(cell, CopyDataOption.all);
	    cellPp.fs.vel.refx += h;
	    foreach(iface; cell.iface) {
		// - perturbation
		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.vel.refx -= h;
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.vel.refx += h;
		}
		
		if(iface.left_cell.id == cellPm.id) {
		    cellR = iface.right_cell;
		    cellL = cellPm;
		}
		else {
		    cellR = cellPm;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.vel.refx -= h;
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.vel.refx += h;
		}
		ifacePm.copy_values_from(iface, CopyDataOption.all);
		// + perturbation
		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.vel.refx += h;
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.vel.refx -= h;
		}
		
		if(iface.left_cell.id == cellPp.id) {
		    cellR = iface.right_cell;
		    cellL = cellPp;
		}
		else {
		    cellR = cellPp;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.vel.refx += h;
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.vel.refx -= h;
		}
		ifacePp.copy_values_from(iface, CopyDataOption.all);
		
		diff = ifacePp.F.mass - ifacePm.F.mass;
		iface.dFdU[0][1] = diff/(2.0*h);
	    
		diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;
		iface.dFdU[1][1] = diff/(2.0*h);

		//diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;
		//iface.dFdU[2][1] = diff/(2.0*h);
	    
		diff = ifacePp.F.total_energy - ifacePm.F.total_energy;
		//writef("iface.id = %d", ifacePp.id);
		//writef(", diff/(2h) = %.16f \n", diff/(2.0*h));
		//writef(", perturbed flux = %.16f", ifacePp.F.total_energy);
		//writef(", unperturbed flux = %.16f \n", ifacePm.F.total_energy);
		iface.dFdU[2][1] = diff/(2.0*h);
	    }
	    // -----------------------------------------------------
	    // perturb v
	    // -----------------------------------------------------
	    /*
	    h = cell.fs.vel.refy * EPSILON + EPSILON;
	    cellPm.copy_values_from(cell, CopyDataOption.all);
	    cellPm.fs.vel.refy -= h;
	    cellPp.copy_values_from(cell, CopyDataOption.all);
	    cellPp.fs.vel.refy += h;

	    foreach(iface; cell.iface) {
		// - perturbation
		// apply bcs
		if (iface.is_on_boundary) {
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		}
		
		ifacePm.copy_values_from(iface, CopyDataOption.all);
		if(iface.left_cell.id == cellPm.id) {
		    cellR = iface.right_cell;
		    cellL = cellPm;
		}
		else {
		    cellR = cellPm;
		    cellL = iface.left_cell;
		}
		
		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, ifacePm, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    ifacePm.copy_values_from(iface, CopyDataOption.all);
		}
		
		// + perturbation
		// apply bcs
		if (iface.is_on_boundary) {
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		}

		ifacePp.copy_values_from(iface, CopyDataOption.all);
		if(iface.left_cell.id == cellPp.id) {
		    cellR = iface.right_cell;
		    cellL = cellPp;
		}
		else {
		    cellR = cellPp;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, ifacePp, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    ifacePp.copy_values_from(iface, CopyDataOption.all);
		}
		
		diff = ifacePp.F.mass - ifacePm.F.mass;
		iface.dFdU[0][2] = diff/(2.0*h);
		
		diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;
		iface.dFdU[1][2] = diff/(2.0*h);

		diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;
		iface.dFdU[2][2] = diff/(2.0*h);
	    
		diff = ifacePp.F.total_energy - ifacePm.F.total_energy;
		iface.dFdU[3][2] = diff/(2.0*h);
	    }
	    */
	    // -----------------------------------------------------
	    // perturb p
	    // -----------------------------------------------------

	    h = cell.fs.gas.p * EPSILON + EPSILON;
	    cellPm.copy_values_from(cell, CopyDataOption.all);
	    cellPm.fs.gas.p -= h;
	    cellPp.copy_values_from(cell, CopyDataOption.all);
	    cellPp.fs.gas.p += h;
	    foreach(iface; cell.iface) {
		// - perturbation
		// update thermo
		blk.myConfig.gmodel.update_thermo_from_rhop(cellPm.fs.gas);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.p -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.p += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		
		if(iface.left_cell.id == cellPm.id) {
		    cellR = iface.right_cell;
		    cellL = cellPm;
		}
		else {
		    cellR = cellPm;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.p -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.p += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		ifacePm.copy_values_from(iface, CopyDataOption.all);
			
		// + perturbation
		// update thermo
		blk.myConfig.gmodel.update_thermo_from_rhop(cellPp.fs.gas);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.p += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.p -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		
		if(iface.left_cell.id == cellPp.id) {
		    cellR = iface.right_cell;
		    cellL = cellPp;
		}
		else {
		    cellR = cellPp;
		    cellL = iface.left_cell;
		}

		blk.Lft.copy_values_from(cellL.fs);
		blk.Rght.copy_values_from(cellR.fs);
		compute_interface_flux(blk.Lft, blk.Rght, iface, blk.myConfig, blk.omegaz);

		// apply bcs
		if (iface.is_on_boundary) {
		    cell.fs.gas.p += h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		    blk.applyPreReconAction(0.0, 0, 0); // assume sim_time = 0.0, gtl = 0, ftl = 0
		    blk.applyPostConvFluxAction(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryFaces(0.0, 0, 0);
		    blk.applyPreSpatialDerivActionAtBndryCells(0.0, 0, 0);
		    blk.applyPostDiffFluxAction(0.0, 0, 0);
		    cell.fs.gas.p -= h;
		    blk.myConfig.gmodel.update_thermo_from_rhop(cell.fs.gas);
		}
		ifacePp.copy_values_from(iface, CopyDataOption.all);

		diff = ifacePp.F.mass - ifacePm.F.mass;
		iface.dFdU[0][2] = diff/(2.0*h);
	    
		diff = ifacePp.F.momentum.x - ifacePm.F.momentum.x;
		iface.dFdU[1][2] = diff/(2.0*h);

		//diff = ifacePp.F.momentum.y - ifacePm.F.momentum.y;
		//iface.dFdU[2][3] = diff/(2.0*h);
	    
		diff = ifacePp.F.total_energy - ifacePm.F.total_energy;
		iface.dFdU[2][2] = diff/(2.0*h);
	    }
	    // -----------------------------------------------------
	    // loop through influenced cells and fill out Jacobian 
	    // -----------------------------------------------------
	    // at this point we can use the cell counter ci to access
	    // the correct stencil
	    foreach(c; stencils[ci]) {
		size_t I, J; // indices in Jacobian matrix
		double integral;
		double volInv = 1.0 / c.volume[0];
		for ( size_t ic = 0; ic < nc; ++ic ) {
		    I = c.id*nc + ic; // row index
		    for ( size_t jc = 0; jc < nc; ++jc ) {
			integral = 0.0;
			J = cell.id*nc + jc; // column index
			//writeln("--------------------------------");
			foreach(fi, iface; c.iface) {
			    //writeln(iface.dFdU);
			    //writeln("iface.id = ", iface.id, ", pcell.id = ", cell.id, ", icell.id = ", c.id, ", ", iface.dFdU[ic][jc], ", 1/vol = ", volInv, ", iface.area", iface.area[0]);
			    integral -= c.outsign[fi] * iface.dFdU[ic][jc] * iface.area[0]; // gtl=0
			}
			//writeln(I, ", ", J, ", ", volInv * integral);
			//writeln("--------------------------------");
			Jac[I][J] = volInv * integral;
		    }
		}
	    }
	    // clear the flux Jacobian entries
	    foreach (iface; cell.iface) {
		foreach (i; 0..iface.dFdU.length) {
		    foreach (j; 0..iface.dFdU[i].length) {
			iface.dFdU[i][j] = 0.0;
		    }
		}
	    }
	}
    }

    writeln(Jac.length);
    //--------------------------------------------------------
    // Transpose Jac
    //--------------------------------------------------------
    /*
    Jac[297][297] = -1386.69;
    Jac[297][298] = -11.5109;
    Jac[297][299] = -0.00748001;
    Jac[298][297] = -192213;
    Jac[298][298] = -8342.65;
    Jac[298][299] = -5.49689;
    Jac[299][297] = -1.34348e+07;
    Jac[299][298] = -6.26195e+06;
    Jac[299][299] = -8717.94;
    */
    writeln("J = ", Jac);
    
    double[][] JacT;
    foreach (i; 0..nc*ncells) {
	double[] row;
	foreach (j; 0..nc*ncells) {
	    row ~= Jac[j][i];
	}
	JacT ~= row;
    }
    //writeln("JT = ", JacT);
    writeln("-------------------------------------------------------------------------------------");
    //writeln("JacT = ", JacT);
    // -----------------------------------------------------
    //  Form cost function sensitvity
    // -----------------------------------------------------
    // Analytically form dJdV by hand differentiation
    // cost function is defined as: J(Q) = 0.5*integral[0->l] (p-p*)^2
    double[] dJdV;
    double[] p_target;

    // target pressure distribution saved in file target.dat
    auto file = File("target.dat", "r");
    foreach(i; 0 .. ncells) {
	auto lineContent = file.readln().strip();
	auto tokens = lineContent.split();
	p_target ~= to!double(tokens[8]);
    }
    writeln("target pressure imported");
    foreach (blk; gasBlocks) {
	foreach(i, cell; blk.cells) {
	    dJdV ~= 0.0;
	    dJdV ~= 0.0;
	    //dJdV ~= 0.0;
	    dJdV ~= 0.5*(2.0*cell.fs.gas.p - 2.0*p_target[i]);
	}
    }

    // -----------------------------------------------------
    // Solve adjoint equations
    // -----------------------------------------------------

    // form augmented matrix aug = [A|B] = [Jac|dJdQ]
    writeln("-------------------------------------------------------------");
    writeln("dJdV = ", dJdV);
    size_t ncols = nc*ncells*2;
    size_t nrows = nc*ncells;
    double[600][300] aug;
    foreach (i; 0 .. JacT.length) {
	foreach (j; 0 .. (Jac[i].length*2) ) {
	    if (j < JacT[i].length) aug[i][j] =  JacT[i][j];
	    else if (j-300 == i) aug[i][j] = 1.0;
	    else aug[i][j] = 0.0;
	}
    }

    //writeln(aug);
    //writeln("aug = ", aug);

    // solve for adjoint variables
    //gaussJordanElimination(aug);
    computeInverse!(300,300)(aug, 1.0e-16);

    //writeln(aug);
    
    //writeln("sol = ", aug);
    
    double[300] psi;
    //foreach (i; 0 .. nrows) {
    //	psi ~= aug[i][ncols-1];
    //}
    foreach (i; 0 .. 100) {
	double temp1 = 0.0; double temp2 = 0.0; double temp3 = 0.0;
	foreach (j; 0 .. JacT.length) {
	    temp1 += aug[i*3][j+300] * -dJdV[j];
	    temp2 += aug[i*3+1][j+300] * -dJdV[j];
	    temp3 += aug[i*3+2][j+300] * -dJdV[j];
	}
	psi[i*3] = temp1;
	psi[i*3+1] = temp2;
	psi[i*3+2] = temp3;
    }
    
    writeln("psi = ", psi);

    double[3*100] temp;
    foreach (i; 0 .. 100) {
	double temp1 = 0.0; double temp2 = 0.0; double temp3 = 0.0;
	foreach (j; 0 .. JacT.length) {
	    temp1 += JacT[i*3][j] * psi[j];
	    temp2 += JacT[i*3+1][j] * psi[j];
	    temp3 += JacT[i*3+2][j] * psi[j];
	}
	temp[i*3] = temp1;
	temp[i*3+1] = temp2;
	temp[i*3+2] = temp3;
    }
    //writeln(temp);
    
    foreach(i; 0 .. ncells) {
	FVCell cell = gasBlocks[0].cells[i];
	auto writer = format("%f %f %f %f \n", cell.pos[0].x, psi[i*3], psi[i*3+1], psi[i*3+2]);
	append("e4_adjoint_vars.dat", writer);
    }
    
    writeln("Done simulation.");
}

