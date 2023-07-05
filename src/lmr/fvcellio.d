/**
 * A class to access data in a finite-volume cell for use with I/O functions.
 *
 * Author: Rowan J. Gollan
 * Date: 2023-06-20
 */

module fvcellio;

import std.stdio;
import std.algorithm : startsWith;
import std.string : split;
import std.conv : to;

import nm.number;

import globalconfig;
import lmrexceptions : LmrException;
import fvcell : FVCell;

/**
 * Build the list of variables for the flow field, based on modelling configuation options.
 *
 * The string names for variables should match those used later in this module in:
 *
 *    FVCellIO.opIndex; and
 *    FVCellIO.opIndexAssign.
 *
 * Authors: RJG
 */
string[] buildFlowVariables()
{
    alias cfg = GlobalConfig;
    string[] variables;
    variables ~= "pos.x";
    variables ~= "pos.y";
    if (cfg.dimensions == 3) variables ~= "pos.z";
    variables ~= "vol";
    variables ~= "rho";
    variables ~= "vel.x";
    variables ~= "vel.y";
    if (cfg.dimensions == 3) variables ~= "vel.z";
    if (cfg.MHD) {
	variables ~= "B.x";
	variables ~= "B.y";
	if (cfg.dimensions == 3) variables ~= "B.z";
    }
    variables ~= "p";
    variables ~= "a";
    if (cfg.viscous) {
	variables ~= "mu";
	variables ~= "k";
	foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
	    variables ~= "k-" ~ cfg.gmodel_master.energy_mode_name(imode);
	}
    }
    if (cfg.turbulence_model_name != "none") {
	variables ~= "mu_t";
	variables ~= "k_t";
	foreach (iturb; 0 .. cfg.turb_model.nturb) {
	    variables ~= "tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
	}
    }
    variables ~= "shock-detector";
    if (cfg.gmodel_master.n_species > 1) {
	foreach (isp; 0 .. cfg.gmodel_master.n_species) {
	    variables ~= "massf-" ~ cfg.gmodel_master.species_name(isp);
	}
    }
    if (cfg.gmodel_master.n_species > 1) variables ~= "dt_subcycle";
    variables ~= "e";
    variables ~= "T";
    foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
	variables ~= "e-" ~ cfg.gmodel_master.energy_mode_name(imode);
	variables ~= "T-" ~ cfg.gmodel_master.energy_mode_name(imode);
    }
    if (cfg.with_local_time_stepping) variables ~= "dt_local";

    return variables;

}

class FVCellIO {
public:

    this(string[] variables)
    {
	alias cfg = GlobalConfig;
	foreach (var; variables) {
	    if (var.startsWith("massf-")) {
		auto spName = var.split("-")[1];
		auto spIdx = cfg.gmodel_master.species_index(spName);
		if (spIdx != -1)
		    mSpecies[var.idup] = spIdx;
		else
		    throw new LmrException("Invalid species name: " ~ spName);
	    }
	    if (var.startsWith("T-")) {
		auto modeName = var.split("-")[1];
		auto modeIdx = cfg.gmodel_master.energy_mode_index(modeName);
		if (modeIdx != -1)
		    mTModes[var.idup] = modeIdx;
		else
		    throw new LmrException("Invalid mode name: " ~ modeName);
	    }
	    if (var.startsWith("e-")) {
		auto modeName = var.split("-")[1];
		auto modeIdx = cfg.gmodel_master.energy_mode_index(modeName);
		if (modeIdx != -1)
		    muModes[var.idup] = modeIdx;
		else
		    throw new LmrException("Invalid mode name: " ~ modeName);
	    }
	    if (var.startsWith("k-")) {
		auto modeName = var.split("-")[1];
		auto modeIdx = cfg.gmodel_master.energy_mode_index(modeName);
		if (modeIdx != -1)
		    mkModes[var.idup] = modeIdx;
		else
		    throw new LmrException("Invalid mode name: " ~ modeName);
	    }
	    if (var.startsWith("tq-")) {
		auto tqName = var.split("-")[1];
		auto tqIdx = cfg.turb_model.primitive_variable_index(tqName);
		if (tqIdx != -1)
		    mTurbQuants[var.idup] = tqIdx;
		else
		    throw new LmrException("Invalid turbulence quantity name: " ~ tqName);
	    }
	}
    }

    const double opIndex(FVCell cell, string var)
    {
	// First handle array-stored values:
	// species, modes and turbulence quantities
	if (var in mSpecies) {
	    return cell.fs.gas.massf[mSpecies[var]].re;
	}
	if (var in mTModes) {
	    return cell.fs.gas.T_modes[mTModes[var]].re;
	}
	if (var in muModes) {
	    return cell.fs.gas.u_modes[muModes[var]].re;
	}
	if (var in mkModes) {
	    return cell.fs.gas.k_modes[mkModes[var]].re;
	}
	if (var in mTurbQuants) {
	    return cell.fs.turb[mTurbQuants[var]].re;
	}

	// For everything else, find by cases
	switch (var) {
	case "pos.x": return cell.pos[0].x.re;
	case "pos.y": return cell.pos[0].y.re;
	case "pos.z": return cell.pos[0].z.re;
	case "vol": return cell.volume[0].re;
	case "rho": return cell.fs.gas.rho.re;
	case "vel.x": return cell.fs.vel.x.re;
	case "vel.y": return cell.fs.vel.y.re;
	case "vel.z": return cell.fs.vel.z.re;
	case "B.x": return cell.fs.B.x.re;
	case "B.y": return cell.fs.B.y.re;
	case "B.z": return cell.fs.B.z.re;
	case "p": return cell.fs.gas.p.re;
	case "a": return cell.fs.gas.a.re;
	case "mu": return cell.fs.gas.mu.re;
	case "k": return cell.fs.gas.k.re;
	case "mu_t": return cell.fs.mu_t.re;
	case "k_t": return cell.fs.k_t.re;
	case "shock-detector": return cell.fs.S.re;
	case "dt_subcycle": return cell.dt_chem.re;
	case "dt_local": return cell.dt_local.re;
	case "e": return cell.fs.gas.u.re;
	case "T": return cell.fs.gas.T.re;
	default:
	    throw new LmrException("Invalid selection for cell variable: " ~ var);
	}
    }

    ref double opIndexAssign(double value, FVCell cell, string var)
    {
	if (var in mSpecies) {
	    cell.fs.gas.massf[mSpecies[var]].re = value;
	    return cell.fs.gas.massf[mSpecies[var]].re;
	}
	if (var in mTModes) {
	    cell.fs.gas.T_modes[mTModes[var]].re = value;
	    return cell.fs.gas.T_modes[mTModes[var]].re;
	}
	if (var in muModes) {
	    cell.fs.gas.u_modes[muModes[var]].re = value;
	    return cell.fs.gas.u_modes[muModes[var]].re;
	}
	if (var in mkModes) {
	    cell.fs.gas.k_modes[mkModes[var]].re = value;
	    return cell.fs.gas.k_modes[mkModes[var]].re;
	}
	if (var in mTurbQuants) {
	    cell.fs.turb[mTurbQuants[var]].re = value;
	    return cell.fs.turb[mTurbQuants[var]].re;
	}

	// For everything else, find by cases
	switch (var) {
	case "pos.x": cell.pos[0].x.re = value; return cell.pos[0].x.re;
	case "pos.y": cell.pos[0].y.re = value; return cell.pos[0].y.re;
	case "pos.z": cell.pos[0].z.re = value; return cell.pos[0].z.re;
	case "vol": cell.volume[0].re = value; return cell.volume[0].re;
	case "rho": cell.fs.gas.rho.re = value; return cell.fs.gas.rho.re;
	case "vel.x": cell.fs.vel.x.re = value; return cell.fs.vel.x.re;
	case "vel.y": cell.fs.vel.y.re = value; return cell.fs.vel.y.re;
	case "vel.z": cell.fs.vel.z.re = value; return cell.fs.vel.z.re;
	case "B.x": cell.fs.B.x.re = value; return cell.fs.B.x.re;
	case "B.y": cell.fs.B.y.re = value; return cell.fs.B.y.re;
	case "B.z": cell.fs.B.z.re = value; return cell.fs.B.z.re;
	case "p": cell.fs.gas.p.re = value; return cell.fs.gas.p.re;
	case "a": cell.fs.gas.a.re = value; return cell.fs.gas.a.re;
	case "mu": cell.fs.gas.mu.re = value; return cell.fs.gas.mu.re;
	case "k": cell.fs.gas.k.re = value; return cell.fs.gas.k.re;
	case "mu_t": cell.fs.mu_t.re = value; return cell.fs.mu_t.re;
	case "k_t": cell.fs.k_t.re = value; return cell.fs.k_t.re;
	case "shock-detector": cell.fs.S.re = value; return cell.fs.S.re;
	case "dt_subcycle": cell.dt_chem.re = value; return cell.dt_chem.re;
	case "dt_local": cell.dt_local.re = value; return cell.dt_local.re;
	case "e": cell.fs.gas.u.re = value; return cell.fs.gas.u.re;
	case "T": cell.fs.gas.T.re = value; return cell.fs.gas.T.re;
	default:
	    throw new LmrException("Invalid selection for cell variable: " ~ var);
	}

    }


private:
    int[string] mSpecies;
    int[string] mTModes;
    int[string] muModes;
    int[string] mkModes;
    int[string] mTurbQuants;
}
