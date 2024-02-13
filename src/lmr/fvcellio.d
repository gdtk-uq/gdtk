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
import std.format : format;
import std.range.primitives : front, popFront;

import nm.number;
import ntypes.complex;

import globalconfig;
import lmrexceptions : LmrException;
import fvcell : FVCell;
import flowstate : FlowState;
import geom : Vector3;

enum FieldVarsType { flow, limiter };

string fieldVarsTypeName(FieldVarsType i)
{
    final switch(i) {
	case FieldVarsType.flow: return "flow";
	case FieldVarsType.limiter: return "limiter";
    }
}

FieldVarsType fieldVarsTypeFromName(string name)
{
    switch (name) {
	case "flow": return FieldVarsType.flow;
	case "limiter": return FieldVarsType.limiter;
	default:
	    string errMsg = "The selection of FieldVarsType is unavailable.\n";
	    errMsg ~= format("You selected '%s'\n", name);
	    errMsg ~= "The available BlockIOTypes are: \n";
	    errMsg ~= "   'flow'\n";
	    errMsg ~= "   'limiter'\n";
	    errMsg ~= "Check the selection or its spelling.\n";
	    throw new Error(errMsg);
    }
}

/**
 * Build the list of variables for the flow field, based on modelling configuation options.
 *
 * The string names for variables should match those used later in this module in:
 *
 *    FVCellFlowIO.opIndex; and
 *    FVCellFlowIO.opIndexAssign.
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
    foreach (isp; 0 .. cfg.gmodel_master.n_species) {
        variables ~= "massf-" ~ cfg.gmodel_master.species_name(isp);
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

abstract class FVCellIO {
public:
    @nogc @property final ref string[] variables() { return mVariables; }
    @nogc @property final FieldVarsType FVT() { return cFVT; }

abstract:
    FVCellIO dup();
    const double opIndex(FVCell cell, string var);
    ref double opIndexAssign(double value, FVCell cell, string var);

private:
    string[] mVariables;
    FieldVarsType cFVT;
}


class FVCellFlowIO : FVCellIO {
public:

    this(string[] variables)
    {
	cFVT = FieldVarsType.flow;
	mVariables = variables.dup;

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

    override
    FVCellFlowIO dup()
    {
	return new FVCellFlowIO(this.mVariables);
    }

    override
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

    override
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



/**
 * Build the list of limiter variables based on modelling configuration.
 *
 * The string names for variables should match those in this module in:
 *
 *   FVCellLimiterIO.opIndex; and
 *   FVCellLimiterIO.opIndexAssign.
 *
 * Authors: RJG
 * Date: 2023-08-13
 */
string[] buildLimiterVariables()
{
    alias cfg = GlobalConfig;
    string[] variables;
    final switch (cfg.thermo_interpolator) {
	case InterpolateOption.pt:
	    variables ~= "p";
	    variables ~= "T";
	    break;
	case InterpolateOption.rhou:
	    variables ~= "rho";
	    variables ~= "e";
	    break;
	case InterpolateOption.rhop:
	    variables ~= "rho";
	    variables ~= "p";
	    break;
	case InterpolateOption.rhot:
	    variables ~= "rho";
	    variables ~= "T";
	    break;
    }
    variables ~= "vel.x";
    variables ~= "vel.y";
    if (cfg.dimensions == 3) variables ~= "vel.z";
    foreach (isp; 0 .. cfg.gmodel_master.n_species) {
	variables ~= "massf-" ~ cfg.gmodel_master.species_name(isp);
    }
    if (cfg.gmodel_master.n_modes > 1) {
        foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
	    if (cfg.thermo_interpolator == InterpolateOption.rhou ||
		cfg.thermo_interpolator == InterpolateOption.rhop ) {
		variables ~= "e-" ~ cfg.gmodel_master.energy_mode_name(imode);
	    }
	    else {
		variables ~= "T-" ~ cfg.gmodel_master.energy_mode_name(imode);
	    }
	}
    }
    if (cfg.turbulence_model_name != "none") {
	foreach (iturb; 0 .. cfg.turb_model.nturb) {
	    variables ~= "tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
	}
    }
    return variables;
}



class FVCellLimiterIO : FVCellIO {
public:

    this(string[] variables)
    {
	cFVT = FieldVarsType.limiter;
	mVariables = variables.dup;

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

    override
    FVCellLimiterIO dup()
    {
	return new FVCellLimiterIO(mVariables);
    }

    override
    const double opIndex(FVCell cell, string var)
    {
	// First handle array-stored values:
	// species, modes and turbulence quantities
	if (var in mSpecies) {
	    return cell.gradients.rho_sPhi[mSpecies[var]].re;
	}
	if (var in mTModes) {
	    return cell.gradients.T_modesPhi[mTModes[var]].re;
	}
	if (var in muModes) {
	    return cell.gradients.u_modesPhi[muModes[var]].re;
	}
	if (var in mTurbQuants) {
	    return cell.gradients.turbPhi[mTurbQuants[var]].re;
	}

	// For everything else, find appropriate case
	switch(var) {
	case "rho": return cell.gradients.rhoPhi.re;
	case "p": return cell.gradients.pPhi.re;
	case "T": return cell.gradients.TPhi.re;
	case "e": return cell.gradients.uPhi.re;
	case "vel.x": return cell.gradients.velxPhi.re;
	case "vel.y": return cell.gradients.velyPhi.re;
	case "vel.z": return cell.gradients.velzPhi.re;
	default:
	    throw new LmrException("Invalid selection for cell gradient variable: " ~ var);
	}
    }

    override
    ref double opIndexAssign(double value, FVCell cell, string var)
    {
	if (var in mSpecies) {
	    cell.gradients.rho_sPhi[mSpecies[var]].re = value;
	    return cell.gradients.rho_sPhi[mSpecies[var]].re;
	}
	if (var in mTModes) {
	    cell.gradients.T_modesPhi[mTModes[var]].re = value;
	    return cell.gradients.T_modesPhi[mTModes[var]].re;
	}
	if (var in muModes) {
	    cell.gradients.u_modesPhi[muModes[var]].re = value;
	    return cell.gradients.u_modesPhi[muModes[var]].re;
	}
	if (var in mTurbQuants) {
	    cell.gradients.turbPhi[mTurbQuants[var]].re = value;
	    return cell.gradients.turbPhi[mTurbQuants[var]].re;
	}

	// For everything else, find appropriate case
	switch(var) {
	case "rho": cell.gradients.rhoPhi.re = value; return cell.gradients.rhoPhi.re;
	case "p": cell.gradients.pPhi.re = value; return cell.gradients.pPhi.re;
	case "T": cell.gradients.TPhi.re = value; return cell.gradients.TPhi.re;
	case "e": cell.gradients.uPhi.re = value; return cell.gradients.uPhi.re;
	case "vel.x": cell.gradients.velxPhi.re = value; return cell.gradients.velxPhi.re;
	case "vel.y": cell.gradients.velyPhi.re = value; return cell.gradients.velyPhi.re;
	case "vel.z": cell.gradients.velzPhi.re = value; return cell.gradients.velzPhi.re;
	default:
	    throw new LmrException("Invalid selection for cell gradient variable: " ~ var);
	}
    }

private:
    int[string] mSpecies;
    int[string] mTModes;
    int[string] muModes;
    int[string] mTurbQuants;
}

FVCellIO createFVCellIO(FieldVarsType fvt, string[] variables)
{
    final switch (fvt) {
	case FieldVarsType.flow: return new FVCellFlowIO(variables);
	case FieldVarsType.limiter: return new FVCellLimiterIO(variables);
    }
}


/**
 * Function moved from Eilmer4. Will need some maintenance.
 *
 * RJG, 2024-02-12
 */
void scan_cell_data_from_fixed_order_string
(string buffer,
 ref Vector3 pos, ref number volume, ref FlowState fs,
 ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
 bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
 bool include_quality, bool MHDflag, bool divergence_cleaning, bool radiation, size_t nturb)
{
    // This function needs to be kept consistent with cell_data_as_string() above.
    auto items = split(buffer);
    version(complex_numbers) {
        // For complex_numbers, we presently set only the real parts.
        // [TODO] Maybe we should read full complex numbers.
        pos.x = Complex!double(items.front); items.popFront();
        pos.y = Complex!double(items.front); items.popFront();
        pos.z = Complex!double(items.front); items.popFront();
        volume = Complex!double(items.front); items.popFront();
        fs.gas.rho = Complex!double(items.front); items.popFront();
        fs.vel.x = Complex!double(items.front); items.popFront();
        fs.vel.y = Complex!double(items.front); items.popFront();
        fs.vel.z = Complex!double(items.front); items.popFront();
        version(MHD) {
            if (MHDflag) {
                fs.B.x = Complex!double(items.front); items.popFront();
                fs.B.y = Complex!double(items.front); items.popFront();
                fs.B.z = Complex!double(items.front); items.popFront();
                fs.divB = Complex!double(items.front); items.popFront();
                if (divergence_cleaning) {
                    fs.psi = Complex!double(items.front); items.popFront();
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fs.gas.quality = Complex!double(items.front); items.popFront();
        } else {
            fs.gas.quality = 1.0;
        }
        fs.gas.p = Complex!double(items.front); items.popFront();
        fs.gas.a = Complex!double(items.front); items.popFront();
        fs.gas.mu = Complex!double(items.front); items.popFront();
        fs.gas.k = Complex!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fs.gas.k_modes[i] = Complex!double(items.front); items.popFront();
            }
        }
        fs.mu_t = Complex!double(items.front); items.popFront();
        fs.k_t = Complex!double(items.front); items.popFront();
        fs.S = Complex!double(items.front); items.popFront();
        if (radiation) {
            Q_rad_org = Complex!double(items.front); items.popFront();
            f_rad_org = Complex!double(items.front); items.popFront();
            Q_rE_rad = Complex!double(items.front); items.popFront();
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(it; 0 .. nturb) {
                fs.turb[it] = Complex!double(items.front); items.popFront();
            }
        }
        version(multi_species_gas) {
            foreach(i; 0 .. fs.gas.massf.length) {
                fs.gas.massf[i] = Complex!double(items.front); items.popFront();
            }
            if (fs.gas.massf.length > 1) {
                dt_chem = to!double(items.front); items.popFront();
            }
        } else {
            items.popFront(); // discard the single-species mass fraction, assumed 1.0
        }
        fs.gas.u = Complex!double(items.front); items.popFront();
        fs.gas.T = Complex!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fs.gas.u_modes[i] = Complex!double(items.front); items.popFront();
                fs.gas.T_modes[i] = Complex!double(items.front); items.popFront();
            }
            if (fs.gas.u_modes.length > 0) {
                dt_therm = to!double(items.front); items.popFront();
            }
        }
        if (with_local_time_stepping) { dt_local = to!double(items.front); items.popFront(); }
    } else {
        // version double_numbers
        pos.x = to!double(items.front); items.popFront();
        pos.y = to!double(items.front); items.popFront();
        pos.z = to!double(items.front); items.popFront();
        volume = to!double(items.front); items.popFront();
        fs.gas.rho = to!double(items.front); items.popFront();
        fs.vel.x = to!double(items.front); items.popFront();
        fs.vel.y = to!double(items.front); items.popFront();
        fs.vel.z = to!double(items.front); items.popFront();
        version(MHD) {
            if (MHDflag) {
                fs.B.x = to!double(items.front); items.popFront();
                fs.B.y = to!double(items.front); items.popFront();
                fs.B.z = to!double(items.front); items.popFront();
                fs.divB = to!double(items.front); items.popFront();
                if (divergence_cleaning) {
                    fs.psi = to!double(items.front); items.popFront();
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fs.gas.quality = to!double(items.front); items.popFront();
        } else {
            fs.gas.quality = 1.0;
        }
        fs.gas.p = to!double(items.front); items.popFront();
        fs.gas.a = to!double(items.front); items.popFront();
        fs.gas.mu = to!double(items.front); items.popFront();
        fs.gas.k = to!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fs.gas.k_modes[i] = to!double(items.front); items.popFront();
            }
        }
        fs.mu_t = to!double(items.front); items.popFront();
        fs.k_t = to!double(items.front); items.popFront();
        fs.S = to!double(items.front); items.popFront();
        if (radiation) {
            Q_rad_org = to!double(items.front); items.popFront();
            f_rad_org = to!double(items.front); items.popFront();
            Q_rE_rad = to!double(items.front); items.popFront();
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(it; 0 .. nturb) {
                fs.turb[it] = to!double(items.front); items.popFront();
            }
        }
        version(multi_species_gas) {
            foreach(i; 0 .. fs.gas.massf.length) {
                fs.gas.massf[i] = to!double(items.front); items.popFront();
            }
            if (fs.gas.massf.length > 1) {
                dt_chem = to!double(items.front); items.popFront();
            }
        } else {
            items.popFront(); // discard single-species mass fraction
        }
        fs.gas.u = to!double(items.front); items.popFront();
        fs.gas.T = to!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fs.gas.u_modes[i] = to!double(items.front); items.popFront();
                fs.gas.T_modes[i] = to!double(items.front); items.popFront();
            }
            if (fs.gas.u_modes.length > 0) {
                dt_therm = to!double(items.front); items.popFront();
            }
        }
        if (with_local_time_stepping) { dt_local = to!double(items.front); items.popFront(); }
    } // end version double_numbers
    version(multi_species_gas) {
        foreach(i; 0 .. fs.gas.massf.length) { fs.gas.rho_s[i] = fs.gas.massf[i] * fs.gas.rho; }
    }
} // end scan_values_from_fixed_order_string()
