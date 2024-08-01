/**
 * A class to access data in a finite-volume cell for use with I/O functions.
 *
 * Author: Rowan J. Gollan
 * Date: 2023-06-20
 */

module fvcellio;

import std.stdio;
import std.algorithm : findSplit, canFind, startsWith;
import std.string : split;
import std.conv : to;
import std.format : format;
import std.range.primitives : front, popFront;

import nm.number;
import ntypes.complex;

import globalconfig;
import lmrexceptions : LmrException;
import lmr.fluidfvcell : FluidFVCell;
import lmr.fvcell : FVCell;
import solidfvcell : SolidFVCell;
import flowstate : FlowState;
import geom : Grid_t, Vector3;

enum FieldVarsType { fluid, solid, limiter, residual, gradient };

string fieldVarsTypeName(FieldVarsType i)
{
    final switch(i) {
    case FieldVarsType.fluid: return "fluid";
    case FieldVarsType.solid: return "solid";
    case FieldVarsType.limiter: return "limiter";
    case FieldVarsType.residual: return "residual";
    case FieldVarsType.gradient: return "gradient";
    }
}

FieldVarsType fieldVarsTypeFromName(string name)
{
    switch (name) {
    case "fluid": return FieldVarsType.fluid;
    case "solid": return FieldVarsType.solid;
    case "limiter": return FieldVarsType.limiter;
    case "residual": return FieldVarsType.residual;
    case "gradient": return FieldVarsType.gradient;
	default:
        string errMsg = "The selection of FieldVarsType is unavailable.\n";
        errMsg ~= format("You selected '%s'\n", name);
        errMsg ~= "The available BlockIOTypes are: \n";
        errMsg ~= "   'fluid'\n";
        errMsg ~= "   'solid'\n";
        errMsg ~= "   'limiter'\n";
        errMsg ~= "   'residual'\n";
        errMsg ~= "   'gradient'\n";
        errMsg ~= "Check the selection or its spelling.\n";
        throw new Error(errMsg);
    }
}

/**
 * Build the list of variables for the flow field, based on modelling configuation options.
 *
 * The string names for variables should match those used later in this module in:
 *
 *    FluidFVCellIO.opIndex; and
 *    FluidFVCellIO.opIndexAssign.
 *
 * Authors: RJG
 */
string[] buildFluidVariables()
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
        variables ~= "B.z";
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
    //
    return variables;
} // end buildFluidVariables()


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
} // end class FVCellIO


class FluidFVCellIO : FVCellIO {
public:

    this()
    {
        this(buildFluidVariables());
    }

    this(string[] variables)
    {
        cFVT = FieldVarsType.fluid;
	    mVariables = variables.dup;

        alias cfg = GlobalConfig;
        foreach (var; variables) {
            if (var.startsWith("massf-")) {
                auto spName = var["massf-".length..$];
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
    } // end constructor from an array of strings

    override
    FluidFVCellIO dup()
    {
        return new FluidFVCellIO(mVariables);
    }

    override
    const double opIndex(FVCell cell, string var)
    {
        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
    	// First handle array-stored values:
    	// species, modes and turbulence quantities
        if (var in mSpecies) {
            return fcell.fs.gas.massf[mSpecies[var]].re;
    	}
        if (var in mTModes) {
            return fcell.fs.gas.T_modes[mTModes[var]].re;
        }
        if (var in muModes) {
            return fcell.fs.gas.u_modes[muModes[var]].re;
        }
        if (var in mkModes) {
            return fcell.fs.gas.k_modes[mkModes[var]].re;
        }
        if (var in mTurbQuants) {
            return fcell.fs.turb[mTurbQuants[var]].re;
        }
        // For everything else, find by cases
    	switch (var) {
        case "pos.x": return fcell.pos[0].x.re;
    	case "pos.y": return fcell.pos[0].y.re;
        case "pos.z": return fcell.pos[0].z.re;
    	case "vol": return fcell.volume[0].re;
        case "rho": return fcell.fs.gas.rho.re;
    	case "vel.x": return fcell.fs.vel.x.re;
        case "vel.y": return fcell.fs.vel.y.re;
    	case "vel.z": return fcell.fs.vel.z.re;
    	case "B.x": return fcell.fs.B.x.re;
        case "B.y": return fcell.fs.B.y.re;
    	case "B.z": return fcell.fs.B.z.re;
    	case "p": return fcell.fs.gas.p.re;
    	case "a": return fcell.fs.gas.a.re;
        case "mu": return fcell.fs.gas.mu.re;
    	case "k": return fcell.fs.gas.k.re;
        case "mu_t": return fcell.fs.mu_t.re;
    	case "k_t": return fcell.fs.k_t.re;
        case "shock-detector": return fcell.fs.S.re;
    	case "dt_subcycle": return fcell.dt_chem.re;
        case "dt_local": return fcell.dt_local.re;
        case "e": return fcell.fs.gas.u.re;
        case "T": return fcell.fs.gas.T.re;
        default:
            throw new LmrException("Invalid selection for cell variable: " ~ var);
        }
    } // end opIndex()

    override
    ref double opIndexAssign(double value, FVCell cell, string var)
    {
        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
        if (var in mSpecies) {
            fcell.fs.gas.massf[mSpecies[var]].re = value;
            return fcell.fs.gas.massf[mSpecies[var]].re;
        }
        if (var in mTModes) {
            fcell.fs.gas.T_modes[mTModes[var]].re = value;
            return fcell.fs.gas.T_modes[mTModes[var]].re;
        }
        if (var in muModes) {
            fcell.fs.gas.u_modes[muModes[var]].re = value;
            return fcell.fs.gas.u_modes[muModes[var]].re;
        }
        if (var in mkModes) {
            fcell.fs.gas.k_modes[mkModes[var]].re = value;
            return fcell.fs.gas.k_modes[mkModes[var]].re;
        }
        if (var in mTurbQuants) {
            fcell.fs.turb[mTurbQuants[var]].re = value;
            return fcell.fs.turb[mTurbQuants[var]].re;
        }
        // For everything else, find by cases
        switch (var) {
        case "pos.x": fcell.pos[0].x.re = value; return fcell.pos[0].x.re;
        case "pos.y": fcell.pos[0].y.re = value; return fcell.pos[0].y.re;
        case "pos.z": fcell.pos[0].z.re = value; return fcell.pos[0].z.re;
        case "vol": fcell.volume[0].re = value; return fcell.volume[0].re;
        case "rho": fcell.fs.gas.rho.re = value; return fcell.fs.gas.rho.re;
        case "vel.x": fcell.fs.vel.x.re = value; return fcell.fs.vel.x.re;
        case "vel.y": fcell.fs.vel.y.re = value; return fcell.fs.vel.y.re;
        case "vel.z": fcell.fs.vel.z.re = value; return fcell.fs.vel.z.re;
        case "B.x": fcell.fs.B.x.re = value; return fcell.fs.B.x.re;
        case "B.y": fcell.fs.B.y.re = value; return fcell.fs.B.y.re;
        case "B.z": fcell.fs.B.z.re = value; return fcell.fs.B.z.re;
        case "p": fcell.fs.gas.p.re = value; return fcell.fs.gas.p.re;
        case "a": fcell.fs.gas.a.re = value; return fcell.fs.gas.a.re;
        case "mu": fcell.fs.gas.mu.re = value; return fcell.fs.gas.mu.re;
        case "k": fcell.fs.gas.k.re = value; return fcell.fs.gas.k.re;
        case "mu_t": fcell.fs.mu_t.re = value; return fcell.fs.mu_t.re;
        case "k_t": fcell.fs.k_t.re = value; return fcell.fs.k_t.re;
        case "shock-detector": fcell.fs.S.re = value; return fcell.fs.S.re;
        case "dt_subcycle": fcell.dt_chem.re = value; return fcell.dt_chem.re;
        case "dt_local": fcell.dt_local.re = value; return fcell.dt_local.re;
        case "e": fcell.fs.gas.u.re = value; return fcell.fs.gas.u.re;
        case "T": fcell.fs.gas.T.re = value; return fcell.fs.gas.T.re;
        default:
            throw new LmrException("Invalid selection for cell variable: " ~ var);
        }
    } // end opIndexAssign()

private:
    int[string] mSpecies;
    int[string] mTModes;
    int[string] muModes;
    int[string] mkModes;
    int[string] mTurbQuants;
} // end class FluidFVCellIO

/**
 * Build the list of solid variables.
 *
 * Authors: RJG
 * Date: 2024-02-25
 */
string[] buildSolidVariables()
{

    alias cfg = GlobalConfig;
    string[] variables;
    variables ~= "pos.x";
    variables ~= "pos.y";
    if (cfg.dimensions == 3) variables ~= "pos.z";
    variables ~= "vol";
    variables ~= "e";
    variables ~= "T";
    variables ~= "rho";
    variables ~= "Cp";
    variables ~= "k";
    return variables;
}

class SolidFVCellIO : FVCellIO {
public:

    this()
    {
        this(buildSolidVariables());
    }

    this(string[] variables)
    {
        cFVT = FieldVarsType.solid;
        mVariables = variables.dup;
    }

    override
    SolidFVCellIO dup()
    {
        return new SolidFVCellIO(mVariables);
    }

    override
    const double opIndex(FVCell cell, string var)
    {
        auto scell = cast(SolidFVCell) cell;
        if (scell is null) {
            throw new LmrException("Invalid cast to SolidFVCell.");
        }
        switch (var) {
        case "pos.x": return scell.pos.x.re;
        case "pos.y": return scell.pos.y.re;
        case "pos.z": return scell.pos.z.re;
        case "vol": return scell.volume.re;
        case "e": return scell.e[0].re;
        case "T": return scell.T.re;
        case "rho": return scell.ss.rho.re;
        case "Cp": return scell.ss.Cp.re;
        case "k": return scell.ss.k.re;
        default:
            throw new LmrException("Invalid selection for solid cell variable: " ~ var);
        }
    } // end opIndex()

    override
    ref double opIndexAssign(double value, FVCell cell, string var)
    {
        auto scell = cast(SolidFVCell) cell;
        if (scell is null) {
            throw new LmrException("Invalid cast to SolidFVCell.");
        }
        switch (var) {
        case "pos.x": scell.pos.x.re = value; return scell.pos.x.re;
        case "pos.y": scell.pos.y.re = value; return scell.pos.y.re;
        case "pos.z": scell.pos.z.re = value; return scell.pos.z.re;
        case "vol": scell.volume.re = value; return scell.volume.re;
        case "e": scell.e[0].re = value; return scell.e[0].re;
        case "T": scell.T.re = value; return scell.T.re;
        case "rho": scell.ss.rho.re = value; return scell.ss.rho.re;
        case "Cp": scell.ss.Cp.re = value; return scell.ss.Cp.re;
        case "k": scell.ss.k.re = value; return scell.ss.k.re;
        default:
            throw new LmrException("Invalid selection for cell variable: " ~ var);
        }
    } // end opIndexAssign()

} // end class SolidFVCellIO


/**
 * Build the list of limiter variables based on modelling configuration.
 *
 * The string names for variables should match those in this module in:
 *
 *   FluidFVCellLimiterIO.opIndex; and
 *   FluidFVCellLimiterIO.opIndexAssign.
 *
 * Authors: RJG
 * Date: 2023-08-13
 */
string[] buildLimiterVariables(Grid_t grid_type)
{
    alias cfg = GlobalConfig;
    bool add_limiter_values = (grid_type == Grid_t.unstructured_grid) && (cfg.interpolation_order > 1);
    string[] variables;

    if (add_limiter_values) {
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
        } // end switch
        variables ~= "vel.x";
        variables ~= "vel.y";
        if (cfg.dimensions == 3) variables ~= "vel.z";
        foreach (isp; 0 .. cfg.gmodel_master.n_species) {
            variables ~= "massf-" ~ cfg.gmodel_master.species_name(isp);
        }
        foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
            if (cfg.thermo_interpolator == InterpolateOption.rhou ||
                cfg.thermo_interpolator == InterpolateOption.rhop ) {
                variables ~= "e-" ~ cfg.gmodel_master.energy_mode_name(imode);
            } else {
                variables ~= "T-" ~ cfg.gmodel_master.energy_mode_name(imode);
            }
        }
        if (cfg.turbulence_model_name != "none") {
            foreach (iturb; 0 .. cfg.turb_model.nturb) {
                variables ~= "tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
            }
        }
    }

    if (variables.length == 0 && cfg.is_master_task) {
        string errMsg = "ERROR: NO LIMITER VARIABLES FOUND.\n";
        errMsg ~= "** You have opted for cell-centered limiter values to be written to disk, however, no limiter variables have been found. \n";
        errMsg ~= "** Note this option is only available for unstructured grid simulations with interpolation_order > 1. \n";
        throw new Error(errMsg);
    }

    return variables;
} // end buildLimiterVariables()



class FluidFVCellLimiterIO : FVCellIO {
public:

    this(Grid_t grid_type)
    {
        this(buildLimiterVariables(grid_type));
    }

    this(string[] variables)
    {
	cFVT = FieldVarsType.limiter;
	mVariables = variables.dup;

	alias cfg = GlobalConfig;
	foreach (var; variables) {
	    if (var.startsWith("massf-")) {
                auto spName = var["massf-".length..$];
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
    } // end constructor

    override
    FluidFVCellLimiterIO dup()
    {
    	return new FluidFVCellLimiterIO(mVariables);
    }

    override
    const double opIndex(FVCell cell, string var)
    {
        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
    	// First handle array-stored values:
    	// species, modes and turbulence quantities
        if (var in mSpecies) {
            return fcell.gradients.rho_sPhi[mSpecies[var]].re;
        }
    	if (var in mTModes) {
    	    return fcell.gradients.T_modesPhi[mTModes[var]].re;
    	}
        if (var in muModes) {
            return fcell.gradients.u_modesPhi[muModes[var]].re;
    	}
        if (var in mTurbQuants) {
            return fcell.gradients.turbPhi[mTurbQuants[var]].re;
        }
    	// For everything else, find appropriate case
        switch(var) {
    	case "rho": return fcell.gradients.rhoPhi.re;
    	case "p": return fcell.gradients.pPhi.re;
    	case "T": return fcell.gradients.TPhi.re;
    	case "e": return fcell.gradients.uPhi.re;
    	case "vel.x": return fcell.gradients.velxPhi.re;
    	case "vel.y": return fcell.gradients.velyPhi.re;
    	case "vel.z": return fcell.gradients.velzPhi.re;
        default:
            throw new LmrException("Invalid selection for cell gradient variable: " ~ var);
        }
    } // end opIndex()

    override
    ref double opIndexAssign(double value, FVCell cell, string var)
    {
        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
    	if (var in mSpecies) {
            fcell.gradients.rho_sPhi[mSpecies[var]].re = value;
            return fcell.gradients.rho_sPhi[mSpecies[var]].re;
    	}
        if (var in mTModes) {
            fcell.gradients.T_modesPhi[mTModes[var]].re = value;
            return fcell.gradients.T_modesPhi[mTModes[var]].re;
    	}
        if (var in muModes) {
            fcell.gradients.u_modesPhi[muModes[var]].re = value;
            return fcell.gradients.u_modesPhi[muModes[var]].re;
        }
    	if (var in mTurbQuants) {
            fcell.gradients.turbPhi[mTurbQuants[var]].re = value;
            return fcell.gradients.turbPhi[mTurbQuants[var]].re;
        }
        // For everything else, find appropriate case
        switch(var) {
    	case "rho": fcell.gradients.rhoPhi.re = value; return fcell.gradients.rhoPhi.re;
    	case "p": fcell.gradients.pPhi.re = value; return fcell.gradients.pPhi.re;
    	case "T": fcell.gradients.TPhi.re = value; return fcell.gradients.TPhi.re;
    	case "e": fcell.gradients.uPhi.re = value; return fcell.gradients.uPhi.re;
    	case "vel.x": fcell.gradients.velxPhi.re = value; return fcell.gradients.velxPhi.re;
    	case "vel.y": fcell.gradients.velyPhi.re = value; return fcell.gradients.velyPhi.re;
    	case "vel.z": fcell.gradients.velzPhi.re = value; return fcell.gradients.velzPhi.re;
    	default:
            throw new LmrException("Invalid selection for cell gradient variable: " ~ var);
    	}
    } // end opIndexAssign()

private:
    int[string] mSpecies;
    int[string] mTModes;
    int[string] muModes;
    int[string] mTurbQuants;
} // end class FluidFVCellLimiterIO

/**
 * Build the list of residual variables based on modelling configuration.
 *
 * The string names for variables should match those in this module in:
 *
 *   FluidFVCellResidualIO.opIndex; and
 *   FluidFVCellResidualIO.opIndexAssign.
 *
 * Authors: KAD and RJG
 * Date: 2024-03-07
 */
string[] buildResidualVariables()
{
    alias cfg = GlobalConfig;
    string[] variables;
    if (cfg.gmodel_master.n_species > 1) {
        // we do not have a (total) mass residual for multi-species simulations
        foreach (isp; 0 .. cfg.gmodel_master.n_species) {
            variables ~= cfg.gmodel_master.species_name(isp);
        }
    } else {
        variables ~= "mass";
    }
    variables ~= "x-momentum";
    variables ~= "y-momentum";
    if (cfg.dimensions == 3) variables ~= "z-momentum";
    variables ~= "total-energy";
    foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
        variables ~= cfg.gmodel_master.energy_mode_name(imode);
    }
    if (cfg.turbulence_model_name != "none") {
	foreach (iturb; 0 .. cfg.turb_model.nturb) {
	    variables ~= cfg.turb_model.primitive_variable_name(iturb);
	}
    }
    return variables;
} // end buildResidualVariables()



class FluidFVCellResidualIO : FVCellIO {
public:

    this()
    {
        this(buildResidualVariables());
    }

    this(string[] variables)
    {
	cFVT = FieldVarsType.residual;
	mVariables = variables.dup;

	alias cfg = GlobalConfig;

        string[] spList;
        foreach (i; 0 .. cfg.gmodel_master.n_species) {
            spList ~= cfg.gmodel_master.species_name(i);
        }
        string[] umList;
        foreach (i; 0 .. cfg.gmodel_master.n_modes) {
            umList ~= cfg.gmodel_master.energy_mode_name(i);
        }
        string[] tbList;
        foreach (i; 0 .. cfg.turb_model.nturb) {
            tbList ~= cfg.turb_model.primitive_variable_name(i);
        }

	foreach (var; variables) {
            if (spList.canFind(var)) {
                mSpecies[var] = cfg.gmodel_master.species_index(var);
            }
            if (umList.canFind(var)) {
                muModes[var] = cfg.gmodel_master.energy_mode_index(var);
            }
            if (tbList.canFind(var)) {
                mTurbQuants[var] = cfg.turb_model.primitive_variable_index(var);
            }
	}
    } // end constructor from string array

    override
    FluidFVCellResidualIO dup()
    {
        return new FluidFVCellResidualIO(mVariables);
    }

    override
    const double opIndex(FVCell cell, string var)
    {
	alias cfg = GlobalConfig;

        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
        // First handle array-stored values:
        // species, modes and turbulence quantities
        if (var in mSpecies) {
            return fcell.dUdt[0][cfg.cqi.species+mSpecies[var]].re;
        }
        if (var in muModes) {
            return fcell.dUdt[0][cfg.cqi.modes+muModes[var]].re;
        }
        if (var in mTurbQuants) {
            return fcell.dUdt[0][cfg.cqi.turb+mTurbQuants[var]].re;
        }

        // For everything else, find appropriate case
        switch(var) {
        case "mass": return fcell.dUdt[0][cfg.cqi.mass].re;
        case "x-momentum": return fcell.dUdt[0][cfg.cqi.xMom].re;
        case "y-momentum": return fcell.dUdt[0][cfg.cqi.yMom].re;
        case "z-momentum": return fcell.dUdt[0][cfg.cqi.zMom].re;
        case "total-energy": return fcell.dUdt[0][cfg.cqi.totEnergy].re;
        default:
            throw new LmrException("Invalid selection for cell dUdt entry: " ~ var);
	}
    } // end opIndex()

    override
    ref double opIndexAssign(double value, FVCell cell, string var)
    {
	alias cfg = GlobalConfig;

        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
        if (var in mSpecies) {
            fcell.dUdt[0][cfg.cqi.species+mSpecies[var]].re = value;
            return fcell.dUdt[0][cfg.cqi.species+mSpecies[var]].re;
        }
        if (var in muModes) {
            fcell.dUdt[0][cfg.cqi.modes+muModes[var]].re = value;
            return fcell.dUdt[0][cfg.cqi.modes+muModes[var]].re;
        }
        if (var in mTurbQuants) {
            fcell.dUdt[0][cfg.cqi.turb+mTurbQuants[var]].re = value;
            return fcell.dUdt[0][cfg.cqi.turb+mTurbQuants[var]].re;
        }

        // For everything else, find appropriate case
        switch(var) {
        case "mass": fcell.dUdt[0][cfg.cqi.mass].re = value; return fcell.dUdt[0][cfg.cqi.mass].re;
        case "x-momentum": fcell.dUdt[0][cfg.cqi.xMom].re = value; return fcell.dUdt[0][cfg.cqi.xMom].re;
        case "y-momentum": fcell.dUdt[0][cfg.cqi.yMom].re = value; return fcell.dUdt[0][cfg.cqi.yMom].re;
        case "z-momentum": fcell.dUdt[0][cfg.cqi.zMom].re = value; return fcell.dUdt[0][cfg.cqi.zMom].re;
        case "total-energy": fcell.dUdt[0][cfg.cqi.totEnergy].re = value; return fcell.dUdt[0][cfg.cqi.totEnergy].re;
        default:
            throw new LmrException("Invalid selection for cell dUdt entry: " ~ var);
        }
    } // end opIndexAssign()

private:
    int[string] mSpecies;
    int[string] muModes;
    int[string] mTurbQuants;
} // end class FluidFVCellResidualIO


/**
 * Build the list of gradient variables based on modelling configuration.
 *
 * The string names for variables should match those in this module in:
 *
 *   FluidFVCellGradientIO.opIndex; and
 *   FluidFVCellGradientIO.opIndexAssign.
 *
 * Authors: KAD and RJG
 * Date: 2024-08-01
 */
string[] buildGradientVariables(Grid_t grid_type)
{
    alias cfg = GlobalConfig;
    bool add_convective_gradients = (grid_type == Grid_t.unstructured_grid) && (cfg.interpolation_order > 1);
    bool add_viscous_gradients = cfg.viscous &&
        (cfg.spatial_deriv_calc == SpatialDerivCalc.least_squares) &&
        (cfg.spatial_deriv_locn == SpatialDerivLocn.cells);
    bool add_z_component = (cfg.dimensions == 3);
    string[] variables;

    // collect the convective gradients
    if (add_convective_gradients) {
        final switch (cfg.thermo_interpolator) {
        case InterpolateOption.pt:
            variables ~= "conv_x_p";
            variables ~= "conv_y_p";
            variables ~= "conv_x_T";
            variables ~= "conv_y_T";
            if (add_z_component) {
                variables ~= "conv_z_p";
                variables ~= "conv_z_T";
            }
            break;
        case InterpolateOption.rhou:
            variables ~= "conv_x_rho";
            variables ~= "conv_y_rho";
            variables ~= "conv_x_e";
            variables ~= "conv_y_e";
            if (add_z_component) {
                variables ~= "conv_z_rho";
                variables ~= "conv_z_e";
            }
            break;
        case InterpolateOption.rhop:
            variables ~= "conv_x_rho";
            variables ~= "conv_y_rho";
            variables ~= "conv_x_p";
            variables ~= "conv_y_p";
            if (add_z_component) {
                variables ~= "conv_z_rho";
                variables ~= "conv_z_p";
            }
            break;
        case InterpolateOption.rhot:
            variables ~= "conv_x_rho";
            variables ~= "conv_y_rho";
            variables ~= "conv_x_T";
            variables ~= "conv_y_T";
            if (add_z_component) {
                variables ~= "conv_z_rho";
                variables ~= "conv_z_T";
            }
            break;
        } // end switch
        if (add_convective_gradients) {
            variables ~= "conv_x_vel.x";
            variables ~= "conv_y_vel.x";
            variables ~= "conv_x_vel.y";
            variables ~= "conv_y_vel.y";
        }
        if (add_z_component) {
            variables ~= "conv_x_vel.z";
            variables ~= "conv_y_vel.z";
            variables ~= "conv_z_vel.x";
            variables ~= "conv_z_vel.y";
            variables ~= "conv_z_vel.z";
        }
        foreach (isp; 0 .. cfg.gmodel_master.n_species) {
            variables ~= "conv_x_massf-" ~ cfg.gmodel_master.species_name(isp);
            variables ~= "conv_y_massf-" ~ cfg.gmodel_master.species_name(isp);
            if (add_z_component) {
                variables ~= "conv_z_massf-" ~ cfg.gmodel_master.species_name(isp);
            }
        }
        foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
            if (cfg.thermo_interpolator == InterpolateOption.rhou ||
                cfg.thermo_interpolator == InterpolateOption.rhop ) {
                variables ~= "conv_x_e-" ~ cfg.gmodel_master.energy_mode_name(imode);
                variables ~= "conv_y_e-" ~ cfg.gmodel_master.energy_mode_name(imode);
                if (add_z_component) {
                    variables ~= "conv_z_e-" ~ cfg.gmodel_master.energy_mode_name(imode);
                }
            } else {
                variables ~= "conv_x_T-" ~ cfg.gmodel_master.energy_mode_name(imode);
                variables ~= "conv_y_T-" ~ cfg.gmodel_master.energy_mode_name(imode);
                if (add_z_component) {
                    variables ~= "conv_z_T-" ~ cfg.gmodel_master.energy_mode_name(imode);
                }
            }
        }
        if (cfg.turbulence_model_name != "none") {
            foreach (iturb; 0 .. cfg.turb_model.nturb) {
                variables ~= "conv_x_tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
                variables ~= "conv_y_tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
                if (add_z_component) {
                    variables ~= "conv_z_tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
                }
            }
        }
    }

    // collect the viscous gradients
    if (add_viscous_gradients) {
        variables ~= "visc_x_vel.x";
        variables ~= "visc_y_vel.x";
        variables ~= "visc_x_vel.y";
        variables ~= "visc_y_vel.y";
        if (add_z_component) {
            variables ~= "visc_x_vel.z";
            variables ~= "visc_y_vel.z";
            variables ~= "visc_z_vel.x";
            variables ~= "visc_z_vel.y";
            variables ~= "visc_z_vel.z";
        }
        variables ~= "visc_x_T";
        variables ~= "visc_y_T";
        if (add_z_component) {
            variables ~= "visc_z_T";
        }
        foreach (isp; 0 .. cfg.gmodel_master.n_species) {
            variables ~= "visc_x_massf-" ~ cfg.gmodel_master.species_name(isp);
            variables ~= "visc_y_massf-" ~ cfg.gmodel_master.species_name(isp);
            if (add_z_component) {
                variables ~= "visc_z_massf-" ~ cfg.gmodel_master.species_name(isp);
            }
        }
        foreach (imode; 0 .. cfg.gmodel_master.n_modes) {
            variables ~= "visc_x_T-" ~ cfg.gmodel_master.energy_mode_name(imode);
            variables ~= "visc_y_T-" ~ cfg.gmodel_master.energy_mode_name(imode);
            if (add_z_component) {
                variables ~= "visc_z_T-" ~ cfg.gmodel_master.energy_mode_name(imode);
            }
        }
        if (cfg.turbulence_model_name != "none") {
            foreach (iturb; 0 .. cfg.turb_model.nturb) {
                variables ~= "visc_x_tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
                variables ~= "visc_y_tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
                if (add_z_component) {
                    variables ~= "visc_z_tq-" ~ cfg.turb_model.primitive_variable_name(iturb);
                }
            }
        }
    }

    if (variables.length == 0 && cfg.is_master_task) {
        string errMsg = "ERROR: NO GRADIENT VARIABLES FOUND.\n";
        errMsg ~= "** You have opted for cell-centered gradient values to be written to disk, however, no gradient variables have been found. \n";
        errMsg ~= "** Note this option is only available for either: \n";
        errMsg ~= "**     1. inviscid unstructured grid simulations with interpolation_order > 1. \n";
        errMsg ~= "** or \n";
        errMsg ~= "**     2. viscous simulations using spatial_deriv_calc set to least_squares and spatial_deriv_locn set to cells. \n";
        throw new Error(errMsg);
    }

    return variables;
} // end buildGradientVariables()

class FluidFVCellGradientIO : FVCellIO {
public:

    this(Grid_t grid_type)
    {
        this(buildGradientVariables(grid_type));
    }

    this(string[] variables)
    {
        cFVT = FieldVarsType.gradient;
        mVariables = variables.dup;

        alias cfg = GlobalConfig;
        foreach (var; variables) {
            if (var.canFind("massf-")) {
                auto spName = var.findSplit("massf-")[2];
                auto spIdx = cfg.gmodel_master.species_index(spName);
                if (spIdx != -1)
                    mSpecies[var.idup] = spIdx;
                else
                    throw new LmrException("Invalid species name: " ~ spName);
            }
            if (var.canFind("T-")) {
                auto modeName = var.split("-")[1];
                auto modeIdx = cfg.gmodel_master.energy_mode_index(modeName);
                if (modeIdx != -1)
                    mTModes[var.idup] = modeIdx;
                else
                    throw new LmrException("Invalid mode name: " ~ modeName);
            }
            if (var.canFind("e-")) {
                auto modeName = var.split("-")[1];
                auto modeIdx = cfg.gmodel_master.energy_mode_index(modeName);
                if (modeIdx != -1)
                    muModes[var.idup] = modeIdx;
                else
                    throw new LmrException("Invalid mode name: " ~ modeName);
            }
            if (var.canFind("tq-")) {
                auto tqName = var.split("-")[1];
                auto tqIdx = cfg.turb_model.primitive_variable_index(tqName);
                if (tqIdx != -1)
                    mTurbQuants[var.idup] = tqIdx;
                else
                    throw new LmrException("Invalid turbulence quantity name: " ~ tqName);
            }
        }
    } // end constructor

    override
    FluidFVCellGradientIO dup()
    {
        return new FluidFVCellGradientIO(mVariables);
    }

    override
    const double opIndex(FVCell cell, string var)
    {
        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }
        // First handle array-stored values:
        // species, modes and turbulence quantities
        if (var in mSpecies) {
            if (var.canFind("conv_x")) { return fcell.gradients.rho_s[mSpecies[var]][0].re; }
            if (var.canFind("conv_y")) { return fcell.gradients.rho_s[mSpecies[var]][1].re; }
            if (var.canFind("conv_z")) { return fcell.gradients.rho_s[mSpecies[var]][2].re; }
            if (var.canFind("visc_x")) { return fcell.grad.massf[mSpecies[var]][0].re; }
            if (var.canFind("visc_y")) { return fcell.grad.massf[mSpecies[var]][1].re; }
            if (var.canFind("visc_z")) { return fcell.grad.massf[mSpecies[var]][2].re; }
        }
        if (var in mTModes) {
            if (var.canFind("conv_x")) { return fcell.gradients.T_modes[mTModes[var]][0].re; }
            if (var.canFind("conv_y")) { return fcell.gradients.T_modes[mTModes[var]][1].re; }
            if (var.canFind("conv_z")) { return fcell.gradients.T_modes[mTModes[var]][2].re; }
            if (var.canFind("visc_x")) { return fcell.grad.T_modes[mTModes[var]][0].re; }
            if (var.canFind("visc_y")) { return fcell.grad.T_modes[mTModes[var]][1].re; }
            if (var.canFind("visc_z")) { return fcell.grad.T_modes[mTModes[var]][2].re; }
        }
        if (var in muModes) {
            if (var.canFind("conv_x")) { return fcell.gradients.u_modes[muModes[var]][0].re; }
            if (var.canFind("conv_y")) { return fcell.gradients.u_modes[muModes[var]][1].re; }
            if (var.canFind("conv_z")) { return fcell.gradients.u_modes[muModes[var]][2].re; }
        }
        if (var in mTurbQuants) {
            if (var.canFind("conv_x")) { return fcell.gradients.turb[mTurbQuants[var]][0].re; }
            if (var.canFind("conv_y")) { return fcell.gradients.turb[mTurbQuants[var]][1].re; }
            if (var.canFind("conv_z")) { return fcell.gradients.turb[mTurbQuants[var]][2].re; }
            if (var.canFind("visc_x")) { return fcell.grad.turb[mTurbQuants[var]][0].re; }
            if (var.canFind("visc_y")) { return fcell.grad.turb[mTurbQuants[var]][1].re; }
            if (var.canFind("visc_z")) { return fcell.grad.turb[mTurbQuants[var]][2].re; }
        }
        // For everything else, find appropriate case
        switch(var) {
        case "conv_x_rho": return fcell.gradients.rho[0].re;
        case "conv_y_rho": return fcell.gradients.rho[1].re;
        case "conv_z_rho": return fcell.gradients.rho[2].re;
        case "conv_x_p": return fcell.gradients.p[0].re;
        case "conv_y_p": return fcell.gradients.p[1].re;
        case "conv_z_p": return fcell.gradients.p[2].re;
        case "conv_x_e": return fcell.gradients.u[0].re;
        case "conv_y_e": return fcell.gradients.u[1].re;
        case "conv_z_e": return fcell.gradients.u[2].re;
        case "conv_x_T": return fcell.gradients.T[0].re;
        case "conv_y_T": return fcell.gradients.T[1].re;
        case "conv_z_T": return fcell.gradients.T[2].re;
        case "visc_x_T": return fcell.grad.T[0].re;
        case "visc_y_T": return fcell.grad.T[1].re;
        case "visc_z_T": return fcell.grad.T[2].re;
        case "conv_x_vel.x": return fcell.gradients.velx[0].re;
        case "conv_y_vel.x": return fcell.gradients.velx[1].re;
        case "conv_z_vel.x": return fcell.gradients.velx[2].re;
        case "visc_x_vel.x": return fcell.grad.vel[0][0].re;
        case "visc_y_vel.x": return fcell.grad.vel[0][1].re;
        case "visc_z_vel.x": return fcell.grad.vel[0][2].re;
        case "conv_x_vel.y": return fcell.gradients.vely[0].re;
        case "conv_y_vel.y": return fcell.gradients.vely[1].re;
        case "conv_z_vel.y": return fcell.gradients.vely[2].re;
        case "visc_x_vel.y": return fcell.grad.vel[1][0].re;
        case "visc_y_vel.y": return fcell.grad.vel[1][1].re;
        case "visc_z_vel.y": return fcell.grad.vel[1][2].re;
        case "conv_x_vel.z": return fcell.gradients.velz[0].re;
        case "conv_y_vel.z": return fcell.gradients.velz[1].re;
        case "conv_z_vel.z": return fcell.gradients.velz[2].re;
        case "visc_x_vel.z": return fcell.grad.vel[2][0].re;
        case "visc_y_vel.z": return fcell.grad.vel[2][1].re;
        case "visc_z_vel.z": return fcell.grad.vel[2][2].re;
        default:
            throw new LmrException("Invalid selection for cell gradient variable: " ~ var);
        }
    } // end opIndex()

    override
    ref double opIndexAssign(double value, FVCell cell, string var)
    {
        auto fcell = cast(FluidFVCell) cell;
        if (fcell is null) {
            throw new LmrException("Invalid cast to FluidFVCell.");
        }

        if (var in mSpecies) {
            if (var.canFind("conv_x")) {
                fcell.gradients.rho_s[mSpecies[var]][0].re = value;
                return fcell.gradients.rho_s[mSpecies[var]][0].re;
            }
            if (var.canFind("conv_y")) {
                fcell.gradients.rho_s[mSpecies[var]][1].re = value;
                return fcell.gradients.rho_s[mSpecies[var]][1].re;
            }
            if (var.canFind("conv_z")) {
                fcell.gradients.rho_s[mSpecies[var]][2].re = value;
                return fcell.gradients.rho_s[mSpecies[var]][2].re;
            }
            if (var.canFind("visc_x")) {
                fcell.grad.massf[mSpecies[var]][0].re = value;
                return fcell.grad.massf[mSpecies[var]][0].re;
            }
            if (var.canFind("visc_y")) {
                fcell.grad.massf[mSpecies[var]][1].re = value;
                return fcell.grad.massf[mSpecies[var]][1].re;
            }
            if (var.canFind("visc_z")) {
                fcell.grad.massf[mSpecies[var]][2].re = value;
                return fcell.grad.massf[mSpecies[var]][2].re;
            }
        }
        if (var in mTModes) {
            if (var.canFind("conv_x")) {
                fcell.gradients.T_modes[mTModes[var]][0].re = value;
                return fcell.gradients.T_modes[mTModes[var]][0].re;
            }
            if (var.canFind("conv_y")) {
                fcell.gradients.T_modes[mTModes[var]][1].re = value;
                return fcell.gradients.T_modes[mTModes[var]][1].re;
            }
            if (var.canFind("conv_z")) {
                fcell.gradients.T_modes[mTModes[var]][2].re = value;
                return fcell.gradients.T_modes[mTModes[var]][2].re;
            }
            if (var.canFind("visc_x")) {
                fcell.grad.T_modes[mTModes[var]][0].re = value;
                return fcell.grad.T_modes[mTModes[var]][0].re;
            }
            if (var.canFind("visc_y")) {
                fcell.grad.T_modes[mTModes[var]][1].re = value;
                return fcell.grad.T_modes[mTModes[var]][1].re;
            }
            if (var.canFind("visc_z")) {
                fcell.grad.T_modes[mTModes[var]][2].re = value;
                return fcell.grad.T_modes[mTModes[var]][2].re;
            }
        }
        if (var in muModes) {
            if (var.canFind("conv_x")) {
                fcell.gradients.u_modes[muModes[var]][0].re = value;
                return fcell.gradients.u_modes[muModes[var]][0].re;
            }
            if (var.canFind("conv_y")) {
                fcell.gradients.u_modes[muModes[var]][1].re = value;
                return fcell.gradients.u_modes[muModes[var]][1].re;
            }
            if (var.canFind("conv_z")) {
                fcell.gradients.u_modes[muModes[var]][2].re = value;
                return fcell.gradients.u_modes[muModes[var]][2].re;
            }
        }
        if (var in mTurbQuants) {
            if (var.canFind("conv_x")) {
                fcell.gradients.turb[mTurbQuants[var]][0].re = value;
                return fcell.gradients.turb[mTurbQuants[var]][0].re;
            }
            if (var.canFind("conv_y")) {
                fcell.gradients.turb[mTurbQuants[var]][1].re = value;
                return fcell.gradients.turb[mTurbQuants[var]][1].re;
            }
            if (var.canFind("conv_z")) {
                fcell.gradients.turb[mTurbQuants[var]][2].re = value;
                return fcell.gradients.turb[mTurbQuants[var]][2].re;
            }
            if (var.canFind("visc_x")) {
                fcell.grad.turb[mTurbQuants[var]][0].re = value;
                return fcell.grad.turb[mTurbQuants[var]][0].re;
            }
            if (var.canFind("visc_y")) {
                fcell.grad.turb[mTurbQuants[var]][1].re = value;
                return fcell.grad.turb[mTurbQuants[var]][1].re;
            }
            if (var.canFind("visc_z")) {
                fcell.grad.turb[mTurbQuants[var]][2].re = value;
                return fcell.grad.turb[mTurbQuants[var]][2].re;
            }
        }
        // For everything else, find appropriate case
        switch(var) {
        case "conv_x_rho": fcell.gradients.rho[0].re = value; return fcell.gradients.rho[0].re;
        case "conv_y_rho": fcell.gradients.rho[1].re = value; return fcell.gradients.rho[1].re;
        case "conv_z_rho": fcell.gradients.rho[2].re = value; return fcell.gradients.rho[2].re;
        case "conv_x_p": fcell.gradients.p[0].re = value; return fcell.gradients.p[0].re;
        case "conv_y_p": fcell.gradients.p[1].re = value; return fcell.gradients.p[1].re;
        case "conv_z_p": fcell.gradients.p[2].re = value; return fcell.gradients.p[2].re;
        case "conv_x_e": fcell.gradients.u[0].re = value; return fcell.gradients.u[0].re;
        case "conv_y_e": fcell.gradients.u[1].re = value; return fcell.gradients.u[1].re;
        case "conv_z_e": fcell.gradients.u[2].re = value; return fcell.gradients.u[2].re;
        case "conv_x_T": fcell.gradients.T[0].re = value; return fcell.gradients.T[0].re;
        case "conv_y_T": fcell.gradients.T[1].re = value; return fcell.gradients.T[1].re;
        case "conv_z_T": fcell.gradients.T[2].re = value; return fcell.gradients.T[2].re;
        case "visc_x_T": fcell.grad.T[0].re = value; return fcell.grad.T[0].re;
        case "visc_y_T": fcell.grad.T[1].re = value; return fcell.grad.T[1].re;
        case "visc_z_T": fcell.grad.T[2].re = value; return fcell.grad.T[2].re;
        case "conv_x_vel.x": fcell.gradients.velx[0].re = value; return fcell.gradients.velx[0].re;
        case "conv_y_vel.x": fcell.gradients.velx[1].re = value; return fcell.gradients.velx[1].re;
        case "conv_z_vel.x": fcell.gradients.velx[2].re = value; return fcell.gradients.velx[2].re;
        case "visc_x_vel.x": fcell.grad.vel[0][0].re = value; return fcell.grad.vel[0][0].re;
        case "visc_y_vel.x": fcell.grad.vel[0][1].re = value; return fcell.grad.vel[0][1].re;
        case "visc_z_vel.x": fcell.grad.vel[0][2].re = value; return fcell.grad.vel[0][2].re;
        case "conv_x_vel.y": fcell.gradients.vely[0].re = value; return fcell.gradients.vely[0].re;
        case "conv_y_vel.y": fcell.gradients.vely[1].re = value; return fcell.gradients.vely[1].re;
        case "conv_z_vel.y": fcell.gradients.vely[2].re = value; return fcell.gradients.vely[2].re;
        case "visc_x_vel.y": fcell.grad.vel[1][0].re = value; return fcell.grad.vel[1][0].re;
        case "visc_y_vel.y": fcell.grad.vel[1][1].re = value; return fcell.grad.vel[1][1].re;
        case "visc_z_vel.y": fcell.grad.vel[1][2].re = value; return fcell.grad.vel[1][2].re;
        case "conv_x_vel.z": fcell.gradients.velz[0].re = value; return fcell.gradients.velz[0].re;
        case "conv_y_vel.z": fcell.gradients.velz[1].re = value; return fcell.gradients.velz[1].re;
        case "conv_z_vel.z": fcell.gradients.velz[2].re = value; return fcell.gradients.velz[2].re;
        case "visc_x_vel.z": fcell.grad.vel[2][0].re = value; return fcell.grad.vel[2][0].re;
        case "visc_y_vel.z": fcell.grad.vel[2][1].re = value; return fcell.grad.vel[2][1].re;
        case "visc_z_vel.z": fcell.grad.vel[2][2].re = value; return fcell.grad.vel[2][2].re;
        default:
            throw new LmrException("Invalid selection for cell gradient variable: " ~ var);
        }
    } // end opIndexAssign()

private:
    int[string] mSpecies;
    int[string] mTModes;
    int[string] muModes;
    int[string] mTurbQuants;
} // end class FluidFVCellLimiterIO

FVCellIO createFVCellIO(FieldVarsType fvt, string[] variables)
{
    final switch (fvt) {
    case FieldVarsType.fluid: return new FluidFVCellIO(variables);
    case FieldVarsType.solid: return new SolidFVCellIO(variables);
    case FieldVarsType.limiter: return new FluidFVCellLimiterIO(variables);
    case FieldVarsType.residual: return new FluidFVCellResidualIO(variables);
    case FieldVarsType.gradient: return new FluidFVCellGradientIO(variables);
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
