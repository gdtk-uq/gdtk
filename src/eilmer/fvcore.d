/**
 * fvcore.d
 * Core definitions for finite-volume cells, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcore;

import std.conv;

class FlowSolverException : Exception {
    this(string message, string file=__FILE__, size_t line=__LINE__,
	 Throwable next=null)
    {
	super(message, file, line, next);
    }
}

// Symbolic names for the time-stepping schemes used to update the gasdynamic eqn.
enum GasdynamicUpdate {
    euler, 
    pc,
    midpoint, 
    classic_rk3,
    tvd_rk3,
    denman_rk3,
    moving_grid_1_stage,
    moving_grid_2_stage
}

string gasdynamic_update_scheme_name(GasdynamicUpdate gdut)
{
    final switch ( gdut ) {
    case GasdynamicUpdate.euler: return "euler";
    case GasdynamicUpdate.pc: return "predictor-corrector";
    case GasdynamicUpdate.midpoint: return "midpoint";
    case GasdynamicUpdate.classic_rk3: return "classic-rk3";
    case GasdynamicUpdate.tvd_rk3: return "tvd-rk3";
    case GasdynamicUpdate.denman_rk3: return "denman-rk3";
    case GasdynamicUpdate.moving_grid_1_stage: return "moving_grid_1_stage";
    case GasdynamicUpdate.moving_grid_2_stage: return "moving_grid_2_stage";
    }
} // end gasdynamic_update_scheme_name()

size_t number_of_stages_for_update_scheme(GasdynamicUpdate gdut)
{
    final switch (gdut) {
    case GasdynamicUpdate.euler: return 1;
    case GasdynamicUpdate.pc: return 2;
    case GasdynamicUpdate.midpoint: return 2;
    case GasdynamicUpdate.classic_rk3: return 3;
    case GasdynamicUpdate.tvd_rk3: return 3;
    case GasdynamicUpdate.denman_rk3: return 3;
    case GasdynamicUpdate.moving_grid_1_stage: return 1;
    case GasdynamicUpdate.moving_grid_2_stage: return 2;
    }
} // end number_of_stages_for_update_scheme()

size_t final_index_for_update_scheme(GasdynamicUpdate gdut)
{
    final switch (gdut) {
    case GasdynamicUpdate.euler: return 1;
    case GasdynamicUpdate.pc: return 2;
    case GasdynamicUpdate.midpoint: return 2;
    case GasdynamicUpdate.classic_rk3: return 3;
    case GasdynamicUpdate.tvd_rk3: return 3;
    case GasdynamicUpdate.denman_rk3: return 3;
    case GasdynamicUpdate.moving_grid_1_stage: return 1;
    case GasdynamicUpdate.moving_grid_2_stage: return 2;
    }
} // end final_index_for_update_scheme()

GasdynamicUpdate update_scheme_from_name(string name)
{
    switch (name) {
    case "euler": return GasdynamicUpdate.euler;
    case "pc": return GasdynamicUpdate.pc;
    case "predictor_corrector": return GasdynamicUpdate.pc;
    case "predictor-corrector": return GasdynamicUpdate.pc;
    case "midpoint": return GasdynamicUpdate.midpoint;
    case "classic_rk3": return GasdynamicUpdate.classic_rk3;
    case "classic-rk3": return GasdynamicUpdate.classic_rk3;
    case "tvd_rk3": return GasdynamicUpdate.tvd_rk3;
    case "tvd-rk3": return GasdynamicUpdate.tvd_rk3;
    case "denman_rk3": return GasdynamicUpdate.denman_rk3;
    case "denman-rk3": return GasdynamicUpdate.denman_rk3;
    case "moving_grid_1_stage": return GasdynamicUpdate.moving_grid_1_stage;
    case "moving-grid-1-stage": return GasdynamicUpdate.moving_grid_1_stage;
    case "moving_grid_2_stage": return GasdynamicUpdate.moving_grid_2_stage;
    case "moving-grid-2-stage": return GasdynamicUpdate.moving_grid_2_stage;
    default:
	string msg = text("Invalid gasdynamic update scheme name:", name);
	throw new FlowSolverException(msg);
    }
}  // end scheme_from_name()

// Symbolic names for grid motion types
enum GridMotion { none, user_defined, shock_fitting }
string grid_motion_name(GridMotion i)
{
    final switch (i) {
    case GridMotion.none: return "none";
    case GridMotion.user_defined: return "user_defined";
    case GridMotion.shock_fitting: return "shock_fitting";
    }
}

GridMotion grid_motion_from_name(string name)
{
    switch (name) {
    case "none": return GridMotion.none;
    case "user_defined": return GridMotion.user_defined;
    case "shock_fitting": return GridMotion.shock_fitting;
    default: return GridMotion.none;
    }
}


// [TODO] think about the following...
enum CopyDataOption { all, minimal_flow, all_flow, grid, cell_lengths_only }

// Minimum values for turbulent kinetic energy (m^2/s^2) and frequency (1/s)
// for applying limiters in the k-omega model.
enum
    small_tke = 0.1,
    small_omega = 1.0;

// Symbolic names for the types of flow-data reconstruction.
enum InterpolateOption { pt, rhoe, rhop, rhot }

string thermo_interpolator_name(InterpolateOption i)
{
    final switch ( i ) {
    case InterpolateOption.pt: return "pT";
    case InterpolateOption.rhoe: return "rhoe";
    case InterpolateOption.rhop: return "rhop";
    case InterpolateOption.rhot: return "rhoT";
    }
} // end thermo_interpolator_name()

InterpolateOption thermo_interpolator_from_name(string name)
{
    switch ( name ) {
    case "pT": return InterpolateOption.pt;
    case "pt": return InterpolateOption.pt;
    case "rhoe": return InterpolateOption.rhoe;
    case "rhop": return InterpolateOption.rhop;
    case "rhoT": return InterpolateOption.rhot;
    case "rhot": return InterpolateOption.rhot;
    default: return InterpolateOption.rhoe;
    }
} // end thermo_interpolator_from_name()

// Symbolic names for the flavours of our flux_calculators.
enum FluxCalculator {
    ausmdv, // Wada and Liou's flux calculator AIAA Paper 94-0083
    efm, // Mike Macrossan's EFM flux calculation
    ausm_plus_up, // Liou's 2006 all-speed flux calculator
    adaptive, // EFM near shocks, AUSMDV otherwise
    hlle // MHD HLLE approximate Riemann solver
}

string flux_calculator_name(FluxCalculator fcalc)
{
    final switch ( fcalc ) {
    case FluxCalculator.ausmdv: return "ausmdv";
    case FluxCalculator.efm: return "efm";
    case FluxCalculator.ausm_plus_up: return "ausm_plus_up";
    case FluxCalculator.adaptive: return "adaptive";
    case FluxCalculator.hlle: return "hlle";
    }
}

FluxCalculator flux_calculator_from_name(string name)
{
    switch ( name ) {
    case "ausmdv": return FluxCalculator.ausmdv;
    case "efm": return FluxCalculator.efm;
    case "ausm_plus_up": return FluxCalculator.ausm_plus_up;
    case "adaptive": return FluxCalculator.adaptive;
    case "hlle": return FluxCalculator.hlle;
    default:
	string msg = text("Invalid flux calculator name:", name);
	throw new FlowSolverException(msg);
    }
}

// Symbolic names for the flavours of spatial-derivative calculators.
enum SpatialDerivCalc {
    least_squares,
    divergence
}

string spatial_deriv_calc_name(SpatialDerivCalc sdc)
{
    final switch ( sdc ) {
    case SpatialDerivCalc.least_squares: return "least_squares";
    case SpatialDerivCalc.divergence: return "divergence";
    }
}

SpatialDerivCalc spatial_deriv_calc_from_name(string name)
{
    switch ( name ) {
    case "least_squares": return SpatialDerivCalc.least_squares;
    case "divergence": return SpatialDerivCalc.divergence;
    default:
	string msg = text("Invalid spatial-derivative calculator name:", name);
	throw new FlowSolverException(msg);
    }
}

enum SpatialDerivLocn {
    faces,
    vertices
}

string spatial_deriv_locn_name(SpatialDerivLocn sdl)
{
    final switch ( sdl ) {
    case SpatialDerivLocn.faces: return "faces";
    case SpatialDerivLocn.vertices: return "vertices";
    }
}

SpatialDerivLocn spatial_deriv_locn_from_name(string name)
{
    switch ( name ) {
    case "faces": return SpatialDerivLocn.faces;
    case "vertices": return SpatialDerivLocn.vertices;
    default:
	string msg = text("Invalid spatial-derivative location name:", name);
	throw new FlowSolverException(msg);
    }
}
