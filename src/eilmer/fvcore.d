/**
 * fvcore.d
 * Core definitions for finite-volume cells, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcore;

class FlowSolverException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

struct FlowStateLimits {
    double max_velocity = 30000.0; // m/s
    double max_tke = 0.01 * double.max;
    double min_tke = 0.0;
    double max_temp = 50000.0; // Kelvin
    double min_temp = 0.0;
    double min_pressure = 0.1; // Pascals
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

@nogc
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

@nogc
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

@nogc
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

@nogc
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
        throw new FlowSolverException("Invalid gasdynamic update scheme name");
    }
}  // end scheme_from_name()

// Symbolic names for grid motion types
enum GridMotion { none, user_defined, shock_fitting }

@nogc
string grid_motion_name(GridMotion i)
{
    final switch (i) {
    case GridMotion.none: return "none";
    case GridMotion.user_defined: return "user_defined";
    case GridMotion.shock_fitting: return "shock_fitting";
    }
}

@nogc
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
enum InterpolateOption { pt, rhou, rhop, rhot }

@nogc
string thermo_interpolator_name(InterpolateOption i)
{
    final switch ( i ) {
    case InterpolateOption.pt: return "pT";
    case InterpolateOption.rhou: return "rhou";
    case InterpolateOption.rhop: return "rhop";
    case InterpolateOption.rhot: return "rhoT";
    }
} // end thermo_interpolator_name()

@nogc
InterpolateOption thermo_interpolator_from_name(string name)
{
    switch ( name ) {
    case "pT": return InterpolateOption.pt;
    case "pt": return InterpolateOption.pt;
    case "rhou": return InterpolateOption.rhou;
    case "rhoe": return InterpolateOption.rhou; // allow the old name
    case "rhop": return InterpolateOption.rhop;
    case "rhoT": return InterpolateOption.rhot;
    case "rhot": return InterpolateOption.rhot;
    default: return InterpolateOption.rhou;
    }
} // end thermo_interpolator_from_name()

// Symbolic names for the flavours of our flux_calculators.
enum FluxCalculator {
    ausmdv, // Wada and Liou's flux calculator AIAA Paper 94-0083
    hanel, // Hanel's flux calculator (details in Wada & Lious's 1997 SIAM paper)
    efm, // Mike Macrossan's EFM flux calculation
    ausm_plus_up, // Liou's 2006 all-speed flux calculator
    adaptive_efm_ausmdv, // EFM near shocks, AUSMDV otherwise
    adaptive_hanel_ausmdv, // Hanel near shocks, AUSMDV otherwise
    adaptive_hlle_roe, // HLLE near shocks, Roe otherwise
    hlle, // MHD HLLE approximate Riemann solver
    roe // Roe approximate Riemann solver
}

@nogc
string flux_calculator_name(FluxCalculator fcalc)
{
    final switch ( fcalc ) {
    case FluxCalculator.ausmdv: return "ausmdv";
    case FluxCalculator.hanel: return "hanel";
    case FluxCalculator.efm: return "efm";
    case FluxCalculator.ausm_plus_up: return "ausm_plus_up";
    case FluxCalculator.adaptive_efm_ausmdv: return "adaptive_efm_ausmdv";
    case FluxCalculator.adaptive_hanel_ausmdv: return "adaptive_hanel_ausmdv";
    case FluxCalculator.adaptive_hlle_roe: return "adaptive_hlle_roe";
    case FluxCalculator.hlle: return "hlle";
    case FluxCalculator.roe: return "roe";
    }
}

@nogc
FluxCalculator flux_calculator_from_name(string name)
{
    switch ( name ) {
    case "ausmdv": return FluxCalculator.ausmdv;
    case "hanel": return FluxCalculator.hanel;
    case "efm": return FluxCalculator.efm;
    case "ausm_plus_up": return FluxCalculator.ausm_plus_up;
    case "adaptive_efm_ausmdv": return FluxCalculator.adaptive_efm_ausmdv;
    case "adaptive_hanel_ausmdv": return FluxCalculator.adaptive_hanel_ausmdv;
    case "adaptive": return FluxCalculator.adaptive_efm_ausmdv;
    case "adaptive_hlle_roe": return FluxCalculator.adaptive_hlle_roe;
    case "hlle": return FluxCalculator.hlle;
    case "roe": return FluxCalculator.roe;
    default:
        throw new FlowSolverException("Invalid flux calculator name");
    }
}

// Symbolic names for the flavours of spatial-derivative calculators.
enum SpatialDerivCalc {
    least_squares,
    divergence,
}

@nogc
string spatial_deriv_calc_name(SpatialDerivCalc sdc)
{
    final switch ( sdc ) {
    case SpatialDerivCalc.least_squares: return "least_squares";
    case SpatialDerivCalc.divergence: return "divergence";
    }
}

@nogc
SpatialDerivCalc spatial_deriv_calc_from_name(string name)
{
    switch ( name ) {
    case "least_squares": return SpatialDerivCalc.least_squares;
    case "divergence": return SpatialDerivCalc.divergence;
    default:
        throw new FlowSolverException("Invalid spatial-derivative calculator name");
    }
}

enum SpatialDerivLocn {
    faces,
    vertices,
    cells
}

@nogc
string spatial_deriv_locn_name(SpatialDerivLocn sdl)
{
    final switch ( sdl ) {
    case SpatialDerivLocn.faces: return "faces";
    case SpatialDerivLocn.vertices: return "vertices";
    case SpatialDerivLocn.cells: return "cells";
    }
}

@nogc
SpatialDerivLocn spatial_deriv_locn_from_name(string name)
{
    switch ( name ) {
    case "faces": return SpatialDerivLocn.faces;
    case "vertices": return SpatialDerivLocn.vertices;
    case "cells": return SpatialDerivLocn.cells;
    default:
        throw new FlowSolverException("Invalid spatial-derivative location name");
    }
}

// Symbolic names for the flavours of unstructured limiters.
enum UnstructuredLimiter {
    van_albada,
    min_mod,
    mlp,
    barth,
    heuristic_van_albada,
    venkat
}

@nogc
string unstructured_limiter_name(UnstructuredLimiter ul)
{
    final switch ( ul ) {
    case UnstructuredLimiter.van_albada: return "van_albada";
    case UnstructuredLimiter.min_mod: return "min_mod";
    case UnstructuredLimiter.mlp: return "mlp";
    case UnstructuredLimiter.barth: return "barth";
    case UnstructuredLimiter.heuristic_van_albada: return "heuristic_van_albada";
    case UnstructuredLimiter.venkat: return "venkat";
    }
}

@nogc
UnstructuredLimiter unstructured_limiter_from_name(string name)
{
    switch ( name ) {
    case "van_albada": return UnstructuredLimiter.van_albada;
    case "min_mod": return UnstructuredLimiter.min_mod;
    case "mlp": return UnstructuredLimiter.mlp;
    case "barth": return UnstructuredLimiter.barth;
    case "heuristic_van_albada": return UnstructuredLimiter.heuristic_van_albada;
    case "venkat": return UnstructuredLimiter.venkat;
    default:
        throw new FlowSolverException("Invalid unstructured limiter name");
    }
}

// Symbolic names for the residual smoothing
enum ResidualSmoothingType {
    explicit,
    implicit
}

@nogc
string residual_smoothing_type_name(ResidualSmoothingType rs)
{
    final switch ( rs ) {
    case ResidualSmoothingType.explicit: return "explicit";
    case ResidualSmoothingType.implicit: return "implicit";
    }
}

@nogc
ResidualSmoothingType residual_smoothing_type_from_name(string name)
{
    switch ( name ) {
    case "explicit": return ResidualSmoothingType.explicit;
    case "implicit": return ResidualSmoothingType.implicit;
    default:
        throw new FlowSolverException("Invalid residual smoothing type name");
    }
}
