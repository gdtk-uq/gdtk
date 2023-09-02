/** globalconfig.d
 * A place to keep the configuration details for the simulation.
 * This is a long file, in several parts:
 *   1. Symbolic names and enums
 *   2. Special zones
 *   3. The structs and classes defining configuration data.
 *   4. Reading from JSON file
 *   5. Lua interaction
 *
 * Authors: Peter J, Rowan G, Kyle D, Nick G and Daryl B.
 * Versions:
 *   2014-07-18: First code
 *   2015-02-05: Resumed work
 *   2021-05-24: Absorbed: fvcore.d and remnants of luaglobalconfig.d
 */

module globalconfig;

import std.conv;
import std.stdio;
import std.string;
import std.typecons;
import core.stdc.stdlib : exit;
import std.json;
import std.file;
import std.array;
import std.format;

import util.lua;
import util.lua_service;
import nm.luabbla;
import nm.schedule;
import lua_helper;
import gas;
import gas.luagas_model;
import kinetics;
import geom;
import geom.luawrap;
version (opencl_gpu_chem) {
    import opencl_gpu_chem;
}
version (cuda_gpu_chem) {
     import cuda_gpu_chem;
}
import json_helper;
import globaldata;
import flowstate;
import conservedquantities;
import fluidblock;
import fluidblockarray;
import fluidblockio_old;
import sfluidblock: SFluidBlock;
import ufluidblock: UFluidBlock;
import ssolidblock;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import grid_motion;
import grid_motion_udf;
import mass_diffusion;
import loads;
import turbulence;

// --------------------------------
// PART 1. Symbolic names and enums
// --------------------------------

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
    rkl1,
    rkl2,
    moving_grid_1_stage,
    moving_grid_2_stage,
    backward_euler,
    implicit_rk1,
    classic_rk4
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
    case GasdynamicUpdate.rkl1: return "rkl1";
    case GasdynamicUpdate.rkl2: return "rkl2";
    case GasdynamicUpdate.moving_grid_1_stage: return "moving_grid_1_stage";
    case GasdynamicUpdate.moving_grid_2_stage: return "moving_grid_2_stage";
    case GasdynamicUpdate.backward_euler: return "backward_euler";
    case GasdynamicUpdate.implicit_rk1: return "implicit_rk1";
    case GasdynamicUpdate.classic_rk4: return "classic_rk4";
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
    case GasdynamicUpdate.rkl1: return 3;
    case GasdynamicUpdate.rkl2: return 3;
    case GasdynamicUpdate.moving_grid_1_stage: return 1;
    case GasdynamicUpdate.moving_grid_2_stage: return 2;
    case GasdynamicUpdate.backward_euler: return 1;
    case GasdynamicUpdate.implicit_rk1: return 1;
    case GasdynamicUpdate.classic_rk4: return 4;
    }
} // end number_of_stages_for_update_scheme()

@nogc
bool is_explicit_update_scheme(GasdynamicUpdate gdut)
{
    final switch (gdut) {
    case GasdynamicUpdate.euler: return true;
    case GasdynamicUpdate.pc: return true;
    case GasdynamicUpdate.midpoint: return true;
    case GasdynamicUpdate.classic_rk3: return true;
    case GasdynamicUpdate.tvd_rk3: return true;
    case GasdynamicUpdate.denman_rk3: return true;
    case GasdynamicUpdate.rkl1: return true;
    case GasdynamicUpdate.rkl2: return true;
    case GasdynamicUpdate.moving_grid_1_stage: return true;
    case GasdynamicUpdate.moving_grid_2_stage: return true;
    case GasdynamicUpdate.backward_euler: return false;
    case GasdynamicUpdate.implicit_rk1: return false;
    case GasdynamicUpdate.classic_rk4: return true;
    }
} // is_explicit_update_scheme()

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
    case GasdynamicUpdate.rkl1: return 3;
    case GasdynamicUpdate.rkl2: return 3;
    case GasdynamicUpdate.moving_grid_1_stage: return 1;
    case GasdynamicUpdate.moving_grid_2_stage: return 2;
    case GasdynamicUpdate.backward_euler: return 1;
    case GasdynamicUpdate.implicit_rk1: return 1;
    case GasdynamicUpdate.classic_rk4: return 4;
    }
} // end final_index_for_update_scheme()

GasdynamicUpdate update_scheme_from_name(string name)
{
    string name_ = name.replace("-", "_");
    switch (name_) {
    case "euler": return GasdynamicUpdate.euler;
    case "pc": return GasdynamicUpdate.pc;
    case "predictor_corrector": return GasdynamicUpdate.pc;
    case "midpoint": return GasdynamicUpdate.midpoint;
    case "classic_rk3": return GasdynamicUpdate.classic_rk3;
    case "tvd_rk3": return GasdynamicUpdate.tvd_rk3;
    case "denman_rk3": return GasdynamicUpdate.denman_rk3;
    case "rkl1": return GasdynamicUpdate.rkl1;
    case "rkl2": return GasdynamicUpdate.rkl2;
    case "moving_grid_1_stage": return GasdynamicUpdate.moving_grid_1_stage;
    case "moving_grid_2_stage": return GasdynamicUpdate.moving_grid_2_stage;
    case "backward_euler": return GasdynamicUpdate.backward_euler;
    case "implicit_rk1": return GasdynamicUpdate.implicit_rk1;
    case "classic_rk4": return GasdynamicUpdate.classic_rk4;
    default:
        string msg = format("Invalid gasdynamic update scheme '%s' selected.", name);
        throw new Error(msg);
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

GridMotion grid_motion_from_name(string name)
{
    string name_ = name.replace("-", "_");
    switch (name_) {
    case "none": return GridMotion.none;
    case "user_defined": return GridMotion.user_defined;
    case "shock_fitting": return GridMotion.shock_fitting;
    default: return GridMotion.none;
    }
}


// [TODO] think about the following...
enum CopyDataOption { all, minimal_flow, all_flow, grid, cell_lengths_only }

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
    hllc,  // HLLC approximate Riemann solver (Details from Toro's textbook on Riemann solvers)
    ldfss0,  // A Low-diffusion flux-splitting scheme (Details in Edward's Computers & Fluids paper 1997)
    ldfss2,  // A Low-diffusion flux-splitting scheme (Details in Edward's Computers & Fluids paper 1997)
    hanel, // Hanel's flux calculator (details in Wada & Lious's 1997 SIAM paper)
    efm, // Mike Macrossan's EFM flux calculation
    ausm_plus_up, // Liou's 2006 all-speed flux calculator
    adaptive_efm_ausmdv, // EFM near shocks, AUSMDV otherwise
    adaptive_hanel_ausmdv, // Hanel near shocks, AUSMDV otherwise
    adaptive_hanel_ausm_plus_up, // Hanel near shocks, AUSM+up otherwise
    adaptive_ldfss0_ldfss2, // LDFSS0 near shocks, LDFSS2 otherwise
    adaptive_hlle_hllc, // Standard HLLE near shocks, HLLC otherwise
    adaptive_hlle_roe, // HLLE near shocks, Roe otherwise
    hlle, // MHD HLLE approximate Riemann solver
    hlle2, // Standard HLLE approximate Riemann solver
    roe, // Roe approximate Riemann solver
    osher, // Osher approximate Riemann solver
    asf, // Alpha-split Flux, high-order central difference. See Fisher et al. "Discete Conservative Finite-Difference Formulations..." 2012.
    adaptive_ausmdv_asf // AUSMDV near shocks, asf otherwise
}

@nogc
string flux_calculator_name(FluxCalculator fcalc)
{
    final switch ( fcalc ) {
    case FluxCalculator.ausmdv: return "ausmdv";
    case FluxCalculator.hllc: return "hllc";
    case FluxCalculator.ldfss0: return "ldfss0";
    case FluxCalculator.ldfss2: return "ldfss2";
    case FluxCalculator.hanel: return "hanel";
    case FluxCalculator.efm: return "efm";
    case FluxCalculator.ausm_plus_up: return "ausm_plus_up";
    case FluxCalculator.adaptive_efm_ausmdv: return "adaptive_efm_ausmdv";
    case FluxCalculator.adaptive_hanel_ausmdv: return "adaptive_hanel_ausmdv";
    case FluxCalculator.adaptive_hanel_ausm_plus_up: return "adaptive_hanel_ausm_plus_up";
    case FluxCalculator.adaptive_ldfss0_ldfss2: return "adaptive_ldfss0_ldfss2";
    case FluxCalculator.adaptive_hlle_hllc: return "adaptive_hlle_hllc";
    case FluxCalculator.adaptive_hlle_roe: return "adaptive_hlle_roe";
    case FluxCalculator.hlle: return "hlle";
    case FluxCalculator.hlle2: return "hlle2";
    case FluxCalculator.roe: return "roe";
    case FluxCalculator.osher: return "osher";
    case FluxCalculator.asf: return "asf";
    case FluxCalculator.adaptive_ausmdv_asf: return "adaptive_ausmdv_asf";
    }
}

FluxCalculator flux_calculator_from_name(string name)
{
    switch ( name ) {
    case "ausmdv": return FluxCalculator.ausmdv;
    case "hllc": return FluxCalculator.hllc;
    case "ldfss0": return FluxCalculator.ldfss0;
    case "ldfss2": return FluxCalculator.ldfss2;
    case "hanel": return FluxCalculator.hanel;
    case "efm": return FluxCalculator.efm;
    case "ausm_plus_up": return FluxCalculator.ausm_plus_up;
    case "adaptive_efm_ausmdv": return FluxCalculator.adaptive_efm_ausmdv;
    case "adaptive_hanel_ausmdv": return FluxCalculator.adaptive_hanel_ausmdv;
    case "adaptive_hanel_ausm_plus_up": return FluxCalculator.adaptive_hanel_ausm_plus_up;
    case "adaptive_ldfss0_ldfss2": return FluxCalculator.adaptive_ldfss0_ldfss2;
    case "adaptive_hlle_hllc": return FluxCalculator.adaptive_hlle_hllc;
    case "adaptive": return FluxCalculator.adaptive_efm_ausmdv;
    case "adaptive_hlle_roe": return FluxCalculator.adaptive_hlle_roe;
    case "hlle": return FluxCalculator.hlle;
    case "hlle2": return FluxCalculator.hlle2;
    case "roe": return FluxCalculator.roe;
    case "osher": return FluxCalculator.osher;
    case "asf": return FluxCalculator.asf;
    case "adaptive_ausmdv_asf" : return FluxCalculator.adaptive_ausmdv_asf;
    default:
        string msg = format("Invalid flux calculator '%s' selected.", name);
        throw new Error(msg);
    }
}

// Symbolic names for mode of chemistry update in transient solver.
enum ChemistryUpdateMode { split, integral }

@nogc
string chemistry_update_mode_name(ChemistryUpdateMode mode)
{
    final switch (mode) {
    case ChemistryUpdateMode.split: return "split";
    case ChemistryUpdateMode.integral: return "integral";
    }
}

@nogc
ChemistryUpdateMode chemistry_update_mode_from_name(string name)
{
    switch (name) {
    case "split": return ChemistryUpdateMode.split;
    case "integral": return ChemistryUpdateMode.integral;
    default:
        throw new FlowSolverException("Invalid chemistry-update mode.");
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
    case "least-squares": return SpatialDerivCalc.least_squares;
    case "leastsquares": return SpatialDerivCalc.least_squares;
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
    svan_albada,
    min_mod,
    barth,
    park,
    hvan_albada,
    van_albada,
    hnishikawa,
    nishikawa,
    hvenkat_mlp,
    venkat_mlp,
    hvenkat,
    venkat
}

@nogc
string unstructured_limiter_name(UnstructuredLimiter ul)
{
    final switch ( ul ) {
    case UnstructuredLimiter.svan_albada: return "svan_albada";
    case UnstructuredLimiter.min_mod: return "min_mod";
    case UnstructuredLimiter.barth: return "barth";
    case UnstructuredLimiter.park: return "park";
    case UnstructuredLimiter.hvan_albada: return "hvan_albada";
    case UnstructuredLimiter.van_albada: return "van_albada";
    case UnstructuredLimiter.hnishikawa: return "hnishikawa";
    case UnstructuredLimiter.nishikawa: return "nishikawa";
    case UnstructuredLimiter.hvenkat_mlp: return "hvenkat_mlp";
    case UnstructuredLimiter.venkat_mlp: return "venkat_mlp";
    case UnstructuredLimiter.hvenkat: return "hvenkat";
    case UnstructuredLimiter.venkat: return "venkat";
    }
}

@nogc
UnstructuredLimiter unstructured_limiter_from_name(string name)
{
    switch ( name ) {
    case "svan_albada": return UnstructuredLimiter.svan_albada;
    case "min_mod": return UnstructuredLimiter.min_mod;
    case "barth": return UnstructuredLimiter.barth;
    case "park": return UnstructuredLimiter.park;
    case "hvan_albada": return UnstructuredLimiter.hvan_albada;
    case "van_albada": return UnstructuredLimiter.van_albada;
    case "hnishikawa": return UnstructuredLimiter.hnishikawa;
    case "nishikawa": return UnstructuredLimiter.nishikawa;
    case "hvenkat_mlp": return UnstructuredLimiter.hvenkat_mlp;
    case "venkat_mlp": return UnstructuredLimiter.venkat_mlp;
    case "hvenkat": return UnstructuredLimiter.hvenkat;
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

// Symbolic names for the shock detectors
enum ShockDetector {
    PJ,
}

@nogc
string shock_detector_name(ShockDetector sd)
{
    final switch ( sd ) {
    case ShockDetector.PJ: return "PJ";
    }
}

@nogc
ShockDetector shock_detector_from_name(string name)
{
    switch ( name ) {
    case "PJ": return ShockDetector.PJ;
    default:
        throw new FlowSolverException("Invalid shock detector name: ");
    }
}

enum StrangSplittingMode { full_T_full_R, half_R_full_T_half_R }
@nogc
string strangSplittingModeName(StrangSplittingMode i)
{
    final switch(i) {
    case StrangSplittingMode.full_T_full_R: return "full_T_full_R";
    case StrangSplittingMode.half_R_full_T_half_R: return "half_R_full_T_half_R";
    }
}

@nogc
StrangSplittingMode strangSplittingModeFromName(string name)
{
    switch(name) {
    case "full_T_full_R": return StrangSplittingMode.full_T_full_R;
    case "half_R_full_T_half_R": return StrangSplittingMode.half_R_full_T_half_R;
    default: return StrangSplittingMode.full_T_full_R;
    }
}

// Symbolic names for JJ Hoste's Turbulence-Chemistry-Interaction Model
enum TCIModel { none, edm, edc }

@nogc
string tci_model_name(TCIModel tcim)
{
    final switch (tcim) {
    case TCIModel.none: return "none";
    case TCIModel.edm: return "edm";
    case TCIModel.edc: return "edc";
    }
} // end tci_model_name()

@nogc
TCIModel tci_model_from_name(string name)
{
    switch (name) {
    case "none": return TCIModel.none;
    case "edm": return TCIModel.edm;
    case "edc": return TCIModel.edc;
    default: return TCIModel.none;
    }
} // end tci_model_from_name()


enum SolidDomainCoupling { tight, loose, lagged, steady_fluid_transient_solid }

@nogc
string solidDomainCouplingName(SolidDomainCoupling i)
{
    final switch (i) {
    case SolidDomainCoupling.tight: return "tight";
    case SolidDomainCoupling.loose: return "loose";
    case SolidDomainCoupling.lagged: return "lagged";
    case SolidDomainCoupling.steady_fluid_transient_solid: return "steady_fluid_transient_solid";
    }
}

@nogc
SolidDomainCoupling solidDomainCouplingFromName(string name)
{
    switch (name) {
    case "tight": return SolidDomainCoupling.tight;
    case "loose": return SolidDomainCoupling.loose;
    case "lagged": return SolidDomainCoupling.lagged;
    case "steady_fluid_transient_solid": return SolidDomainCoupling.steady_fluid_transient_solid;
    default: return SolidDomainCoupling.tight;
    }
}

enum PreconditionMatrixType { diagonal, jacobi, sgs, ilu }

string preconditionMatrixTypeName(PreconditionMatrixType i)
{
    final switch (i) {
    case PreconditionMatrixType.diagonal: return "diagonal";
    case PreconditionMatrixType.jacobi: return "jacobi";
    case PreconditionMatrixType.sgs: return "sgs";
    case PreconditionMatrixType.ilu: return "ilu";
    }
} // end preconditionMatrixTypeName()

PreconditionMatrixType preconditionMatrixTypeFromName(string name)
{
    switch (name) {
    case "diagonal": return PreconditionMatrixType.diagonal;
    case "jacobi": return PreconditionMatrixType.jacobi;
    case "sgs": return PreconditionMatrixType.sgs;
    case "ilu": return PreconditionMatrixType.ilu;
    default:
        string errMsg = "The selected 'preconditioner' is unavailable.\n";
        errMsg ~= format("You selected: '%s'\n", name);
        errMsg ~= "The available strategies are: \n";
        errMsg ~= "   'lu_sgs'\n";
        errMsg ~= "   'diagonal'\n";
        errMsg ~= "   'jacobi'\n";
        errMsg ~= "   'sgs'\n";
        errMsg ~= "   'ilu'\n";
        errMsg ~= "Check your selection or its spelling in the input file.\n";
        throw new Error(errMsg);
    }
} // end preconditionMatrixTypeFromName()

enum EtaStrategy { constant, geometric, adaptive, adaptive_capped }
string etaStrategyName(EtaStrategy i)
{
    final switch (i) {
    case EtaStrategy.constant: return "constant";
    case EtaStrategy.geometric: return "geometric";
    case EtaStrategy.adaptive: return "adaptive";
    case EtaStrategy.adaptive_capped: return "adaptive_capped";
    }
} // end etaStrategyName()

EtaStrategy etaStrategyFromName(string name)
{
    switch (name) {
    case "constant": return EtaStrategy.constant;
    case "geometric": return EtaStrategy.geometric;
    case "adaptive": return EtaStrategy.adaptive;
    case "adaptive_capped": return EtaStrategy.adaptive_capped;
    default:
        string errMsg = "The selected 'eta_strategy' is unavailable.\n";
        errMsg ~= format("You selected: '%s'\n", name);
        errMsg ~= "The available strategies are: \n";
        errMsg ~= "   'constant'\n";
        errMsg ~= "   'geometric'\n";
        errMsg ~= "   'adaptive'\n";
        errMsg ~= "   'adaptive_capped'\n";
        errMsg ~= "Check your selection or its spelling in the input file.\n";
        throw new Error(errMsg);
    }
} // end etaStrategyFromName()


// ---------------------
// PART 2. Special zones
// ---------------------

class BlockZone {
    // Note that these data are not shared across threads
    // because we will want to access them frequently at the lower levels of the code.
    Vector3 p0, p1;
    this(in Vector3 p0, in Vector3 p1) {
        this.p0 = p0;
        this.p1 = p1;
    }
    this(const BlockZone other)
    {
        p0 = other.p0;
        p1 = other.p1;
    }
    override string toString() {
        return text("BlockZone(p0=", to!string(p0), ", p1=", to!string(p1), ")");
    }
    @nogc
    bool is_inside(ref const(Vector3) p, int dimensions) {
        if ( p.x >= p0.x && p.x <= p1.x &&
             p.y >= p0.y && p.y <= p1.y ) {
            if ( dimensions == 2 ) {
                return true;
            } else if ( p.z >= p0.z && p.z <= p1.z ) {
                return true;
            }
        }
        return false;
    } // end is_inside()
}

class IgnitionZone : BlockZone {
    double Tig; // temperature to apply within reaction_update ensure ignition
    this(in Vector3 p0, in Vector3 p1, double Tig) {
        super(p0, p1);
        this.Tig = Tig;
    }
    this(const IgnitionZone other)
    {
        super(other.p0, other.p1);
        Tig = other.Tig;
    }
    override string toString() {
        return text("IgnitionZone(p0=", to!string(p0), ", p1=", to!string(p1),
                    ", Tig=", Tig, ")");
    }
} // end class IgnitionZone



// -----------------------------------------------------------
// PART 3. The structs and classes defining configuration data
// -----------------------------------------------------------

struct SolidDomainLooseUpdateOptions {
    int maxNewtonIterations = 10;
    double toleranceNewtonUpdate = 1.0e-2;
    int maxGMRESIterations = 10;
    double toleranceGMRESSolve = 1.0e-3;
    double perturbationSize = 1.0e-2;
}

struct ShapeSensitivityCalculatorOptions {
    bool pseudotime = false;
    int pseudotime_lhs_jacobian_order = 1;
    int adjoint_precondition_matrix_order = 0;
    bool read_frozen_limiter_values_from_file = false;
    // sensitivity parameters
    double epsilon = 1.0e-30;
    // GMRES parameters
    int maxOuterIterations = 10;
    int maxRestarts = 10;
    double cfl0 = 1.0; // initial CFL
    double eta = 0.1; // inner loop tolerance
    double stopOnRelativeGlobalResidual = 1.0e-16; // outer loop tolerance
    // bezier curve fit parameters
    double tolBezierCurveFit = 1.0e-06;
    int maxStepsBezierCurveFit = 10000;
    string userDefinedObjectiveFile = "";
}

struct SteadyStateSolverOptions {
    int nConserved = 4;
    bool usePreconditioner = true;
    int frozenPreconditionerCount = 1;
    int startPreconditioning = 1;
    int iluFill = 0;
    PreconditionMatrixType preconditionMatrixType = PreconditionMatrixType.jacobi;
    FluxCalculator preconditionMatrixFluxCalculator = FluxCalculator.adaptive_hanel_ausmdv;
    double preconditionerSigma = 1.0e-30;
    bool frozenLimiterOnLHS = false;
    bool useAdaptivePreconditioner = false;
    bool usePhysicalityCheck = false;
    double physicalityCheckTheta = 0.2;
    bool useLineSearch = false;
    bool inviscidCFL = false;
    bool useScaling = true;
    bool useComplexMatVecEval = false;
    int temporalIntegrationMode = 0;
    int nPreSteps = 10;
    int nTotalSteps = 100;
    int maxNumberAttempts = 3; // at taking a Newton step.
    double stopOnRelGlobalResid = 1.0e-99;
    double stopOnAbsGlobalResid = 1.0e-99;
    double stopOnMassBalance = -1.0;
    // Restarted preconditioned FGMRES settings
    int maxSubIterations = 1;
    int maxOuterIterations = 10;
    int maxRestarts = 10;
    int nInnerIterations = 5;
    double cfl_max = 1e8;
    double cfl_min = 1e-02;
    bool include_turb_quantities_in_residual = true;
    bool residual_based_cfl_scheduling = true;
    int cfl_schedule_length = 0;
    double[] cfl_schedule_value_list;
    int[] cfl_schedule_iter_list;
    // Options for start-up phase
    int nStartUpSteps = 5;
    int LHSeval0 = 1;
    int RHSeval0 = 1;
    double p0 = 0.75;
    double cfl0 = 1.0;
    double eta0 = 0.5;
    double tau0 = 0.1;
    double sigma0 = 1.0e-8;
    // Options for inexact Newton phase
    int LHSeval1 = 2;
    int RHSeval1 = 2;
    double p1 = 1.0;
    double cfl1 = 10.0;
    double tau1 = 0.1;
    double sigma1 = 1.0e-8;
    EtaStrategy etaStrategy = EtaStrategy.constant;
    double eta1 = 0.5;
    double eta1_max = 0.9;
    double eta1_min = 0.01;
    double etaRatioPerStep = 0.9;
    double gamma = 0.9;
    double alpha = 2.0;
    double limiterFreezingResidReduction = 1e-99;
    int limiterFreezingCount = 50;
    // Options related to writing out snapshots and diagnostics
    int snapshotsCount = 10;
    int nTotalSnapshots = 5;
    int writeDiagnosticsCount = 20;
    int writeLoadsCount = 20;
} // end struct SteadyStateSolverOptions


final class GlobalConfig {
    shared static bool in_mpi_context = false; // Usual context is thread-parallel only.
    shared static int mpi_size = 0; // Number of MPI tasks, overall.
    shared static int mpi_rank_for_local_task = 0;
    shared static bool is_master_task = true; // In an MPI run, only one task will be master.
    shared static int[] mpi_rank_for_block; // To know where each block has been assigned.
    shared static int[] localFluidBlockIds; // We will search this array to see if the fluid block is local.
    shared static int[] localSolidBlockIds; // We will search this array to see if the solid block is local.
    //
    shared static string base_file_name = "job"; // Change this to suit at run time.
    shared static string grid_format = "gziptext"; // alternative is "rawbinary"
    shared static string flow_format = "gziptext";
    shared static bool new_flow_format = false; // enable/disable new flow format (to be deprecated)
    // Depending on the format of the contained data, grid and solution files will have
    // a particular file extension.
    shared static string gridFileExt = "gz";
    shared static string flowFileExt = "gz";
    shared static string title = "Eilmer4 simulation"; // Change this to suit at run time.
    shared static string gas_model_file = "gas-model.lua";
    // Note that the following reference to the GasModel is NOT shared.
    // Do not use this one from within a thread.
    static GasModel gmodel_master; // to provide GasModel access to the Lua domain
    // Presumably, we won't be accessing this particular gas model from the
    // individual block computations, so that parallel computations for the blocks
    // don't trip over each other.
    static uint n_species;
    static uint n_heavy;
    static uint n_modes;
    static bool sticky_electrons = false; // Default to electrons being a separate species.
    // Setting sticky_electrons=true will cause the electron mass fraction to be reset after every
    // timestep to preserve charge neutrality in each cell. Note that a transport equation for the
    // electrons is still solved, even though its output is effectively discarded.
    //
    // Customization of the simulation is via user-defined actions.
    shared static string udf_supervisor_file; // empty to start
    // A scratch-pad area for the user-defined functions.
    // This will allow use to persist some arbitrary user data
    // between calls to the master Lua interpreter.
    // The meaning of the data items is user-defined.
    shared static double[] userPad;
    shared static int user_pad_length = 0;
    //
    shared static bool include_quality = false; // if true, we include quality in the solution file
    //
    shared static int nBlocks = 0; // Number of blocks in the overall simulation (nFluidBlocks + nSolidBlocks).
    shared static int nFluidBlocks = 0; // Number of fluid blocks in the overall simulation.
    shared static int nSolidBlocks = 0; // Number of solid blocks in the overall simulation.
    shared static int dimensions = 2; // default is 2, other valid option is 3
    // Approximating the centroids of 3D cells via a simple averaging of vertex positions
    // has shown to be more robust when importing an unstructurd grid and appears to provide
    // a better distribution of points for the least-squares gradient estimation [KAD 2022-08-03]
    // So, by default, we use the average points.
    shared static bool true_centroids = false;
    shared static bool axisymmetric = false;
    shared static bool gravity_non_zero = false;
    shared static double gravity_x = 0.0, gravity_y = 0.0, gravity_z = 0.0; // N/kg
    static ConservedQuantitiesIndices* cqi;
    shared static int nFluidBlockArrays = 0;
    //
    // Parameters controlling update
    shared static GasdynamicUpdate gasdynamic_update_scheme = GasdynamicUpdate.pc;
    shared static bool eval_udf_source_terms_at_each_stage = false;
    shared static size_t n_flow_time_levels = 3;
    shared static bool residual_smoothing = false;
    shared static double residual_smoothing_weight = 0.2;
    shared static ResidualSmoothingType residual_smoothing_type = ResidualSmoothingType.explicit;
    shared static int residual_smoothing_iterations = 2;
    shared static bool with_local_time_stepping = false;
    shared static int local_time_stepping_limit_factor = 10000;
    shared static bool with_super_time_stepping = false;
    shared static bool with_super_time_stepping_flexible_stages = false;
    shared static int max_attempts_for_step = 3; // 3 for resilience, 1 for early fail
    shared static double perturbation_for_real_differences = 1.0e-6; // When forming the Jacobian with real finite differences.
    //
    // Parameter controlling Strang-splitting mode when simulating reacting flows
    shared static StrangSplittingMode strangSplitting = StrangSplittingMode.full_T_full_R;
    //
    // Parameters controlling solid domain update
    shared static SolidDomainCoupling coupling_with_solid_domains = SolidDomainCoupling.tight;
    shared static SolidDomainLooseUpdateOptions sdluOptions;
    shared static bool solid_has_isotropic_properties = true;
    shared static bool solid_has_homogeneous_properties = true;
    shared static bool solid_domain_augmented_deriv_avg = true;
    shared static bool fluid_solid_bc_use_heat_transfer_coeff = false;
    shared static double solid_domain_cfl = 0.85;
    //
    // Parameters related to possible motion of the grid.
    shared static grid_motion = GridMotion.none;
    shared static bool write_vertex_velocities = false;
    shared static string udf_grid_motion_file; // empty to start
    static lua_State* master_lua_State; // null to start
    shared static size_t n_grid_time_levels = 1;
    //
    // The number of ghost-cell layers is adjustable for structured-grid blocks.
    // Ghost-cell layers surround the active cells of a block.
    // For the high-order reconstruction right to the edge on
    // structured-grids, we will need a minimum number of ghost cells.
    // 2 is the classic number of layers (since 1991), allowing piece-wise-parabolic reconstruction.
    // 3 will allow higher-order reconstruction for Lachlan's work.
    shared static int n_ghost_cell_layers = 2;
    // The number of ghost-cell layers for unstructured-grid blocks is always assumed to be 1.
    //
    // Shock-fitting
    //
    // The amount of time by which to delay any grid movement for shock fitting.
    // We'll often be doing shock-fitting of a strong bow shock over a blunt body.
    // To get the simulation started, we may operate with a fixed grid,
    // with the flow solver in shock-capturing mode.
    // Once the bow shock has formed, we then allow the inflow boundary
    // to be moved toward the captured shock.
    shared static double shock_fitting_delay = 0.0;
    // We may want to override reconstruction at the shock fitting inflow boundary.
    shared static bool shock_fitting_allow_flow_reconstruction = true;
    // Scaling factor applied to vertices in shock fitting simulations for stability.
    shared static double shock_fitting_scale_factor = 0.5;
    // Shock-fitting corrugation/kink filter.
    shared static double shock_fitting_filter_velocity_scale = 0.0; // default is none
    shared static bool shock_fitting_assume_symmetry_at_first_point = false;
    //
    // Some of the user-defined functionality depends on having access to all blocks
    // from a single thread.  For safety, in those cases, do not use parallel loops.
    shared static bool apply_bcs_in_parallel = true;
    //
    // When decoding the array of conserved quantities,
    // the temperature or the density may try to go negative.
    //
    // If the density is OK but the update fails to find a valid temperature,
    // it is possible that the internal energy is erroneously small and
    // it may be reasonable to ignore the failure, resetting a low temperature.
    shared static bool ignore_low_T_thermo_update_failure = true;
    shared static double suggested_low_T_value = 200.0;
    //
    /// If the cell still has invalid data and adjust_invalid_cell_data == true,
    /// the cell data are adjusted to make them reasonable.
    shared static bool adjust_invalid_cell_data = false;
    shared static bool report_invalid_cells = true;
    // The maximum number of bad cells (per block) that will be tolerated
    // without throwing an exception.
    shared static int max_invalid_cells = 0;
    shared static FlowStateLimits flowstate_limits;
    // The velocities in the user-specified FlowStates are assumed to be in the
    // non-rotating frame of reference.
    shared static bool user_specified_velocities_are_in_non_rotating_frame = true;
    //
    // Convective flux calculation can be via either:
    // (1) low-order, one-dimensional flux calculators with local flow-field
    //     reconstruction done as a preprocessing step, or
    // (2) a high-order flux calculator.
    shared static bool high_order_flux_calculator = false;
    //
    // Parameters controlling the low-order flux calculators and the
    // preprocessing reconstruction/interpolation procedure.
    // This is the "classic" arrangement for Eilmer so there are a lot of parameters.
    //
    // Default flux calculator is the adaptive mix of (diffusive) Hanel and AUSMDV.
    shared static FluxCalculator flux_calculator = FluxCalculator.adaptive_hanel_ausmdv;
    //
    // We use interpolation and reconstruction to mean the same.
    // It comes in a number of flavours:
    // 1. Low order uses just the cell-centre data as left- and right-
    //    flow properties in the convective flux calculation.
    // 2. adds a correction term to the cell-centre values, to approach something like
    //    a piecewise-quadratic interpolation between the cell centres for structured-grids
    //    and a linear model across a cloud of cell centres for unstructured-grids.
    // 3. high-order reconstruction on structured-grids using Lagrangian interpolation
    //    across a 6-cell stencil.  Must be used with ghost-cell-based boundary conditions.
    shared static int interpolation_order = 2;
    // We have the option to start a calculation without high-order reconstruction
    // and later activate it, presumably once the difficult flow has passed.
    shared static double interpolation_delay = 0.0;
    // For the implicit-update schemes (backward_euler, implicit_rk1), a quick and dirty
    // calculation of the coefficients in the sensitivity matrix will suppress reconstruction
    // for each call of evaluation of R(U).  This leads to some inconsistency, so the user
    // may wish to pay the computational expense and do the reconstruction while assembling
    // the matrix.
    shared static bool allow_interpolation_for_sensitivity_matrix = false;
    // We may elect to suppress reconstruction in particular zones always.
    static BlockZone[] suppress_reconstruction_zones;
    // For structured-grid blocks in axisymmetric simulations,
    // we may wish to suppress the reconstruction as we approach the x-axis.
    // Faces along the x-axis and the first row off the axis are flagged.
    // Faces normal to the x-axis are not flagged.
    shared static bool suppress_radial_reconstruction_at_xaxis = false;
    // We will activate the shock detector if selected features need it.
    shared static bool do_shock_detect = false;
    shared static bool damped_outflow = false;
    // enforce strict usage of the shock detector, if either the interface or
    // a touching cell is marked assume need for only diffusive flux
    shared static bool strict_shock_detector = true;
    // We might optionally want to suppress reconstruction at faces at
    // shocks and boundaries. Presently, we need to opt-in to these features.
    shared static bool suppress_reconstruction_at_shocks = false;
    shared static bool suppress_reconstruction_at_boundaries = false;
    // Default flow-data reconstruction includes interpolation of density
    // and internal energy.  Other options for the thermodunamic properties
    // to be interpolated are pressure+temperature, density+temperature and
    // density+pressure.
    shared static InterpolateOption thermo_interpolator = InterpolateOption.rhou;
    shared static bool allow_reconstruction_for_species = true;
    shared static bool allow_reconstruction_for_energy_modes = true;
    shared static bool allow_reconstruction_for_turbulent_variables = true;
    shared static bool apply_limiter = true;
    shared static bool extrema_clipping = true;
    shared static double epsilon_van_albada = 1e-12;
    shared static bool apply_heuristic_pressure_based_limiting = false;
    shared static bool interpolate_in_local_frame = true; // only for structured-grid
    // The unstructured solver has a selection of limiters available
    shared static UnstructuredLimiter unstructured_limiter = UnstructuredLimiter.venkat;
    shared static int freeze_limiter_on_step = 1_000_000_000;
    shared static bool frozen_limiter = false;
    // Allow the AUSMDV entropy fix to be switched off
    // Note: this is an experimental feature that will probably be removed in a later revision [KAD 20-12-2021]
    shared static bool apply_entropy_fix = true;
    //
    // Switch for enforcing strict positivity on the species densities in fvcell's decode routine
    shared static bool enforce_species_density_positivity = false;
    // A switch for calling scale_mass_fractions in onedinterp.d
    shared static bool scale_species_after_reconstruction = true;
    //
    // Allow the least-squares cloud of points (used to compute a cell-center gradient for
    // reconstruction in the unstructured solver) to grow.
    shared static bool use_extended_stencil = false;
    shared static double smooth_limiter_coeff = 0.3;
    // There are another couple of reconstruction-control parameters
    // further down in the viscous effects parameters.
    //
    shared static int nsteps_of_chemistry_ramp = -1;
    //
    // Set the tolerance to shear when applying the adaptive flux calculator.
    // We don't want EFM to be applied in situations of significant shear.
    // The shear value is computed as the tangential-velocity difference across an interface
    // normalised by the local sound speed.
    shared static double shear_tolerance = 0.20;
    //
    // Reference free-stream Mach number, for use in the ausm_plus_up flux calculator.
    // Choose a value for M_inf that is good for low Mach numbers.
    // To be strictly correct, we should set this at run time
    // if an M_inf value is easily defined.
    shared static double M_inf = 0.01;
    //
    // Set the tolerance in relative velocity change for the shock detector.
    // This value is expected to be a negative number (for compression)
    // and not too large in magnitude.
    // We have been using a value of -0.05 for years, based on some
    // early experiments with the sod and cone20 test cases, however,
    // the values may need to be tuned for other cases, especially where
    // viscous effects are important.
    shared static double compression_tolerance = -0.30;
    shared static ShockDetector shock_detector = ShockDetector.PJ;
    //
    // How many iterations to perform shock detector averaging
    shared static int shock_detector_smoothing = 0;
    //
    // Shock detector freezing
    shared static bool frozen_shock_detector = false;
    shared static int shock_detector_freeze_step = 1_000_000_000;
    //
    // With this flag on, the energy equation is modified such that
    // an artificial compressibility form of equations is solved.
    shared static bool artificial_compressibility = false;
    shared static double ac_alpha = 0.1;
    //
    // With this flag on, finite-rate evolution of the vibrational energies
    // (and in turn the total energy) is computed.
    shared static bool thermal_energy_exchange = false;
    //
    // For including radiation energy exchange.
    //
    shared static bool radiation = false;
    shared static int radiation_update_frequency; // = 1 for every time-step
    shared static bool radiation_scaling = false;
    shared static bool halt_on_large_flow_change = false;
    // Set to true to halt simulation when any monitor point sees a large flow change.
    shared static double tolerance_in_T;   // Temperature change for the flow change.
    //
    shared static bool electric_field_work;
    //
    // For Daryl Bond and Vince Wheatley's Single-fluid MHD additions.
    //
    shared static bool MHD = false;
    shared static bool MHD_static_field = false; // A value of false allows the internal update.
    shared static bool MHD_resistive = false;
    //
    // Lachlan Whyborn's Divergence cleaning to go with MHD.
    shared static bool divergence_cleaning = false;
    shared static double c_h = 0.0;
    shared static double divB_damping_length = 1.0;
    // Activate the electric field solver by Nick Gibbons
    shared static int electric_field_count = 1000000000;
    shared static bool solve_electric_field = false;
    shared static string field_conductivity_model = "none";

    // Parameters controlling viscous/molecular transport
    //
    shared static bool viscous = false;
    // If true, viscous effects are included in the gas-dynamic update.
    shared static bool use_viscosity_from_cells = false;
    // Proper treatment of viscous effects at the walls implies using the viscous
    // transport coefficients evaluated at the wall temperature, however,
    // Eilmer3 was set up to use the cell-average values of transport coefficients
    // because a full gas state was not available at the cell interfaces.
    // Setting use_viscosity_from_cells to true should give Eilmer3-like behaviour.
    shared static SpatialDerivCalc spatial_deriv_calc = SpatialDerivCalc.least_squares;
    shared static SpatialDerivLocn spatial_deriv_locn = SpatialDerivLocn.cells;
    // 2021-07-14, PJ, Change default settings for gradient calculations
    // for viscous fluxes. The old, Eilmer3-compatible settings were
    // SpatialDerivCalc.divergence and SpatialDerivLocn.vertices.
    shared static bool include_ghost_cells_in_spatial_deriv_clouds = true;
    shared static bool upwind_vertex_gradients = true;
    // We may elect to suppress the calculation of gradients in particular zones.
    static BlockZone[] suppress_viscous_stresses_zones;
    //
    // save the gradients used in the viscous calculations to file
    shared static bool save_convective_gradients = false;
    // save the cell-centered limiter values used in the flowstate reconstruction to file
    shared static bool save_viscous_gradients = false;
    // save the gradients used in the unstructured reconstruction calculations to file
    shared static bool save_limiter_values = false;
    // save the cell residual values to file
    shared static bool save_residual_values = false;
    // save the cell timestep values to file
    shared static bool save_timestep_values = false;
    // how often in space to save the state to disk
    shared static int nic_write = 1;
    shared static int njc_write = 1;
    shared static int nkc_write = 1;
    //
    // A factor to scale the viscosity in order to achieve a soft start.
    // The soft-start for viscous effects may be handy for impulsively-started flows.
    // A value of 1.0 means that the viscous effects are fully applied.
    shared static double viscous_factor = 1.0;
    // The amount by which to increment the viscous factor during soft-start.
    shared static double viscous_factor_increment = 0.01;
    shared static double viscous_delay = 0.0;
    //
    // When things go wrong with the least-squares estimates of flow gradients
    // it might make the flow-solver updates more stable if the consequent shear
    // stresses were limited to physically-more-reasonable values.
    // Also, for leading edged of sharp plates, the continuum model leads to
    // unreasonably large shear stresses that cannot be true in the physical world.
    shared static double shear_stress_relative_limit = 1.0;
    shared static bool apply_shear_stress_relative_limit = false;
    //
    shared static MassDiffusionModel mass_diffusion_model = MassDiffusionModel.none;
    static MassDiffusion massDiffusion;
    shared static string diffusion_coefficient_type = "none";
    shared static double lewis_number = 1.0;
    //
    shared static string turbulence_model_name = "none";
    shared static double turbulence_prandtl_number = 0.89;
    shared static double turbulence_schmidt_number = 0.75;
    shared static double max_mu_t_factor = 3000.0;
    shared static double transient_mu_t_factor = 1.0;
    shared static double freestream_turbulent_intensity = 0.01;
    static TurbulenceModel turb_model;
    static BlockZone[] turbulent_zones;
    //
    // Indicate presence of user-defined source terms
    shared static string udf_source_terms_file; // empty to start
    shared static bool udf_source_terms = false;
    //
    // Parameters controlling thermochemistry
    //
    // Turning on the reactions activates the chemical update function calls.
    // Chemical equilibrium simulations (via Look-Up Table) does not use this
    // chemical update function call.
    shared static ChemistryUpdateMode chemistry_update = ChemistryUpdateMode.split; // or integral
    shared static bool reacting = false;
    shared static string reactions_file; // empty to start
    shared static double reaction_time_delay = 0.0;
    static Schedule!double reaction_fraction_schedule;
    shared static double T_frozen = 300.0; // temperature (in K) below which reactions are frozen
    shared static double T_frozen_energy = 300.0; // temperature (in K) below which energy exchanges are skipped
    static BlockZone[] reaction_zones;
    shared static double ignition_time_start = 0.0;
    shared static double ignition_time_stop = 0.0;
    static IgnitionZone[] ignition_zones;
    shared static bool ignition_zone_active = false;
    shared static string energy_exchange_file; // empty to start
    // JJ Hoste's Turbulence-Chemistry Interaction model
    shared static TCIModel tci_model = TCIModel.none;
    shared static bool radiation_energy_dump_allowed = false;
    shared static double radiation_energy_dump_temperature_limit = 30000.0;
    //
    // Parameters controlling other simulation options
    //
    shared static int max_step = 100;      // iteration limit
    shared static int t_level;             // time level within update
    shared static int halt_now = 0;        // flag for premature halt
    shared static int write_loads_at_step = -1; // flag for premature writing of loads files
    shared static int print_count = 20; // Number of steps between writing messages to console.
    shared static int control_count = 10; // Number of steps between rereading .control file.
    //
    shared static int verbosity_level = 1;
    // Messages have a hierarchy:
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose boundary condition messages
    //
    shared static bool report_residuals; // indicate if residuals are computed and reported
    //                                   // to a file for time-integrated simulations
    //
    shared static double start_time = 0.0; // Initial solution time, in seconds.
    shared static double max_time = 1.0e-3; // final solution time, in seconds, set by user
    shared static double dt_init = 1.0e-3; // initial time step, set by user
    shared static double dt_max = 1.0e-3; // Maximum allowable time-step, after all other considerations.
    //
    // In the input script, the user may specify a single cfl_value of a schedule of values.
    // Either spec will be written into the JSON .config file as a pair of tables.
    // These will specify the target CFL number, interpolated from (time, value) pairs.
    shared static double cfl_value = 0.5;
    static Schedule!double cfl_schedule;
    shared static cfl_scale_factor = 1.0; // You may edit this factor in the .control file to modulate cfl.
    shared static bool stringent_cfl = false;
    // If true, assume the worst with respect to cell geometry and wave speed.
    shared static double viscous_signal_factor = 1.0; // can reduce the viscous influence in CFL condition
    shared static double turbulent_signal_factor = 1.0; // can reduce the turbulent omega influence in CFL condition
    shared static int cfl_count = 10;  // steps between checking time step size
    shared static bool fixed_time_step = false; // set true to fix dt_allow
    //
    shared static double dt_plot = 1.0e-3; // interval for writing soln
    shared static int write_flow_solution_at_step = -1; // flag for premature writing of flow solution files
    shared static double dt_history = 1.0e-3; // interval for writing sample
    shared static double dt_loads = 1.0e-3; // interval for writing loads on boundary groups
    // For controlling the writing of snapshots
    shared static int snapshotCount = 1_000_000_000; // Set to something very large so that default behaviour
    //                                               // does not attempt to write snapshots.
    shared static int nTotalSnapshots = 0; // By default, do not write any snapshots
    // The following initialization preserves old behaviour
    // where only one group called loads was expected.
    shared static string boundary_groups_for_loads = "loads";
    shared static string[] group_names_for_loads = ["loads"];
    shared static bool write_loads = false;
    shared static bool compute_run_time_loads = false;
    shared static int run_time_loads_count = 100;
    shared static Tuple!(size_t, size_t)[] hcells;
    shared static Tuple!(size_t, size_t)[] solid_hcells;
    //
    shared static double energy_residual;      // to be monitored for steady state
    shared static Vector3 energy_residual_loc; // location of largest value
    shared static double mass_residual;
    shared static Vector3 mass_residual_loc;
    //
    // Parameters related to special block initialisation
    shared static bool diffuseWallBCsOnInit = false;
    shared static int nInitPasses = 30;
    shared static double initTWall = -1.0; // negative indicates that initTWall is NOT used.
    //
    // Block-marching parameters
    shared static bool block_marching = false;
    shared static int nib = 1;
    shared static int njb = 1;
    shared static int nkb = 1;
    shared static bool propagate_inflow_data = false;
    shared static bool save_intermediate_results = false;
    //
    // Parameters related to the solid domain and conjugate coupling
    shared static bool udfSolidSourceTerms = false;
    shared static string udfSolidSourceTermsFile = "dummy-solid-source-terms.txt";
    //
    // Parameters for the on-the-go DFT
    shared static bool do_temporal_DFT = false;
    shared static int DFT_n_modes = 5;
    shared static int DFT_step_interval = 10;
    //
    shared static bool do_flow_average = false;
    //
    // Parameters related to the gpu chemistry mode
    version (gpu_chem) {
        static GPUChem gpuChem;
    }
    version (nk_accelerator) {
        static SteadyStateSolverOptions sssOptions;
    }
    version (shape_sensitivity) {
        static ShapeSensitivityCalculatorOptions sscOptions;
    }
    static void finalize()
    {
        lua_close(master_lua_State);
    }
    shared static string[] flow_variable_list;
} // end class GlobalConfig

// A place to store configuration parameters that are in memory that can be
// dedicated to a thread.   We will want to access them frequently
// at the lower levels of the code without having to guard them with memory barriers.
class LocalConfig {
public:
    bool in_mpi_context;
    int universe_blk_id;
    string grid_format;
    string flow_format;
    bool new_flow_format;
    //
    int dimensions;
    bool true_centroids;
    bool axisymmetric;
    bool gravity_non_zero;
    Vector3 gravity;
    ConservedQuantitiesIndices* cqi;
    GasdynamicUpdate gasdynamic_update_scheme;
    size_t n_flow_time_levels;
    bool residual_smoothing;
    bool with_local_time_stepping;
    int local_time_stepping_limit_factor;
    double perturbation_for_real_differences;
    bool with_super_time_stepping;
    bool with_super_time_stepping_flexible_stages;
    GridMotion grid_motion;
    string udf_grid_motion_file;
    size_t n_grid_time_levels;
    //
    bool shock_fitting_allow_flow_reconstruction;
    double shock_fitting_scale_factor;
    //
    bool solid_has_isotropic_properties;
    bool solid_has_homogeneous_properties;
    bool solid_domain_augmented_deriv_avg;
    bool fluid_solid_bc_use_heat_transfer_coeff;
    //
    bool ignore_low_T_thermo_update_failure;
    double suggested_low_T_value;
    bool adjust_invalid_cell_data;
    bool report_invalid_cells;
    FlowStateLimits flowstate_limits;
    bool user_specified_velocities_are_in_non_rotating_frame;
    //
    bool high_order_flux_calculator;
    FluxCalculator flux_calculator;
    int interpolation_order;
    double interpolation_delay;
    BlockZone[] suppress_reconstruction_zones;
    bool suppress_radial_reconstruction_at_xaxis;
    bool suppress_reconstruction_at_shocks;
    bool suppress_reconstruction_at_boundaries;
    InterpolateOption thermo_interpolator;
    bool allow_reconstruction_for_species;
    bool allow_reconstruction_for_energy_modes;
    bool allow_reconstruction_for_turbulent_variables;
    bool apply_limiter;
    double epsilon_van_albada;
    bool extrema_clipping;
    bool apply_heuristic_pressure_based_limiting;
    bool interpolate_in_local_frame;
    bool apply_entropy_fix;
    bool enforce_species_density_positivity;
    bool scale_species_after_reconstruction;
    UnstructuredLimiter unstructured_limiter;
    int freeze_limiter_on_step;
    bool use_extended_stencil;
    double smooth_limiter_coeff;
    int nsteps_of_chemistry_ramp;
    double shear_tolerance;
    double M_inf;
    double compression_tolerance;
    int shock_detector_smoothing;
    bool frozen_shock_detector;
    int shock_detector_freeze_step;
    ShockDetector shock_detector;
    bool do_shock_detect;
    bool damped_outflow;
    bool strict_shock_detector;
    bool artificial_compressibility;
    double ac_alpha;
    //
    bool radiation;
    bool electric_field_work;
    bool MHD;
    bool MHD_static_field;
    bool MHD_resistive;
    bool divergence_cleaning;
    double c_h;
    double divB_damping_length;
    int electric_field_count;
    bool solve_electric_field;
    string field_conductivity_model;
    //
    bool viscous;
    bool use_viscosity_from_cells;
    bool save_convective_gradients;
    bool save_viscous_gradients;
    bool save_limiter_values;
    bool save_residual_values;
    bool save_timestep_values;
    //
    int nic_write;
    int njc_write;
    int nkc_write;
    //
    SpatialDerivCalc spatial_deriv_calc;
    SpatialDerivLocn spatial_deriv_locn;
    bool include_ghost_cells_in_spatial_deriv_clouds;
    bool upwind_vertex_gradients;
    BlockZone[] suppress_viscous_stresses_zones;
    double viscous_factor;
    double shear_stress_relative_limit;
    bool apply_shear_stress_relative_limit;
    MassDiffusionModel mass_diffusion_model;
    MassDiffusion massDiffusion;
    string diffusion_coefficient_type;
    double lewis_number;
    //
    bool stringent_cfl;
    double viscous_signal_factor;
    double turbulent_signal_factor;
    //
    string turbulence_model_name;
    double turbulence_prandtl_number;
    double turbulence_schmidt_number;
    double max_mu_t_factor;
    double transient_mu_t_factor;
    double freestream_turbulent_intensity;
    TurbulenceModel turb_model;
    BlockZone[] turbulent_zones;
    //
    bool udf_source_terms;
    //
    ChemistryUpdateMode chemistry_update;
    bool reacting;
    double reaction_time_delay;
    double T_frozen;
    double T_frozen_energy;
    BlockZone[] reaction_zones;
    TCIModel tci_model;
    bool radiation_energy_dump_allowed;
    double radiation_energy_dump_temperature_limit;
    //
    double ignition_time_start;
    double ignition_time_stop;
    IgnitionZone[] ignition_zones;
    bool ignition_zone_active;
    //
    GasModel gmodel;
    uint n_species;
    uint n_heavy;
    uint n_modes;
    bool sticky_electrons;
    bool include_quality;
    ThermochemicalReactor thermochemUpdate;
    //
    int verbosity_level;
    //
    bool do_temporal_DFT;
    int DFT_n_modes;
    int DFT_step_interval;
    //
    bool do_flow_average;
    //
    version (nk_accelerator) {
        SteadyStateSolverOptions sssOptions;
    }
    version (shape_sensitivity) {
        ShapeSensitivityCalculatorOptions sscOptions;
    }
    string[] flow_variable_list;
    //
    this(int universe_blk_id)
    {
        alias cfg = GlobalConfig;
        in_mpi_context = cfg.in_mpi_context;
        this.universe_blk_id = universe_blk_id;
        grid_format = cfg.grid_format;
        flow_format = cfg.flow_format;
        new_flow_format = cfg.new_flow_format;
        dimensions = cfg.dimensions;
        true_centroids = cfg.true_centroids;
        axisymmetric = cfg.axisymmetric;
        gravity_non_zero = cfg.gravity_non_zero;
        if (gravity_non_zero) { gravity.set(cfg.gravity_x, cfg.gravity_y, cfg.gravity_z); }
        // Sometimes this constructor is called and GlobalConfig will not have been completely initialized.
        // Presumably these times, not all of the GlobalConfig is needed, so press on doing what can be done.
        if (cfg.cqi) { cqi = new ConservedQuantitiesIndices(*(cfg.cqi)); }
        gasdynamic_update_scheme = cfg.gasdynamic_update_scheme;
        n_flow_time_levels = cfg.n_flow_time_levels;
        residual_smoothing = cfg.residual_smoothing;
        with_local_time_stepping = cfg.with_local_time_stepping;
        local_time_stepping_limit_factor = cfg.local_time_stepping_limit_factor;
        perturbation_for_real_differences = cfg.perturbation_for_real_differences;
	with_super_time_stepping = cfg.with_super_time_stepping;
        with_super_time_stepping_flexible_stages = cfg.with_super_time_stepping_flexible_stages;
        grid_motion = cfg.grid_motion;
        udf_grid_motion_file = cfg.udf_grid_motion_file;
        n_grid_time_levels = cfg.n_grid_time_levels;
        //
        shock_fitting_allow_flow_reconstruction = cfg.shock_fitting_allow_flow_reconstruction;
        shock_fitting_scale_factor = cfg.shock_fitting_scale_factor;
        //
        solid_has_isotropic_properties = cfg.solid_has_isotropic_properties;
        solid_has_homogeneous_properties = cfg.solid_has_homogeneous_properties;
        solid_domain_augmented_deriv_avg = cfg.solid_domain_augmented_deriv_avg;
        fluid_solid_bc_use_heat_transfer_coeff = cfg.fluid_solid_bc_use_heat_transfer_coeff;
        //
        ignore_low_T_thermo_update_failure = cfg.ignore_low_T_thermo_update_failure;
        suggested_low_T_value = cfg.suggested_low_T_value;
        adjust_invalid_cell_data = cfg.adjust_invalid_cell_data;
        report_invalid_cells = cfg.report_invalid_cells;
        flowstate_limits = cfg.flowstate_limits;
        user_specified_velocities_are_in_non_rotating_frame = cfg.user_specified_velocities_are_in_non_rotating_frame;
        //
        high_order_flux_calculator = cfg.high_order_flux_calculator;
        flux_calculator = cfg.flux_calculator;
        interpolation_order = cfg.interpolation_order;
        interpolation_delay = cfg.interpolation_delay;
        foreach (bz; cfg.suppress_reconstruction_zones) {
            suppress_reconstruction_zones ~= new BlockZone(bz);
        }
        suppress_radial_reconstruction_at_xaxis = cfg.suppress_radial_reconstruction_at_xaxis;
        suppress_reconstruction_at_shocks = cfg.suppress_reconstruction_at_shocks;
        suppress_reconstruction_at_boundaries = cfg.suppress_reconstruction_at_boundaries;
        thermo_interpolator = cfg.thermo_interpolator;
        allow_reconstruction_for_species = cfg.allow_reconstruction_for_species;
        allow_reconstruction_for_energy_modes = cfg.allow_reconstruction_for_energy_modes;
        allow_reconstruction_for_turbulent_variables = cfg.allow_reconstruction_for_turbulent_variables;
        apply_limiter = cfg.apply_limiter;
        epsilon_van_albada = cfg.epsilon_van_albada;
        extrema_clipping = cfg.extrema_clipping;
        apply_heuristic_pressure_based_limiting = cfg.apply_heuristic_pressure_based_limiting;
        interpolate_in_local_frame = cfg.interpolate_in_local_frame;
        apply_entropy_fix = cfg.apply_entropy_fix;
        enforce_species_density_positivity = cfg.enforce_species_density_positivity;
        scale_species_after_reconstruction = cfg.scale_species_after_reconstruction;
        unstructured_limiter = cfg.unstructured_limiter;
        freeze_limiter_on_step = cfg.freeze_limiter_on_step;
        use_extended_stencil = cfg.use_extended_stencil;
        smooth_limiter_coeff = cfg.smooth_limiter_coeff;
        nsteps_of_chemistry_ramp = cfg.nsteps_of_chemistry_ramp;
        shear_tolerance = cfg.shear_tolerance;
        M_inf = cfg.M_inf;
        compression_tolerance = cfg.compression_tolerance;
        shock_detector = cfg.shock_detector;
        shock_detector_smoothing = cfg.shock_detector_smoothing;
        frozen_shock_detector = cfg.frozen_shock_detector;
        shock_detector_freeze_step = cfg.shock_detector_freeze_step;
        do_shock_detect = cfg.do_shock_detect;
        damped_outflow = cfg.damped_outflow;
        strict_shock_detector = cfg.strict_shock_detector;
        //
        artificial_compressibility = cfg.artificial_compressibility;
        ac_alpha = cfg.ac_alpha;
        //
        radiation = cfg.radiation;
        electric_field_work = cfg.electric_field_work;
        MHD = cfg.MHD;
        MHD_static_field = cfg.MHD_static_field;
        MHD_resistive = cfg.MHD_resistive;
        divergence_cleaning = cfg.divergence_cleaning;
        c_h = cfg.c_h;
        divB_damping_length = cfg.divB_damping_length;
        electric_field_count = cfg.electric_field_count;
        solve_electric_field = cfg.solve_electric_field;
        field_conductivity_model = cfg.field_conductivity_model;
        //
        viscous = cfg.viscous;
        use_viscosity_from_cells = cfg.use_viscosity_from_cells;
        spatial_deriv_calc = cfg.spatial_deriv_calc;
        spatial_deriv_locn = cfg.spatial_deriv_locn;
        include_ghost_cells_in_spatial_deriv_clouds =
            cfg.include_ghost_cells_in_spatial_deriv_clouds;
        upwind_vertex_gradients = cfg.upwind_vertex_gradients;
        foreach (bz; cfg.suppress_viscous_stresses_zones) {
            suppress_viscous_stresses_zones ~= new BlockZone(bz);
        }
        save_convective_gradients = cfg.save_convective_gradients;
        save_viscous_gradients = cfg.save_viscous_gradients;
        save_limiter_values = cfg.save_limiter_values;
        save_residual_values = cfg.save_residual_values;
        save_timestep_values = cfg.save_timestep_values;
        nic_write = cfg.nic_write;
        njc_write = cfg.njc_write;
        nkc_write = cfg.nkc_write;
        shear_stress_relative_limit = cfg.shear_stress_relative_limit;
        apply_shear_stress_relative_limit = cfg.apply_shear_stress_relative_limit;
        viscous_factor = cfg.viscous_factor;
        mass_diffusion_model = cfg.mass_diffusion_model;
        diffusion_coefficient_type = cfg.diffusion_coefficient_type;
        lewis_number = cfg.lewis_number;
        //
        stringent_cfl = cfg.stringent_cfl;
        viscous_signal_factor = cfg.viscous_signal_factor;
        turbulent_signal_factor = cfg.turbulent_signal_factor;
        //
        turbulence_model_name = cfg.turbulence_model_name;
        turbulence_prandtl_number = cfg.turbulence_prandtl_number;
        turbulence_schmidt_number = cfg.turbulence_schmidt_number;
        max_mu_t_factor = cfg.max_mu_t_factor;
        transient_mu_t_factor = cfg.transient_mu_t_factor;
        freestream_turbulent_intensity = cfg.freestream_turbulent_intensity;
        turb_model = cfg.turb_model.dup;
        foreach (bz; cfg.turbulent_zones) { turbulent_zones ~= new BlockZone(bz); }
        //
        udf_source_terms = cfg.udf_source_terms;
        //
        chemistry_update = cfg.chemistry_update;
        reacting = cfg.reacting;
        reaction_time_delay = cfg.reaction_time_delay;
        T_frozen = cfg.T_frozen;
        T_frozen_energy = cfg.T_frozen_energy;
        foreach (rz; cfg.reaction_zones) { reaction_zones ~= new BlockZone(rz); }
        tci_model = cfg.tci_model;
        radiation_energy_dump_allowed = cfg.radiation_energy_dump_allowed;
        radiation_energy_dump_temperature_limit = cfg.radiation_energy_dump_temperature_limit;
        //
        ignition_time_start = cfg.ignition_time_start;
        ignition_time_stop = cfg.ignition_time_stop;
        ignition_zone_active = cfg.ignition_zone_active;
        foreach (iz; cfg.ignition_zones) { ignition_zones ~= new IgnitionZone(iz); }
        //
        verbosity_level = cfg.verbosity_level;
        //
        do_temporal_DFT = cfg.do_temporal_DFT;
        DFT_n_modes = cfg.DFT_n_modes;
        DFT_step_interval = cfg.DFT_step_interval;
        //
        do_flow_average = cfg.do_flow_average;
        //
        version (nk_accelerator) { sssOptions = cfg.sssOptions; }
        version (shape_sensitivity) { sscOptions = cfg.sscOptions; }
        foreach (varName; cfg.flow_variable_list) { flow_variable_list ~= varName; }
    } // end constructor

    void init_gas_model_bits()
    {
        // Some gas models are storage intensive, so we do not want to initialize them
        // unless we really want to use them from the LocalConfig object.
        // This function should be called for local blocks, only, within an MPI context.
        alias cfg = GlobalConfig;
        gmodel = init_gas_model(cfg.gas_model_file);
        n_species = gmodel.n_species;
        n_heavy = gmodel.n_heavy;
        n_modes = gmodel.n_modes;
        if (mass_diffusion_model != MassDiffusionModel.none) {
            massDiffusion = initMassDiffusion(gmodel, cfg.diffusion_coefficient_type, mass_diffusion_model, cfg.lewis_number);
        }
        include_quality = cfg.include_quality;
        if (cfg.reacting) {
            thermochemUpdate = init_thermochemical_reactor(gmodel, cfg.reactions_file, cfg.energy_exchange_file);
        }
    }

    void update_control_parameters()
    // to be used after reading job.control file.
    {
        alias cfg = GlobalConfig;
        stringent_cfl = cfg.stringent_cfl;
        viscous_signal_factor = cfg.viscous_signal_factor;
        turbulent_signal_factor = cfg.turbulent_signal_factor;
        version(nk_accelerator) { sssOptions = cfg.sssOptions; }
    }
} // end class LocalConfig


//-------------------------------
// PART 4. Reading from JSON file
//-------------------------------

// Utility functions to condense the following code.
//
string update_string(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONstring(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_bool(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONbool(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_int(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONint(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_double(string key, string field)
{
    return "GlobalConfig."~field~" = getJSONdouble(jsonData, \""~key~"\", GlobalConfig."~field~");";
}
string update_enum(string key, string field, string enum_from_name)
{
    return "
    try {
        string name = jsonData[\""~key~"\"].str;
        GlobalConfig."~field~" = "~enum_from_name~"(name);
    } catch (JSONException e) {
        writeln(\"WARNING! Missing config variable: '"~field~"'.\");
    }";
}

// To check if we are starting up with a suitable gas model,
// we need to know about all of the gas-model modules that are in play.
import gas.ideal_gas;
import gas.cea_gas;
import gas.therm_perf_gas;
import gas.very_viscous_air;
import gas.uniform_lut;
import gas.adaptive_lut_CEA;
import gas.ideal_air_proxy;
import gas.ideal_gas_ab;
import gas.pseudo_species_gas;
import gas.two_temperature_reacting_argon;
import gas.ideal_dissociating_gas;
import gas.two_temperature_air;
import gas.two_temperature_nitrogen;
import gas.vib_specific_nitrogen;
import gas.fuel_air_mix;
import gas.equilibrium_gas;
import gas.electronically_specific_gas: ElectronicallySpecificGas;
import gas.two_temperature_gasgiant: TwoTemperatureGasGiant;

JSONValue read_config_file()
{
    alias cfg = GlobalConfig;
    if (cfg.verbosity_level > 1) writeln("Read config file.");
    JSONValue jsonData = readJSONfile("config/"~cfg.base_file_name~".config");
    set_config_for_core(jsonData);
    set_config_for_blocks(jsonData);
    // We send this config-file data back because we have not yet finished
    // the configuration activities for the blocks and their boundary conditions.
    return jsonData;
}

void set_config_for_core(JSONValue jsonData)
{
    // Now that we have parsed JSON data, proceed to update those config values.
    // Note that some of the lines below are much longer than PJ would normally tolerate.
    // The trade-off for ease of reading with one line per entry won out...
    //
    alias cfg = GlobalConfig;
    mixin(update_string("title", "title"));
    mixin(update_double("start_time", "start_time"));
    mixin(update_string("grid_format", "grid_format"));
    mixin(update_string("flow_format", "flow_format"));
    mixin(update_bool("new_flow_format", "new_flow_format"));
    mixin(update_string("gas_model_file", "gas_model_file"));
    // The gas model may have been initialized earlier, possibly by a setGasModel call
    // in the user's Lua script.
    if (!cfg.gmodel_master) {
        auto gm = init_gas_model(cfg.gas_model_file);
        // The following checks on gas model will need to be maintained
        // as new gas models are added.
        bool multiSpecies = true; // assumed
        bool multiT = false; // assumed
        if (cast(IdealGas)gm) { multiSpecies = false; }
        if (cast(CEAGas)gm) { multiSpecies = false; }
        if (cast(VeryViscousAir)gm) { multiSpecies = false; }
        if (cast(UniformLUT)gm) { multiSpecies = false; }
        if (cast(AdaptiveLUT)gm) { multiSpecies = false; }
        if (cast(TwoTemperatureReactingArgon)gm) { multiT = true; }
        if (cast(TwoTemperatureAir)gm) { multiT = true; }
        if (cast(TwoTemperatureNitrogen)gm) { multiT = true; }
        if (cast(TwoTemperatureGasGiant)gm) { multiT = true; }
        version(multi_species_gas) {
            // Any gas model will do.
        } else {
            if (multiSpecies) { throw new Error("Cannot accommodate multi-species gas model."); }
        }
        version(multi_T_gas) {
            // Any gas model will do.
        } else {
            if (multiT) { throw new Error("Cannot accommodate multi-temperature gas model."); }
        }
        cfg.gmodel_master = gm;
    }
    cfg.n_species = cfg.gmodel_master.n_species;
    cfg.n_heavy = cfg.gmodel_master.n_heavy;
    cfg.n_modes = cfg.gmodel_master.n_modes;
    mixin(update_string("udf_supervisor_file", "udf_supervisor_file"));
    mixin(update_int("user_pad_length", "user_pad_length"));
    cfg.userPad.length = cfg.user_pad_length;
    double[] default_user_pad_data;
    foreach (i; 0 .. cfg.userPad.length) { default_user_pad_data ~= 0.0; }
    auto user_pad_data = getJSONdoublearray(jsonData, "user_pad_data", default_user_pad_data);
    foreach (i; 0 .. cfg.userPad.length) {
        cfg.userPad[i] = (i < user_pad_data.length) ? user_pad_data[i] : default_user_pad_data[i];
    }
    mixin(update_bool("sticky_electrons", "sticky_electrons"));
    mixin(update_bool("include_quality", "include_quality"));
    mixin(update_int("dimensions", "dimensions"));
    mixin(update_bool("true_centroids", "true_centroids"));
    mixin(update_bool("axisymmetric", "axisymmetric"));
    double[] components_default = [0.0, 0.0, 0.0];
    double[] components = getJSONdoublearray(jsonData, "gravity", components_default);
    cfg.gravity_x = components[0]; cfg.gravity_y = components[1]; cfg.gravity_z = components[2];
    cfg.gravity_non_zero = !((cfg.gravity_x == 0.0) && (cfg.gravity_y == 0.0) && (cfg.gravity_z == 0.0));
    if (cfg.verbosity_level > 1) {
        writeln("  grid_format: ", to!string(cfg.grid_format));
        writeln("  flow_format: ", to!string(cfg.flow_format));
        writeln("  new_flow_format: ", to!string(cfg.new_flow_format));
        writeln("  title: ", to!string(cfg.title));
        writeln("  gas_model_file: ", to!string(cfg.gas_model_file));
        writeln("  udf_supervisor_file: ", to!string(cfg.udf_supervisor_file));
        writeln("  user_pad_length: ", cfg.user_pad_length);
        writeln("  user_pad_data: ", to!string(cfg.userPad));
        writeln("  sticky_electrons: ", cfg.sticky_electrons);
        writeln("  include_quality: ", cfg.include_quality);
        writeln("  dimensions: ", cfg.dimensions);
        writeln("  true_centroids: ", cfg.true_centroids);
        writeln("  axisymmetric: ", cfg.axisymmetric);
        writeln("  gravity_non_zero: ", cfg.gravity_non_zero);
        writeln(format("  gravity: [%e, %e, %e]", cfg.gravity_x, cfg.gravity_y, cfg.gravity_z));
    }
    //
    // Parameter controlling Strang splitting mode
    mixin(update_enum("strang_splitting", "strangSplitting", "strangSplittingModeFromName"));
    //
    // Parameters controlling convective update and size of storage arrays
    //
    mixin(update_enum("gasdynamic_update_scheme", "gasdynamic_update_scheme", "update_scheme_from_name"));
    version(nk_accelerator) {
        // We need temporalIntegrationMode to fill n_flow_time_levels for the Newton-Krylov solver,
        // this parameter is sitting in the steady-state solver options in the control file,
        // which at this point in the initialisation hasn't been read in yet. So we need to dip
        // into the control file and pull this information out here. The remainder of the control
        // file will be imported at a later stage.
        // TODO: should we consider reading the control file earlier or setting this variable later? KAD 23-08-2023
        JSONValue jsonCntrlData = readJSONfile("config/"~cfg.base_file_name~".control");
        auto sssOptions = jsonCntrlData["steady_state_solver_options"];
        auto ssso = &(cfg.sssOptions);
        ssso.temporalIntegrationMode = getJSONint(sssOptions, "temporal_integration_mode", ssso.temporalIntegrationMode);
        cfg.n_flow_time_levels = 3 + cfg.sssOptions.temporalIntegrationMode;
    } else {
        cfg.n_flow_time_levels = 1 + number_of_stages_for_update_scheme(cfg.gasdynamic_update_scheme);
    }
    mixin(update_bool("eval_udf_source_terms_at_each_stage", "eval_udf_source_terms_at_each_stage"));
    // The CFL schedule arrives as a pair of tables that should have at least one entry each.
    int cfl_schedule_length = getJSONint(jsonData, "cfl_schedule_length", 1);
    double[] cfl_schedule_values_default;
    double[] cfl_schedule_times_default;
    foreach (i; 0 .. cfl_schedule_length) {
        cfl_schedule_times_default ~= 0.0;
        cfl_schedule_values_default ~= 0.5;
    }
    double[] cfl_schedule_times = getJSONdoublearray(jsonData, "cfl_schedule_times", cfl_schedule_times_default);
    double[] cfl_schedule_values = getJSONdoublearray(jsonData, "cfl_schedule_values", cfl_schedule_values_default);
    cfg.cfl_schedule = new Schedule!double(cfl_schedule_times, cfl_schedule_values);
    //
    mixin(update_bool("residual_smoothing", "residual_smoothing"));
    mixin(update_bool("with_local_time_stepping", "with_local_time_stepping"));
    mixin(update_int("local_time_stepping_limit_factor", "local_time_stepping_limit_factor"));
    mixin(update_bool("with_super_time_stepping", "with_super_time_stepping"));
    mixin(update_bool("with_super_time_stepping_flexible_stages", "with_super_time_stepping_flexible_stages"));
    mixin(update_int("max_attempts_for_step", "max_attempts_for_step"));
    mixin(update_double("perturbation_for_real_differences", "perturbation_for_real_differences"));
    mixin(update_enum("grid_motion", "grid_motion", "grid_motion_from_name"));
    if (cfg.grid_motion == GridMotion.none) {
        cfg.n_grid_time_levels = 1;
    } else {
        cfg.n_grid_time_levels = 1 + number_of_stages_for_update_scheme(cfg.gasdynamic_update_scheme);
    }
    mixin(update_bool("write_vertex_velocities", "write_vertex_velocities"));
    mixin(update_string("udf_grid_motion_file", "udf_grid_motion_file"));
    // If we have user-defined grid motion, we'll need to initialise
    // the lua_State that holds the user's function. But we can't
    // do that initialisation just yet. We have to wait until all
    // of the blocks are configured since we set that information
    // as global information available to the user. Hence, you'll
    // find that step at the very end of this function.
    //
    // Shock fitting involves grid motion.
    mixin(update_double("shock_fitting_delay", "shock_fitting_delay"));
    mixin(update_bool("shock_fitting_allow_flow_reconstruction", "shock_fitting_allow_flow_reconstruction"));
    mixin(update_double("shock_fitting_scale_factor", "shock_fitting_scale_factor"));
    mixin(update_double("shock_fitting_filter_velocity_scale", "shock_fitting_filter_velocity_scale"));
    mixin(update_bool("shock_fitting_assume_symmetry_at_first_point", "shock_fitting_assume_symmetry_at_first_point"));
    //
    mixin(update_double("solid_domain_cfl", "solid_domain_cfl"));
    mixin(update_enum("coupling_with_solid_domains", "coupling_with_solid_domains", "solidDomainCouplingFromName"));
    mixin(update_bool("solid_has_isotropic_properties", "solid_has_isotropic_properties"));
    mixin(update_bool("solid_has_homogeneous_properties", "solid_has_homogeneous_properties"));
    mixin(update_bool("solid_domain_augmented_deriv_avg", "solid_domain_augmented_deriv_avg"));
    mixin(update_bool("fluid_solid_bc_use_heat_transfer_coeff", "fluid_solid_bc_use_heat_transfer_coeff"));

    // Parameters controlling convective update in detail
    //
    mixin(update_bool("apply_bcs_in_parallel", "apply_bcs_in_parallel"));
    mixin(update_double("flowstate_limits_max_velocity", "flowstate_limits.max_velocity"));
    mixin(update_double("flowstate_limits_max_tke", "flowstate_limits.max_tke"));
    mixin(update_double("flowstate_limits_min_tke", "flowstate_limits.min_tke"));
    mixin(update_double("flowstate_limits_max_temp", "flowstate_limits.max_temp"));
    mixin(update_double("flowstate_limits_min_temp", "flowstate_limits.min_temp"));
    mixin(update_double("flowstate_limits_min_pressure", "flowstate_limits.min_pressure"));
    mixin(update_bool("user_specified_velocities_are_in_non_rotating_frame", "user_specified_velocities_are_in_non_rotating_frame"));
    mixin(update_bool("ignore_low_T_thermo_update_failure", "ignore_low_T_thermo_update_failure"));
    mixin(update_double("suggested_low_T_value", "suggested_low_T_value"));
    mixin(update_bool("adjust_invalid_cell_data", "adjust_invalid_cell_data"));
    mixin(update_bool("report_invalid_cells", "report_invalid_cells"));
    mixin(update_int("max_invalid_cells", "max_invalid_cells"));
    //
    mixin(update_int("n_ghost_cell_layers", "n_ghost_cell_layers"));
    mixin(update_bool("high_order_flux_calculator", "high_order_flux_calculator"));
    mixin(update_enum("flux_calculator", "flux_calculator", "flux_calculator_from_name"));
    mixin(update_int("interpolation_order", "interpolation_order"));
    mixin(update_double("interpolation_delay", "interpolation_delay"));
    mixin(update_bool("allow_interpolation_for_sensitivity_matrix", "allow_interpolation_for_sensitivity_matrix"));
    mixin(update_bool("suppress_radial_reconstruction_at_xaxis", "suppress_radial_reconstruction_at_xaxis"));
    mixin(update_bool("suppress_reconstruction_at_shocks", "suppress_reconstruction_at_shocks"));
    mixin(update_bool("suppress_reconstruction_at_captured_shocks", "suppress_reconstruction_at_shocks")); // old name
    mixin(update_bool("suppress_reconstruction_at_boundaries", "suppress_reconstruction_at_boundaries"));
    mixin(update_enum("thermo_interpolator", "thermo_interpolator", "thermo_interpolator_from_name"));
    mixin(update_bool("allow_reconstruction_for_species", "allow_reconstruction_for_species"));
    mixin(update_bool("allow_reconstruction_for_energy_modes", "allow_reconstruction_for_energy_modes"));
    mixin(update_bool("allow_reconstruction_for_turbulent_variables", "allow_reconstruction_for_turbulent_variables"));
    mixin(update_bool("apply_limiter", "apply_limiter"));
    mixin(update_double("epsilon_van_albada", "epsilon_van_albada"));
    mixin(update_bool("extrema_clipping", "extrema_clipping"));
    mixin(update_bool("apply_heuristic_pressure_based_limiting", "apply_heuristic_pressure_based_limiting"));
    mixin(update_bool("interpolate_in_local_frame", "interpolate_in_local_frame"));
    mixin(update_bool("apply_entropy_fix", "apply_entropy_fix"));
    mixin(update_bool("enforce_species_density_positivity", "enforce_species_density_positivity"));
    mixin(update_bool("scale_species_after_reconstruction", "scale_species_after_reconstruction"));
    mixin(update_enum("unstructured_limiter", "unstructured_limiter", "unstructured_limiter_from_name"));
    mixin(update_int("freeze_limiter_on_step", "freeze_limiter_on_step"));
    mixin(update_bool("use_extended_stencil", "use_extended_stencil"));
    mixin(update_double("smooth_limiter_coeff", "smooth_limiter_coeff"));
    mixin(update_int("nsteps_of_chemistry_ramp", "nsteps_of_chemistry_ramp"));
    mixin(update_double("shear_tolerance", "shear_tolerance"));
    mixin(update_int("shock_detector_smoothing", "shock_detector_smoothing"));
    mixin(update_bool("frozen_shock_detector", "frozen_shock_detector"));
    mixin(update_int("shock_detector_freeze_step", "shock_detector_freeze_step"));
    mixin(update_double("M_inf", "M_inf"));
    mixin(update_double("compression_tolerance", "compression_tolerance"));
    mixin(update_enum("shock_detector", "shock_detector", "shock_detector_from_name"));
    mixin(update_bool("do_shock_detect", "do_shock_detect"));
    mixin(update_bool("damped_outflow", "damped_outflow"));
    mixin(update_bool("strict_shock_detector", "strict_shock_detector"));
    mixin(update_enum("flux_calculator", "flux_calculator", "flux_calculator_from_name"));
    mixin(update_bool("artificial_compressibility", "artificial_compressibility"));
    mixin(update_double("ac_alpha", "ac_alpha"));
    mixin(update_bool("radiation", "radiation"));
    mixin(update_bool("MHD", "MHD"));
    version(MHD) {
        // no problem
    } else {
        if (cfg.MHD) { throw new Error("MHD capability has not been enabled."); }
    }
    mixin(update_bool("MHD_static_field", "MHD_static_field"));
    mixin(update_bool("MHD_resistive", "MHD_resistive"));
    mixin(update_bool("divergence_cleaning", "divergence_cleaning"));
    mixin(update_double("divB_damping_length", "divB_damping_length"));
    mixin(update_int("electric_field_count", "electric_field_count"));
    mixin(update_bool("solve_electric_field", "solve_electric_field"));
    mixin(update_string("field_conductivity_model", "field_conductivity_model"));

    // Checking of constraints.
    // The following checks/overrides must happen after the relevant config elements
    // have been set.  This is the first such check.  For details, see the function below.
    configCheckPoint1();
    //
    if (cfg.verbosity_level > 1) {
        writeln("  gasdynamic_update_scheme: ", gasdynamic_update_scheme_name(cfg.gasdynamic_update_scheme));
        writeln("  eval_udf_source_terms_at_each_stage: ", cfg.eval_udf_source_terms_at_each_stage);
        writeln("  cfl_schedule: ", cfg.cfl_schedule);
        writeln("  residual_smoothing: ", cfg.residual_smoothing);
        writeln("  with_local_time_stepping: ", cfg.with_local_time_stepping);
        writeln("  local_time_stepping_limit_factor: ", cfg.local_time_stepping_limit_factor);
	writeln("  with_super_time_stepping: ", cfg.with_super_time_stepping);
        writeln("  with_super_time_stepping_flexible_stages: ", cfg.with_super_time_stepping_flexible_stages);
        writeln("  max_attempts_for_step: ", cfg.max_attempts_for_step);
        writeln("  perturbation_for_real_differences: ", cfg.perturbation_for_real_differences);
	writeln("  grid_motion: ", grid_motion_name(cfg.grid_motion));
        writeln("  write_vertex_velocities: ", cfg.write_vertex_velocities);
        writeln("  udf_grid_motion_file: ", to!string(cfg.udf_grid_motion_file));
        writeln("  shock_fitting_delay: ", cfg.shock_fitting_delay);
        writeln("  shock_fitting_allow_flow_reconstruction: ", cfg.shock_fitting_allow_flow_reconstruction);
        writeln("  shock_fitting_scale_factor: ", cfg.shock_fitting_scale_factor);
        writeln("  solid_domain_cfl: ", cfg.solid_domain_cfl);
        writeln("  coupling_with_solid_domains: ", cfg.coupling_with_solid_domains);
        writeln("  solid_has_isotropic_properties: ", cfg.solid_has_isotropic_properties);
        writeln("  solid_has_homogeneous_properties: ", cfg.solid_has_homogeneous_properties);
        writeln("  solid_domain_augmented_deriv_avg: ", cfg.solid_domain_augmented_deriv_avg);
        writeln("  fluid_solid_bc_use_heat_transfer_coeff: ", cfg.fluid_solid_bc_use_heat_transfer_coeff);
        writeln("  apply_bcs_in_parallel: ", cfg.apply_bcs_in_parallel);
        writeln("  flowstate_limits_max_velocity: ", cfg.flowstate_limits.max_velocity);
        writeln("  flowstate_limits_max_tke: ", cfg.flowstate_limits.max_tke);
        writeln("  flowstate_limits_min_tke: ", cfg.flowstate_limits.min_tke);
        writeln("  flowstate_limits_max_temp: ", cfg.flowstate_limits.max_temp);
        writeln("  flowstate_limits_min_temp: ", cfg.flowstate_limits.min_temp);
        writeln("  ignore_low_T_thermo_update_failure: ", cfg.ignore_low_T_thermo_update_failure);
        writeln("  suggested_low_T_value: ", cfg.suggested_low_T_value);
        writeln("  adjust_invalid_cell_data: ", cfg.adjust_invalid_cell_data);
        writeln("  report_invalid_cells: ", cfg.report_invalid_cells);
        writeln("  max_invalid_cells: ", cfg.max_invalid_cells);
        //
        writeln("  n_ghost_cell_layers: ", cfg.n_ghost_cell_layers);
        writeln("  high_order_flux_calculator: ", cfg.high_order_flux_calculator);
        writeln("  flux_calculator: ", flux_calculator_name(cfg.flux_calculator));
        writeln("  interpolation_order: ", cfg.interpolation_order);
        writeln("  interpolation_delay: ", cfg.interpolation_delay);
        writeln("  allow_interpolation_for_sensitivity_matrix: ", cfg.allow_interpolation_for_sensitivity_matrix);
        writeln("  suppress_radial_reconstruction_at_xaxis: ", cfg.suppress_radial_reconstruction_at_xaxis);
        writeln("  do_shock_detect: ", cfg.do_shock_detect);
        writeln("  damped_outflow: ", cfg.damped_outflow);
        writeln("  strict_shock_detector: ", cfg.strict_shock_detector);
        writeln("  suppress_reconstruction_at_shocks: ", cfg.suppress_reconstruction_at_shocks);
        writeln("  suppress_reconstruction_at_boundaries: ", cfg.suppress_reconstruction_at_boundaries);
        writeln("  thermo_interpolator: ", thermo_interpolator_name(cfg.thermo_interpolator));
        writeln("  allow_reconstruction_for_species: ", cfg.allow_reconstruction_for_species);
        writeln("  allow_reconstruction_for_energy_modes: ", cfg.allow_reconstruction_for_energy_modes);
        writeln("  allow_reconstruction_for_turbulent_variables: ", cfg.allow_reconstruction_for_turbulent_variables);
        writeln("  apply_limiter: ", cfg.apply_limiter);
        writeln("  epsilon_van_albada: ", cfg.epsilon_van_albada);
        writeln("  apply_entropy_fix: ", cfg.apply_entropy_fix);
        writeln("  enforce_species_density_positivity: ", cfg.enforce_species_density_positivity);
        writeln("  scale_species_after_reconstruction: ", cfg.scale_species_after_reconstruction);
        writeln("  unstructured_limiter: ", unstructured_limiter_name(cfg.unstructured_limiter));
        writeln("  freeze_limiter_on_step: ", cfg.freeze_limiter_on_step);
        writeln("  use_extended_stencil: ", cfg.use_extended_stencil);
        writeln("  smooth_limiter_coeff: ", cfg.smooth_limiter_coeff);
        writeln("  nsteps_of_chemistry_ramp: ", cfg.nsteps_of_chemistry_ramp);
        writeln("  extrema_clipping: ", cfg.extrema_clipping);
        writeln("  apply_heuristic_pressure_based_limiting: ", cfg.apply_heuristic_pressure_based_limiting);
        writeln("  interpolate_in_local_frame: ", cfg.interpolate_in_local_frame);
        writeln("  shear_tolerance: ", cfg.shear_tolerance);
        writeln("  M_inf: ", cfg.M_inf);
        writeln("  compression_tolerance: ", cfg.compression_tolerance);
        writeln("  shock_detector: ", shock_detector_name(cfg.shock_detector));
        writeln("  shock_detector_smoothing: ", cfg.shock_detector_smoothing);
        writeln("  frozen_shock_detector: ", cfg.frozen_shock_detector);
        writeln("  shock_detector_freeze_step: ", cfg.shock_detector_freeze_step);
        //
        writeln("  MHD: ", cfg.MHD);
        writeln("  MHD_static_field: ", cfg.MHD_static_field);
        writeln("  MHD_resistive: ", cfg.MHD_resistive);
        writeln("  divergence_cleaning: ", cfg.divergence_cleaning);
        writeln("  divB_damping_length: ", cfg.divB_damping_length);
        writeln("  electric_field_count: ", cfg.electric_field_count);
        writeln("  solve_electric_field: ", cfg.solve_electric_field);
        writeln("  field_conductivity_model: ", cfg.field_conductivity_model);
    }
    configCheckPoint2();
    //
    // Parameters controlling viscous/molecular transport
    //
    mixin(update_bool("viscous", "viscous"));
    mixin(update_bool("use_viscosity_from_cells", "use_viscosity_from_cells"));
    mixin(update_enum("spatial_deriv_calc", "spatial_deriv_calc", "spatial_deriv_calc_from_name"));
    mixin(update_enum("spatial_deriv_locn", "spatial_deriv_locn", "spatial_deriv_locn_from_name"));
    mixin(update_bool("include_ghost_cells_in_spatial_deriv_clouds", "include_ghost_cells_in_spatial_deriv_clouds"));
    mixin(update_bool("upwind_vertex_gradients", "upwind_vertex_gradients"));
    mixin(update_bool("save_convective_gradients", "save_convective_gradients"));
    mixin(update_bool("save_viscous_gradients", "save_viscous_gradients"));
    mixin(update_bool("save_limiter_values", "save_limiter_values"));
    mixin(update_bool("save_residual_values", "save_residual_values"));
    mixin(update_bool("save_timestep_values", "save_timestep_values"));
    mixin(update_int("nic_write", "nic_write"));
    mixin(update_int("njc_write", "njc_write"));
    mixin(update_int("nkc_write", "nkc_write"));
    mixin(update_double("viscous_delay", "viscous_delay"));
    mixin(update_double("viscous_factor_increment", "viscous_factor_increment"));
    mixin(update_double("shear_stress_relative_limit", "shear_stress_relative_limit"));
    mixin(update_bool("apply_shear_stress_relative_limit", "apply_shear_stress_relative_limit"));
    mixin(update_enum("mass_diffusion_model", "mass_diffusion_model", "massDiffusionModelFromName"));
    mixin(update_string("diffusion_coefficient_type", "diffusion_coefficient_type"));
    mixin(update_double("lewis_number", "lewis_number"));
    mixin(update_string("turbulence_model", "turbulence_model_name"));
    mixin(update_double("turbulence_prandtl_number", "turbulence_prandtl_number"));
    mixin(update_double("turbulence_schmidt_number", "turbulence_schmidt_number"));
    mixin(update_double("max_mu_t_factor", "max_mu_t_factor"));
    mixin(update_double("transient_mu_t_factor", "transient_mu_t_factor"));
    mixin(update_double("freestream_turbulent_intensity", "freestream_turbulent_intensity"));
    cfg.turb_model = init_turbulence_model(cfg.turbulence_model_name, jsonData);
    if (cfg.verbosity_level > 1) {
        writeln("  viscous: ", cfg.viscous);
        writeln("  use_viscosity_from_cells: ", cfg.use_viscosity_from_cells);
        writeln("  spatial_deriv_calc: ", spatial_deriv_calc_name(cfg.spatial_deriv_calc));
        writeln("  spatial_deriv_locn: ", spatial_deriv_locn_name(cfg.spatial_deriv_locn));
        writeln("  include_ghost_cells_in_spatial_deriv_clouds: ", cfg.include_ghost_cells_in_spatial_deriv_clouds);
        writeln("  upwind_vertex_gradients: ", cfg.upwind_vertex_gradients);
        writeln("  save_convective_gradients: ", cfg.save_convective_gradients);
        writeln("  save_viscous_gradients: ", cfg.save_viscous_gradients);
        writeln("  save_limiter_values: ", cfg.save_limiter_values);
        writeln("  save_residual_values: ", cfg.save_residual_values);
        writeln("  save_timestep_values: ", cfg.save_timestep_values);
        writeln("  viscous_delay: ", cfg.viscous_delay);
        writeln("  viscous_factor_increment: ", cfg.viscous_factor_increment);
        writeln("  shear_stress_relative_limit: ", cfg.shear_stress_relative_limit);
        writeln("  apply_shear_stress_relative_limit: ", cfg.apply_shear_stress_relative_limit);
        writeln("  mass_diffusion_model: ", massDiffusionModelName(cfg.mass_diffusion_model));
        writeln("  diffusion_coefficient_type: ", cfg.diffusion_coefficient_type);
        writeln("  lewis_number: ", cfg.lewis_number);
        writeln("  turbulence_model: ", cfg.turbulence_model_name);
        writeln("  turbulence_prandtl_number: ", cfg.turbulence_prandtl_number);
        writeln("  turbulence_schmidt_number: ", cfg.turbulence_schmidt_number);
        writeln("  max_mu_t_factor: ", cfg.max_mu_t_factor);
        writeln("  transient_mu_t_factor: ", cfg.transient_mu_t_factor);
        writeln("  nturb equations: ", cfg.turb_model.nturb);
        writeln("  freestream_turbulent_intensity: ", cfg.freestream_turbulent_intensity);
    }
    //
    configCheckPoint3();
    //
    if (cfg.mass_diffusion_model != MassDiffusionModel.none) {
        cfg.massDiffusion = initMassDiffusion(cfg.gmodel_master, cfg.diffusion_coefficient_type,
                                              cfg.mass_diffusion_model, cfg.lewis_number);
    }
    // User-defined source terms
    mixin(update_bool("udf_source_terms", "udf_source_terms"));
    mixin(update_string("udf_source_terms_file", "udf_source_terms_file"));
    if (cfg.verbosity_level > 1) {
        writeln("  udf_source_terms: ", cfg.udf_source_terms);
        writeln("  udf_source_terms_file: ", to!string(cfg.udf_source_terms_file));
    }
    //
    // Parameters controlling thermochemistry
    //
    mixin(update_enum("chemistry_update", "chemistry_update", "chemistry_update_mode_from_name"));
    mixin(update_bool("reacting", "reacting"));
    mixin(update_string("reactions_file", "reactions_file"));
    mixin(update_double("reaction_time_delay", "reaction_time_delay"));
    // The reaction_fractionschedule arrives as a pair of tables that should have at least one entry each.
    int rf_schedule_length = getJSONint(jsonData, "reaction_fraction_schedule_length", 1);
    double[] rf_schedule_values_default;
    double[] rf_schedule_times_default;
    foreach (i; 0 .. rf_schedule_length) {
        rf_schedule_times_default ~= 0.0;
        rf_schedule_values_default ~= 1.0;
    }
    double[] rf_schedule_times = getJSONdoublearray(jsonData, "reaction_fraction_schedule_times", rf_schedule_times_default);
    double[] rf_schedule_values = getJSONdoublearray(jsonData, "reaction_fraction_schedule_values", rf_schedule_values_default);
    cfg.reaction_fraction_schedule = new Schedule!double(rf_schedule_times, rf_schedule_values);
    mixin(update_double("T_frozen", "T_frozen"));
    mixin(update_double("T_frozen_energy", "T_frozen_energy"));
    mixin(update_enum("tci_model", "tci_model", "tci_model_from_name"));
    mixin(update_double("ignition_time_start", "ignition_time_start"));
    mixin(update_double("ignition_time_stop", "ignition_time_stop"));
    mixin(update_string("energy_exchange_file", "energy_exchange_file"));
    mixin(update_bool("radiation_energy_dump_allowed", "radiation_energy_dump_allowed"));
    mixin(update_double("radiation_energy_dump_temperature_limit", "radiation_energy_dump_temperature_limit"));
    if (cfg.verbosity_level > 1) {
        writeln("  chemistry_mode: ", chemistry_update_mode_name(cfg.chemistry_update));
        writeln("  reacting: ", cfg.reacting);
        writeln("  reactions_file: ", to!string(cfg.reactions_file));
        writeln("  reaction_time_delay: ", cfg.reaction_time_delay);
        writeln("  reaction_fraction_schedule: ", cfg.reaction_fraction_schedule);
        writeln("  T_frozen: ", cfg.T_frozen);
        writeln("  T_frozen_energy: ", cfg.T_frozen_energy);
        writeln("  tci_model: ", tci_model_name(cfg.tci_model));
        writeln("  ignition_time_start: ", cfg.ignition_time_start);
        writeln("  ignition_time_stop: ", cfg.ignition_time_start);
        writeln("  energy_exchange_file: ", cfg.energy_exchange_file);
        writeln("  radiation_energy_dump_allowed: ", cfg.radiation_energy_dump_allowed);
        writeln("  radiation_energy_dump_temperature_limit: ", cfg.radiation_energy_dump_temperature_limit);
    }
    //
    // Parameters controlling other simulation options
    //
    mixin(update_bool("diffuse_wall_bcs_on_init", "diffuseWallBCsOnInit"));
    mixin(update_int("number_init_passes", "nInitPasses"));
    mixin(update_double("wall_temperature_on_init", "initTWall"));
    mixin(update_int("control_count", "control_count"));
    mixin(update_bool("block_marching", "block_marching"));
    mixin(update_int("nib", "nib"));
    mixin(update_int("njb", "njb"));
    mixin(update_int("nkb", "nkb"));
    mixin(update_bool("propagate_inflow_data", "propagate_inflow_data"));
    mixin(update_bool("save_intermediate_results", "save_intermediate_results"));
    mixin(update_string("boundary_groups_for_loads", "boundary_groups_for_loads"));
    string[] group_names = cfg.boundary_groups_for_loads.replace(" ","").split(",");
    foreach (name; group_names) {
        if (name.length > 0) { cfg.group_names_for_loads ~= name; }
    }
    mixin(update_bool("write_loads", "write_loads"));
    mixin(update_bool("compute_run_time_loads", "compute_run_time_loads"));
    mixin(update_int("run_time_loads_count", "run_time_loads_count"));
    mixin(update_bool("do_temporal_DFT", "do_temporal_DFT"));
    mixin(update_int("DFT_n_modes", "DFT_n_modes"));
    mixin(update_int("DFT_step_interval", "DFT_step_interval"));
    mixin(update_bool("do_flow_average", "do_flow_average"));
    if (cfg.verbosity_level > 1) {
        writeln("  diffuse_wall_bcs_on_init: ", cfg.diffuseWallBCsOnInit);
        writeln("  number_init_passes: ", cfg.nInitPasses);
        writeln("  wall_temperature_on_init: ", cfg.initTWall);
        writeln("  control_count: ", cfg.control_count);
        writeln("  block_marching: ", cfg.block_marching);
        writeln("  nib: ", cfg.nib);
        writeln("  njb: ", cfg.njb);
        writeln("  nkb: ", cfg.nkb);
        writeln("  propagate_inflow_data: ", cfg.propagate_inflow_data);
        writeln("  save_intermediate_results: ", cfg.save_intermediate_results);
        writeln("  boundary_groups_for_loads: ", cfg.boundary_groups_for_loads);
        writeln("  group_names_for_loads: ", cfg.group_names_for_loads);
        writeln("  write_loads: ", cfg.write_loads);
        writeln("  do_temporal_DFT: ", cfg.do_temporal_DFT);
        writeln("  DFT_n_modes: ", cfg.DFT_n_modes);
        writeln("  DFT_step_interval: ", cfg.DFT_step_interval);
        writeln("  do_flow_average: ", cfg.do_flow_average);
    }
    //
    configCheckPoint4();
    //
    int nhcell = getJSONint(jsonData, "nhcell", 0);
    foreach (i; 0 .. nhcell) {
        string jsonKey = format("history-cell-%d", i);
        auto hcell = getJSONintarray(jsonData, jsonKey, [0, 0]);
        cfg.hcells ~= tuple(cast(size_t) hcell[0], cast(size_t) hcell[1]);
    }
    int nsolidhcell = getJSONint(jsonData, "nsolidhcell", 0);
    foreach (i; 0 .. nsolidhcell) {
        string jsonKey = format("solid-history-cell-%d", i);
        auto hcell = getJSONintarray(jsonData, jsonKey, [0, 0]);
        cfg.solid_hcells ~= tuple(cast(size_t) hcell[0], cast(size_t) hcell[1]);
    }
    int n_reaction_zones = getJSONint(jsonData, "n-reaction-zones", 0);
    foreach (i; 0 .. n_reaction_zones) {
        string jsonKey = format("reaction-zone-%d", i);
        auto zone_data = getJSONdoublearray(jsonData, jsonKey, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        Vector3 p0 = Vector3(zone_data[0], zone_data[1], zone_data[2]);
        Vector3 p1 = Vector3(zone_data[3], zone_data[4], zone_data[5]);
        cfg.reaction_zones ~= new BlockZone(p0, p1);
    }
    int n_ignition_zones = getJSONint(jsonData, "n-ignition-zones", 0);
    foreach (i; 0 .. n_ignition_zones) {
        string jsonKey = format("ignition-zone-%d", i);
        auto zone_data = getJSONdoublearray(jsonData, jsonKey, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0]);
        Vector3 p0 = Vector3(zone_data[0], zone_data[1], zone_data[2]);
        Vector3 p1 = Vector3(zone_data[3], zone_data[4], zone_data[5]);
        double Tig = zone_data[6];
        cfg.ignition_zones ~= new IgnitionZone(p0, p1, Tig);
    }
    int n_turbulent_zones = getJSONint(jsonData, "n-turbulent-zones", 0);
    foreach (i; 0 .. n_turbulent_zones) {
        string jsonKey = format("turbulent-zone-%d", i);
        auto zone_data = getJSONdoublearray(jsonData, jsonKey, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        Vector3 p0 = Vector3(zone_data[0], zone_data[1], zone_data[2]);
        Vector3 p1 = Vector3(zone_data[3], zone_data[4], zone_data[5]);
        cfg.turbulent_zones ~= new BlockZone(p0, p1);
    }
    int n_suppress_reconstruction_zones = getJSONint(jsonData, "n-suppress-reconstruction-zones", 0);
    foreach (i; 0 .. n_suppress_reconstruction_zones) {
        string jsonKey = format("suppress-reconstruction-zone-%d", i);
        auto zone_data = getJSONdoublearray(jsonData, jsonKey, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        Vector3 p0 = Vector3(zone_data[0], zone_data[1], zone_data[2]);
        Vector3 p1 = Vector3(zone_data[3], zone_data[4], zone_data[5]);
        cfg.suppress_reconstruction_zones ~= new BlockZone(p0, p1);
    }
    int n_suppress_viscous_stresses_zones = getJSONint(jsonData, "n-suppress-viscous-stresses-zones", 0);
    foreach (i; 0 .. n_suppress_viscous_stresses_zones) {
        string jsonKey = format("suppress-viscous-stresses-zone-%d", i);
        auto zone_data = getJSONdoublearray(jsonData, jsonKey, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        Vector3 p0 = Vector3(zone_data[0], zone_data[1], zone_data[2]);
        Vector3 p1 = Vector3(zone_data[3], zone_data[4], zone_data[5]);
        cfg.suppress_viscous_stresses_zones ~= new BlockZone(p0, p1);
    }
    //
    auto sdluOptions = jsonData["solid_domain_loose_update_options"];
    auto sdluo = &(cfg.sdluOptions);
    sdluo.maxNewtonIterations = getJSONint(sdluOptions, "max_newton_iterations", sdluo.maxNewtonIterations);
    sdluo.toleranceNewtonUpdate = getJSONdouble(sdluOptions, "tolerance_newton_update", sdluo.toleranceNewtonUpdate);
    sdluo.maxGMRESIterations = getJSONint(sdluOptions, "max_gmres_iterations", sdluo.maxGMRESIterations);
    sdluo.toleranceGMRESSolve = getJSONdouble(sdluOptions, "tolerance_gmres_solve", sdluo.toleranceGMRESSolve);
    sdluo.perturbationSize = getJSONdouble(sdluOptions, "perturbation_size", sdluo.perturbationSize);
    //
    version (shape_sensitivity) {
        auto sscOptions = jsonData["shape_sensitivity_calculator_options"];
        auto ssco = &(cfg.sscOptions);
        ssco.pseudotime = getJSONbool(sscOptions, "pseudotime", ssco.pseudotime);
        ssco.pseudotime_lhs_jacobian_order = getJSONint(sscOptions, "pseudotime_lhs_jacobian_order", ssco.pseudotime_lhs_jacobian_order);
        ssco.adjoint_precondition_matrix_order = getJSONint(sscOptions, "adjoint_precondition_matrix_order", ssco.adjoint_precondition_matrix_order);
        ssco.read_frozen_limiter_values_from_file = getJSONbool(sscOptions, "read_frozen_limiter_values_from_file", ssco.read_frozen_limiter_values_from_file);
        ssco.epsilon = getJSONdouble(sscOptions, "epsilon", ssco.epsilon);
        ssco.maxOuterIterations = getJSONint(sscOptions, "maxOuterIterations", ssco.maxOuterIterations);
        ssco.maxRestarts = getJSONint(sscOptions, "maxRestarts", ssco.maxRestarts);
        ssco.cfl0 = getJSONdouble(sscOptions, "cfl0", ssco.cfl0);
        ssco.eta = getJSONdouble(sscOptions, "eta", ssco.eta);
        ssco.stopOnRelativeGlobalResidual = getJSONdouble(sscOptions, "stop_on_relative_global_residual", ssco.stopOnRelativeGlobalResidual);
        ssco.tolBezierCurveFit = getJSONdouble(sscOptions, "tol_bezier_curve_fit", ssco.tolBezierCurveFit);
        ssco.maxStepsBezierCurveFit = getJSONint(sscOptions, "max_steps_bezier_curve_fit", ssco.maxStepsBezierCurveFit);
        ssco.userDefinedObjectiveFile = sscOptions["user_defined_objective_file"].str;
    }
    //
    // Enough configuration should be known, such we can build a list of variable names
    // for which data will be written into the flow data files.
    // This list needs to be built before the block-local config is copied.
    foreach (vname; build_flow_variable_list()) { cfg.flow_variable_list ~= vname; }
    //
    // We have enough information here to create the ConservedQuantitiesIndices struct.
    cfg.cqi = new ConservedQuantitiesIndices(cfg.dimensions, cfg.turb_model.nturb, cfg.MHD, cfg.n_species, cfg.n_modes);
} // end set_config_for_core()

void set_config_for_blocks(JSONValue jsonData)
{
    // Now, configure blocks that make up the flow domain.
    //
    // This is done in phases.  The blocks need valid references to LocalConfig objects
    // and the boundary conditions need valid references to FluidBlock objects.
    alias cfg = GlobalConfig;
    mixin(update_int("nfluidblock", "nFluidBlocks"));
    mixin(update_int("nfluidblockarrays", "nFluidBlockArrays"));
    foreach (i; 0 .. cfg.nFluidBlockArrays) {
        auto jsonDataForBlockArray = jsonData["fluid_block_array_" ~ to!string(i)];
        fluidBlockArrays ~= new FBArray(jsonDataForBlockArray);
    }
    if (cfg.verbosity_level > 1) {
        writeln("  nFluidBlocks: ", cfg.nFluidBlocks);
        writeln("  nFluidBlockArrays: ", cfg.nFluidBlockArrays);
        foreach (i, fba; fluidBlockArrays) {
            writefln("  fluid_block_array_%d: %s", i, fba);
        }
    }
    //
    // Set up dedicated copies of the configuration parameters for the threads.
    foreach (i; 0 .. cfg.nFluidBlocks) {
        dedicatedConfig ~= new LocalConfig(i);
    }
    foreach (i; 0 .. cfg.nFluidBlocks) {
        auto jsonDataForBlock = jsonData["block_" ~ to!string(i)];
        string gridType = getJSONstring(jsonDataForBlock, "grid_type", "");
        switch (gridType) {
        case "structured_grid":
            globalBlocks ~= new SFluidBlock(i, jsonDataForBlock);
            break;
        case "unstructured_grid":
            globalBlocks ~= new UFluidBlock(i, jsonDataForBlock);
            dedicatedConfig[i].stringent_cfl = true; // for signal_frequency calc in FVCell.
            break;
        default:
            throw new Error(format("Construction of fluidblock[%d], unknown grid type: %s", i, gridType));
        } // end switch gridType
    }
    // Defer the remaining configuration of the FluidBlocks until they have
    // been assigned to their MPI tasks out in the main part of init_simulation()
    // in simcore.d
    //
    // Read in any blocks in the solid domain.
    cfg.udfSolidSourceTerms = getJSONbool(jsonData, "udf_solid_source_terms", false);
    cfg.udfSolidSourceTermsFile = jsonData["udf_solid_source_terms_file"].str;
    cfg.nSolidBlocks = getJSONint(jsonData, "nsolidblock", 0);
    cfg.nBlocks = cfg.nFluidBlocks + cfg.nSolidBlocks;
    if (cfg.verbosity_level > 1) {
        writeln("  nSolidBlocks: ", cfg.nSolidBlocks);
        writeln("  udf_solid_source_terms: ", cfg.udfSolidSourceTerms);
        writeln("  udf_solid_source_terms_file: ", to!string(cfg.udfSolidSourceTermsFile));
    }
    // Set up dedicated copies of the configuration parameters for SolidBlocks, for the threads.
    foreach (int i; cfg.nFluidBlocks .. cfg.nBlocks) {
        dedicatedConfig ~= new LocalConfig(i);
    }
    // Note that we want the solidblock ids to continue on from the fluidblock ids.
    foreach (int i; cfg.nFluidBlocks .. cfg.nBlocks) {
        globalBlocks ~= new SSolidBlock(i, jsonData["solid_block_" ~ to!string(i)]);
        if (cfg.verbosity_level > 1) {
            writeln("  SolidBlock[", i, "]: ", globalBlocks[i]);
        }
    }
    foreach (int i; cfg.nFluidBlocks .. cfg.nBlocks) {
        // Note that we want the solidblock ids to continue on from the fluidblocks.
        auto sblk = cast(SSolidBlock) globalBlocks[i];
        assert(sblk !is null, "Oops, this should be a SolidBlock object.");
        sblk.initLuaGlobals();
        sblk.initBoundaryConditions(jsonData["solid_block_" ~ to!string(sblk.id)]);
        if (cfg.udfSolidSourceTerms) {
            initUDFSolidSourceTerms(sblk.myL, cfg.udfSolidSourceTermsFile);
        }
    }
    foreach (blk; globalBlocks) {
        auto myblk = cast(FluidBlock) blk;
        if (myblk) {
            myblk.init_boundary_conditions(jsonData["block_" ~ to!string(myblk.id)]);
        }
    }
    // Now that the blocks are partly configured, we can initialize
    // the lua_State that holds the user's functions
    // for simulation supervision and for defining grid motion.
    init_master_lua_State();
    if (cfg.user_pad_length > 0) {
        push_array_to_Lua(cfg.master_lua_State, cfg.userPad, "userPad");
    }
    if (cfg.udf_supervisor_file.length > 0) {
        doLuaFile(cfg.master_lua_State, cfg.udf_supervisor_file);
    }
    if (cfg.grid_motion == GridMotion.user_defined) {
        doLuaFile(cfg.master_lua_State, cfg.udf_grid_motion_file);
    }
} // end set_config_for_blocks()

void read_control_file()
{
    alias cfg = GlobalConfig;
    if (cfg.verbosity_level > 1) writeln("read_control_file()");
    JSONValue jsonData = readJSONfile("config/"~cfg.base_file_name~".control");
    mixin(update_double("dt_init", "dt_init"));
    mixin(update_double("dt_max", "dt_max"));
    // 2021-05-21 PJ: We no longer read cfl_value from .control file.
    // 2021-05-25 PJ: but we do allow a cfl scaling factor.
    mixin(update_double("cfl_scale_factor", "cfl_scale_factor"));
    mixin(update_bool("stringent_cfl", "stringent_cfl"));
    mixin(update_double("viscous_signal_factor", "viscous_signal_factor"));
    mixin(update_double("turbulent_signal_factor", "turbulent_signal_factor"));
    mixin(update_enum("residual_smoothing_type", "residual_smoothing_type", "residual_smoothing_type_from_name"));
    mixin(update_double("residual_smoothing_weight", "residual_smoothing_weight"));
    mixin(update_int("residual_smoothing_iterations", "residual_smoothing_iterations"));
    mixin(update_bool("fixed_time_step", "fixed_time_step"));
    mixin(update_int("print_count", "print_count"));
    mixin(update_int("cfl_count", "cfl_count"));
    mixin(update_double("max_time", "max_time"));
    mixin(update_int("max_step", "max_step"));
    mixin(update_double("dt_plot", "dt_plot"));
    mixin(update_double("dt_history", "dt_history"));
    mixin(update_double("dt_loads", "dt_loads"));
    mixin(update_int("write_loads_at_step", "write_loads_at_step"));
    mixin(update_int("write_flow_solution_at_step", "write_flow_solution_at_step"));
    mixin(update_int("snapshot_count", "snapshotCount"));
    mixin(update_int("number_total_snapshots", "nTotalSnapshots"));
    //
    mixin(update_int("halt_now", "halt_now"));
    //
    if (cfg.verbosity_level > 1) {
        writeln("  dt_init: ", cfg.dt_init);
        writeln("  dt_max: ", cfg.dt_max);
        writeln("  cfl_scale_factor: ", cfg.cfl_scale_factor);
        writeln("  stringent_cfl: ", cfg.stringent_cfl);
        writeln("  viscous_signal_factor: ", cfg.viscous_signal_factor);
        writeln("  turbulent_signal_factor: ", cfg.turbulent_signal_factor);
        writeln("  residual_smoothing_type: ", cfg.residual_smoothing_type);
        writeln("  residual_smoothing_weight: ", cfg.residual_smoothing_weight);
        writeln("  residual_smoothing_iterations: ", cfg.residual_smoothing_iterations);
        writeln("  fixed_time_step: ", cfg.fixed_time_step);
        writeln("  print_count: ", cfg.print_count);
        writeln("  cfl_count: ", cfg.cfl_count);
        writeln("  max_time: ", cfg.max_time);
        writeln("  max_step: ", cfg.max_step);
        writeln("  dt_plot: ", cfg.dt_plot);
        writeln("  dt_history: ", cfg.dt_history);
        writeln("  dt_loads: ", cfg.dt_loads);
        writeln("  write_loads_at_step: ", cfg.write_loads_at_step);
        writeln("  write_flow_solution_at_step: ", cfg.write_flow_solution_at_step);
        writeln("  snapshot_count: ", cfg.snapshotCount);
        writeln("  number_total_snapshots: ", cfg.nTotalSnapshots);
        writeln("  halt_now: ", cfg.halt_now);
    }
    //
    version(nk_accelerator) {
        auto sssOptions = jsonData["steady_state_solver_options"];
        auto ssso = &(cfg.sssOptions);
        ssso.usePreconditioner = getJSONbool(sssOptions, "use_preconditioner", ssso.usePreconditioner);
        ssso.frozenPreconditionerCount = getJSONint(sssOptions, "frozen_preconditioner_count", ssso.frozenPreconditionerCount);
        ssso.startPreconditioning = getJSONint(sssOptions, "start_preconditioning", ssso.startPreconditioning);
        ssso.iluFill = getJSONint(sssOptions, "ilu_fill", ssso.iluFill);
        auto mySaveValue1 = ssso.preconditionMatrixType;
        try {
            string name = sssOptions["precondition_matrix_type"].str;
            ssso.preconditionMatrixType = preconditionMatrixTypeFromName(name);
        } catch (Exception e) {
            ssso.preconditionMatrixType = mySaveValue1;
        }
        auto mySaveValue2 = ssso.preconditionMatrixFluxCalculator;
        try {
            string name = sssOptions["precondition_matrix_flux_calculator"].str;
            // if the user hasn't specified a specific flux calculator for the precondition matrix then we will just set it to the main flux calculator
            if ( name == "same as config.flux_calculator" ) {
                ssso.preconditionMatrixFluxCalculator = cfg.flux_calculator;
            } else {
                ssso.preconditionMatrixFluxCalculator = flux_calculator_from_name(name);
            }
        } catch (Exception e) {
            ssso.preconditionMatrixFluxCalculator = mySaveValue2;
        }
        ssso.preconditionerSigma = getJSONdouble(sssOptions, "preconditioner_sigma", ssso.preconditionerSigma);
        ssso.frozenLimiterOnLHS = getJSONbool(sssOptions, "frozen_limiter_on_lhs", ssso.frozenLimiterOnLHS);
        ssso.useAdaptivePreconditioner = getJSONbool(sssOptions, "use_adaptive_preconditioner", ssso.useAdaptivePreconditioner);
        ssso.usePhysicalityCheck = getJSONbool(sssOptions, "use_physicality_check", ssso.usePhysicalityCheck);
        ssso.physicalityCheckTheta = getJSONdouble(sssOptions, "physicality_check_theta", ssso.physicalityCheckTheta);
        ssso.useLineSearch = getJSONbool(sssOptions, "use_line_search", ssso.useLineSearch);
        ssso.inviscidCFL = getJSONbool(sssOptions, "inviscid_cfl", ssso.inviscidCFL);
        ssso.useScaling = getJSONbool(sssOptions, "use_scaling", ssso.useScaling);
        ssso.useComplexMatVecEval = getJSONbool(sssOptions, "use_complex_matvec_eval", ssso.useComplexMatVecEval);
        ssso.temporalIntegrationMode = getJSONint(sssOptions, "temporal_integration_mode", ssso.temporalIntegrationMode);
        ssso.nPreSteps = getJSONint(sssOptions, "number_pre_steps", ssso.nPreSteps);
        ssso.nTotalSteps = getJSONint(sssOptions, "number_total_steps", ssso.nTotalSteps);
        ssso.maxNumberAttempts = getJSONint(sssOptions, "max_number_attempts", ssso.maxNumberAttempts);
        ssso.stopOnRelGlobalResid = getJSONdouble(sssOptions, "stop_on_relative_global_residual", ssso.stopOnRelGlobalResid);
        ssso.stopOnAbsGlobalResid = getJSONdouble(sssOptions, "stop_on_absolute_global_residual", ssso.stopOnAbsGlobalResid);
        ssso.stopOnMassBalance    = getJSONdouble(sssOptions, "stop_on_mass_balance", ssso.stopOnMassBalance);
        ssso.maxSubIterations = getJSONint(sssOptions, "max_sub_iterations", ssso.maxSubIterations);
        ssso.maxOuterIterations = getJSONint(sssOptions, "max_outer_iterations", ssso.maxOuterIterations);
        ssso.maxRestarts = getJSONint(sssOptions, "max_restarts", ssso.maxRestarts);
        ssso.nInnerIterations = getJSONint(sssOptions, "number_inner_iterations", ssso.nInnerIterations);
        // Settings for start-up phase
        ssso.nStartUpSteps = getJSONint(sssOptions, "number_start_up_steps", ssso.nStartUpSteps);
        //
        ssso.cfl_schedule_length = getJSONint(sssOptions, "cfl_schedule_length", ssso.cfl_schedule_length);
        //
        double[] default_cfl_schedule_value_data;
        foreach (i; 0 .. ssso.cfl_schedule_length) { default_cfl_schedule_value_data ~= 0.0; }
        auto cfl_schedule_value_data = getJSONdoublearray(sssOptions, "cfl_schedule_value_list", default_cfl_schedule_value_data);
        ssso.cfl_schedule_value_list.length = ssso.cfl_schedule_length;
        foreach (i; 0 .. ssso.cfl_schedule_length) {
            ssso.cfl_schedule_value_list[i] = (i < cfl_schedule_value_data.length) ? cfl_schedule_value_data[i] : default_cfl_schedule_value_data[i];
        }
        //
        int[] default_cfl_schedule_iter_data;
        foreach (i; 0 .. ssso.cfl_schedule_length) { default_cfl_schedule_iter_data ~= 0; }
        auto cfl_schedule_iter_data = getJSONintarray(sssOptions, "cfl_schedule_iter_list", default_cfl_schedule_iter_data);
        ssso.cfl_schedule_iter_list.length = ssso.cfl_schedule_length;
        foreach (i; 0 .. ssso.cfl_schedule_length) {
            ssso.cfl_schedule_iter_list[i] = (i < cfl_schedule_iter_data.length) ? cfl_schedule_iter_data[i] : default_cfl_schedule_iter_data[i];
        }
        //
        ssso.include_turb_quantities_in_residual = getJSONbool(sssOptions, "include_turb_quantities_in_residual", ssso.include_turb_quantities_in_residual);
        ssso.residual_based_cfl_scheduling = getJSONbool(sssOptions, "residual_based_cfl_scheduling", ssso.residual_based_cfl_scheduling);
        ssso.cfl_max = getJSONdouble(sssOptions, "cfl_max", ssso.cfl_max);
        ssso.cfl_min = getJSONdouble(sssOptions, "cfl_min", ssso.cfl_min);
        ssso.LHSeval0 = getJSONint(sssOptions, "LHSeval0", ssso.LHSeval0);
        ssso.RHSeval0 = getJSONint(sssOptions, "RHSeval0", ssso.RHSeval0);
        ssso.cfl0 = getJSONdouble(sssOptions, "cfl0", ssso.cfl0);
        ssso.eta0 = getJSONdouble(sssOptions, "eta0", ssso.eta0);
        ssso.tau0 = getJSONdouble(sssOptions, "tau0", ssso.tau0);
        ssso.sigma0 = getJSONdouble(sssOptions, "sigma0", ssso.sigma0);
        ssso.p0 = getJSONdouble(sssOptions, "p0", ssso.p0);
        // Setting for inexact Newton phase
        ssso.LHSeval1 = getJSONint(sssOptions, "LHSeval1", ssso.LHSeval1);
        ssso.RHSeval1 = getJSONint(sssOptions, "RHSeval1", ssso.RHSeval1);
        ssso.cfl1 = getJSONdouble(sssOptions, "cfl1", ssso.cfl1);
        ssso.tau1 = getJSONdouble(sssOptions, "tau1", ssso.tau1);
        ssso.sigma1 = getJSONdouble(sssOptions, "sigma1", ssso.sigma1);
        ssso.p1 = getJSONdouble(sssOptions, "p1", ssso.p1);
        auto mySaveValue3 = ssso.etaStrategy;
        try {
            string name = sssOptions["eta_strategy"].str;
            ssso.etaStrategy = etaStrategyFromName(name);
        } catch (Exception e) {
            ssso.etaStrategy = mySaveValue3;
        }
        ssso.eta1 = getJSONdouble(sssOptions, "eta1", ssso.eta1);
        ssso.eta1_max = getJSONdouble(sssOptions, "eta1_max", ssso.eta1_max);
        ssso.eta1_min = getJSONdouble(sssOptions, "eta1_min", ssso.eta1_min);
        ssso.etaRatioPerStep = getJSONdouble(sssOptions, "eta_ratio_per_step", ssso.etaRatioPerStep);
        ssso.gamma = getJSONdouble(sssOptions, "gamma", ssso.gamma);
        ssso.alpha = getJSONdouble(sssOptions, "alpha", ssso.alpha);
        ssso.limiterFreezingResidReduction = getJSONdouble(sssOptions, "limiter_freezing_residual_reduction", ssso.limiterFreezingResidReduction);
        ssso.limiterFreezingCount = getJSONint(sssOptions, "limiter_freezing_count", ssso.limiterFreezingCount);
        // Settings for writing out snapshots and diagnostics
        ssso.snapshotsCount = getJSONint(sssOptions, "snapshots_count", ssso.snapshotsCount);
        ssso.nTotalSnapshots = getJSONint(sssOptions, "number_total_snapshots", ssso.nTotalSnapshots);
        ssso.writeDiagnosticsCount = getJSONint(sssOptions, "write_diagnostics_count", ssso.writeDiagnosticsCount);
        ssso.writeLoadsCount = getJSONint(sssOptions, "write_loads_count", ssso.writeLoadsCount);
    } // end version(nk_accelerator)
    //
    // Propagate new values to the local copies of config.
    foreach (localConfig; dedicatedConfig) {
        localConfig.update_control_parameters();
    }
} // end read_control_file()

//
// Functions to check the compatibility of the GlobalConfig parameters.
// The individual checks are called at particular locations during the
// reading of the config parameters from the JSON .config file.
//
// Sometimes, we can make a suitable change and continue (after emitting a warning)
// but other times it is not sensible to continue, but throw an exception.
//
void configCheckPoint1()
{
    alias cfg = GlobalConfig;
    if (cfg.suppress_reconstruction_at_shocks) {
        cfg.do_shock_detect = true;
    }
    if (cfg.flux_calculator == FluxCalculator.adaptive_hanel_ausmdv ||
        cfg.flux_calculator == FluxCalculator.adaptive_hanel_ausm_plus_up ||
        cfg.flux_calculator == FluxCalculator.adaptive_ldfss0_ldfss2 ||
        cfg.flux_calculator == FluxCalculator.adaptive_hlle_hllc ||
        cfg.flux_calculator == FluxCalculator.adaptive_hlle_roe ||
        cfg.flux_calculator == FluxCalculator.adaptive_efm_ausmdv ||
        cfg.flux_calculator == FluxCalculator.adaptive_ausmdv_asf) {
        cfg.do_shock_detect = true;
    }
    if (cfg.do_shock_detect && cfg.is_master_task) {
        writeln("NOTE: shock detector is on.");
    }
    if (!cfg.high_order_flux_calculator) {
        if (cfg.flux_calculator == FluxCalculator.adaptive_ausmdv_asf) {
            cfg.high_order_flux_calculator = true;
        }
    }
    return;
} // end configCheckPoint1()

void configCheckPoint2()
{
    alias cfg = GlobalConfig;
    // More checking of constraints on the config parameters.
    version(nk_accelerator) {
        if (cfg.grid_motion != GridMotion.none) {
            throw new Error("Grid motion is not compatible e4-nk-dist.");
        }
    }
    if (cfg.grid_motion != GridMotion.none) {
        if (cfg.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_1_stage ||
            cfg.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_2_stage ||
            cfg.gasdynamic_update_scheme == GasdynamicUpdate.backward_euler ||
            cfg.gasdynamic_update_scheme == GasdynamicUpdate.implicit_rk1) {
            // pass, we have a consistent selection.
        } else {
            string msg = "We have some grid_motion but not a valid GasdynamicUpdate scheme for grid motion.";
            throw new FlowSolverException(msg);
        }
    } else {
        if (cfg.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_1_stage ||
            cfg.gasdynamic_update_scheme == GasdynamicUpdate.moving_grid_2_stage) {
            string msg = "We have no grid_motion but have asked for a GasdynamicUpdate scheme for grid motion.";
            throw new FlowSolverException(msg);
        }
    }
    if (cfg.high_order_flux_calculator) {
        if (cfg.n_ghost_cell_layers < 3) {
            if (cfg.is_master_task) {
                writeln("NOTE: Increasing n_ghost_cell_layers to 3.");
            }
            cfg.n_ghost_cell_layers = 3;
        }
    }
    return;
} // end configCheckPoint2()

void configCheckPoint3()
{
    alias cfg = GlobalConfig;
    // Check the compatibility of the gas model if mass diffusion is selected.
    if (cfg.mass_diffusion_model != MassDiffusionModel.none) {
        if (cfg.gmodel_master.n_species == 1) {
            string msg = format("The selected mass diffusion model '%s'",
                                massDiffusionModelName(cfg.mass_diffusion_model));
            msg ~= " makes no sense when number of species = 1.\n";
            throw new FlowSolverException(msg);
        }
    }
    // Check the compatibility of turbulence model selection and flux calculator.
    if (cfg.turb_model && cfg.turb_model.isTurbulent) {
        if (cfg.flux_calculator == FluxCalculator.hlle) {
            string msg = format("The selected flux calculator '%s'",
                                flux_calculator_name(cfg.flux_calculator));
            msg ~= " is incompatible with an active turbulence model.";
            throw new FlowSolverException(msg);
        }
    }
    //
    // Check the compatibility of update scheme and viscous flag.
    if (cfg.viscous == false &&
        (cfg.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 ||
         cfg.gasdynamic_update_scheme == GasdynamicUpdate.rkl2)) {
        string msg = format("The selected gas dynamic update scheme '%s'",
                            gasdynamic_update_scheme_name(cfg.gasdynamic_update_scheme));
        msg ~= " is incompatible with an inviscid simulation.";
        throw new FlowSolverException(msg);
    }
    // The super_time_stepping is associated with two particular update schemes:
    // rkl1, rkl2.
    if (cfg.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 ||
        cfg.gasdynamic_update_scheme == GasdynamicUpdate.rkl2) {
        cfg.with_super_time_stepping = true;
    }
    //
    // Check for compatbility between viscous flag and turbulence model
    if (cfg.turb_model && cfg.turb_model.isTurbulent &&
        (cfg.viscous==false)) {
        string msg = format("The selected turbulence model '%s'",
                            cfg.turbulence_model_name);
        msg ~= " is incompatible with an inviscid simulation.";
        throw new FlowSolverException(msg);
    }
    if (cfg.turb_model is null) {
        // PJ 2020-10-08 Have relaxed the Exception that was thrown here to just a warning.
        // This allows an essentially empty input script to be processed by e4shared --prep.
        writeln("Warning: GlobalConfig does not have a turbulence model.");
    }
    //
    if (cfg.spatial_deriv_calc == SpatialDerivCalc.divergence) {
        // The divergence method is the old default for the type of gradient calculation
        // while 'cells' is the new default for gradient location and together they are
        // incomplete.  For least bother, let's alter the calculation type.
        if (cfg.dimensions == 3 ||
            cfg.spatial_deriv_locn == SpatialDerivLocn.cells) {
            if (cfg.is_master_task) {
                writeln("NOTE: Changing spatial_deriv_calc to least_squares.");
            }
            cfg.spatial_deriv_calc = SpatialDerivCalc.least_squares;
        }
    }
    return;
} // end configCheckPoint3()

void configCheckPoint4()
{
    alias cfg = GlobalConfig;
    // The shape sensitivity calculator shouldn't apply diffuse_bcs_on_init_flag.
    version(shape_sensitivity) {
        cfg.n_grid_time_levels = 3;
    }
}

void checkGlobalConfig()
// Bundle all of the checks together so that they may be conveniently applied
// at the end of processing the user's Lua input script.
{
    configCheckPoint1();
    configCheckPoint2();
    configCheckPoint3();
    configCheckPoint4();
    return;
}


// -----------------------
// PART 5. Lua interaction
// -----------------------

void init_master_lua_State()
{
    alias cfg = GlobalConfig;
    cfg.master_lua_State = init_lua_State();
    // Give me a conveniently-named pointer for use in this function.
    auto L = cfg.master_lua_State;
    registerVector3(L);
    registerBBLA(L);
    // Load some Lua modules using 'require'.
    // There is no convenient C API expression to do the equivalent of "require"
    luaL_dostring(L, "require 'lua_helper'");
    // Set some globally available constants for the Lua state.
    lua_pushboolean(L, SimState.is_restart); lua_setglobal(L, "is_restart");
    lua_pushnumber(L, cfg.nFluidBlocks); lua_setglobal(L, "nFluidBlocks");
    lua_pushnumber(L, cfg.n_ghost_cell_layers); lua_setglobal(L, "n_ghost_cell_layers"); // interpreters for blocks use this name
    lua_pushnumber(L, cfg.n_ghost_cell_layers); lua_setglobal(L, "nGhostCellLayers"); // keep both names
    // Give the user a table that holds information about
    // all of the blocks in the full simulation.
    // Note that not all of these blocks may be fully present
    // in an MPI-parallel simulation.
    lua_newtable(L);
    foreach (i; 0..cfg.nFluidBlocks) {
        auto fluidblk = cast(FluidBlock) globalBlocks[i];
        assert(fluidblk !is null, "Oops, this should be a FluidBlock object.");
        lua_newtable(L);
        lua_pushnumber(L, fluidblk.cells.length); lua_setfield(L, -2, "nCells");
        lua_pushnumber(L, fluidblk.vertices.length); lua_setfield(L, -2, "nVertices");
        if (fluidblk.grid_type == Grid_t.structured_grid) {
            auto sblk = cast(SFluidBlock) fluidblk;
            assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
            lua_pushnumber(L, sblk.nic); lua_setfield(L, -2, "niCells");
            lua_pushnumber(L, sblk.njc); lua_setfield(L, -2, "njCells");
            lua_pushnumber(L, sblk.nkc); lua_setfield(L, -2, "nkCells");
            lua_pushnumber(L, 0); lua_setfield(L, -2, "vtxImin");
            lua_pushnumber(L, sblk.nic); lua_setfield(L, -2, "vtxImax");
            lua_pushnumber(L, 0); lua_setfield(L, -2, "vtxJmin");
            lua_pushnumber(L, sblk.njc); lua_setfield(L, -2, "vtxJmax");
            lua_pushnumber(L, 0); lua_setfield(L, -2, "vtxKmin");
            lua_pushnumber(L, sblk.nkc); lua_setfield(L, -2, "vtxKmax");
        }
        lua_rawseti(L, -2, to!int(i));
    }
    lua_setglobal(L, "blockData");
    //
    setSampleHelperFunctions(L);
    setGridMotionHelperFunctions(L);
} // end init_master_lua_State()


// Functions related to the managed gas model.

extern(C) int setGasModel(lua_State* L)
{
    if (lua_isstring(L, 1)) {
        alias cfg = GlobalConfig;
        string fname = to!string(luaL_checkstring(L, 1));
        cfg.gas_model_file = fname;
        // 2021-05-22 Also, to set config.gas_model_name in Lua world.
        lua_getglobal(L, "config".toStringz);
        if (lua_istable(L, -1)) {
            lua_pushstring(L, fname.toStringz);
            lua_setfield(L, -2, "gas_model_file".toStringz);
        }
        lua_pop(L, 1); // dispose of config
        try {
            cfg.gmodel_master = init_gas_model(fname);
        } catch (GasModelException e) {
            string msg = "\nThere is a problem in call to setGasModel. Reported errors are:\n";
            msg ~= e.msg;
            msg ~= "\n---------------------------\n";
            msg ~= "The preparation stage cannot proceed. Exiting without completing.\n";
            writeln(msg);
            exit(1);
        }
        lua_settop(L, 0); // clear the stack
        lua_pushinteger(L, cfg.gmodel_master.n_species);
        lua_pushinteger(L, cfg.gmodel_master.n_modes);
        GasModelStore ~= pushObj!(GasModel, GasModelMT)(L, cfg.gmodel_master);
        return 3;
    } else {
        string msg = "setGasModel expects a string as the name of the gas model file.";
        luaL_error(L, msg.toStringz);
        return 0;
    }
}

extern(C) int getGasModel(lua_State* L)
{
    alias cfg = GlobalConfig;
    if (cfg.gmodel_master is null) {
        string msg = "The master gas model appears to be uninitialized. "~
            "You should initialize it with setGasModel().";
        luaL_error(L, msg.toStringz);
    }
    lua_settop(L, 0); // clear the stack
    pushObj!(GasModel, GasModelMT)(L, cfg.gmodel_master);
    return 1;
}

//-----------------------------------------------------------------------
// Call the following function from the main program to get the
// functions appearing in the Lua interpreter.

void registerGlobalConfig(lua_State* L)
{
    // Register global functions related to the managed gas model.
    lua_pushcfunction(L, &setGasModel);
    lua_setglobal(L, "setGasModel");
    lua_pushcfunction(L, &getGasModel);
    lua_setglobal(L, "getGasModel");
} // end registerGlobalConfig()
