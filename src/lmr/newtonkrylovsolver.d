/** 
 * Core set of functions used in the Newton-Krylov updates for steady-state convergence.
 *
 * Authors: RJG and KAD
 * Date: 2022-02-28
 * History:
 *   2022-02-28 Fairly aggressive refactor of code that was in: steadystate_core.d
 */

module newtonkrylovsolver;

import core.stdc.stdlib : exit;
import core.memory : GC;
import std.algorithm : min;
import std.algorithm.searching : countUntil;
import std.datetime : Clock;
import std.parallelism : parallel, defaultPoolThreads;
import std.stdio : File, writeln, writefln, stdout;
import std.file : rename;
import std.array : appender;
import std.format : formattedWrite;
import std.json : JSONValue;
import std.math;
import std.string;
import std.conv : to;
import std.typecons : Tuple, tuple;

import nm.complex;
import nm.number : number;
import nm.smla;
import nm.bbla;
import util.time_utils : timeStringToSeconds;
import util.lua;
import util.lua_service;
import lua_helper;
import json_helper;
import geom;

import lmrconfig;
import conservedquantities : ConservedQuantities;
import fileutil : ensure_directory_is_present;

import globalconfig;
import globaldata;
import simcore : synchronize_corner_coords_for_all_blocks, compute_wall_distances;
import simcore_exchange;
import simcore_gasdynamic_step : detect_shocks;
import bc;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;
import user_defined_source_terms : getUDFSourceTermsForCell;
import fluidblockio_new : read_zip_solution;

version(mpi_parallel) {
    import mpi;
}

/*---------------------------------------------------------------------
 * Module globals
 * Typically for diagnostics and/or debugging.
 *---------------------------------------------------------------------
 */

File diagnostics;
string diagnosticsDir = "diagnostics";
string diagnosticsFilename = "lmr-nk-diagnostics.dat";

/*---------------------------------------------------------------------
 * Exception class to signal N-K specific exceptions.
 *---------------------------------------------------------------------
 */
class NewtonKrylovException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

/*---------------------------------------------------------------------
 * Enums for preconditioners
 *---------------------------------------------------------------------
 */
enum PreconditionerType { lusgs, diagonal, jacobi, sgs, ilu }

string preconditionerTypeName(PreconditionerType i)
{
    final switch (i) {
    case PreconditionerType.lusgs: return "lusgs";
    case PreconditionerType.diagonal: return "diagonal";
    case PreconditionerType.jacobi: return "jacobi";
    case PreconditionerType.sgs: return "sgs";
    case PreconditionerType.ilu: return "ilu";
    }
} 

PreconditionerType preconditionerTypeFromName(string name)
{
    switch (name) {
    case "lusgs", "lu_sgs": return PreconditionerType.lusgs;
    case "diagonal": return PreconditionerType.diagonal;
    case "jacobi": return PreconditionerType.jacobi;
    case "sgs": return PreconditionerType.sgs;
    case "ilu": return PreconditionerType.ilu;
    default:
        string errMsg = "The selected 'preconditioner' is unavailable.\n";
        errMsg ~= format("You selected: '%s'\n", name);
        errMsg ~= "The available strategies are: \n";
        errMsg ~= "   'lusgs'\n";
        errMsg ~= "   'diagonal'\n";
        errMsg ~= "   'jacobi'\n";
        errMsg ~= "   'sgs'\n";
        errMsg ~= "   'ilu'\n";
        errMsg ~= "Check your selection or its spelling in the input file.\n";
        throw new Error(errMsg);
    }
} 


/*---------------------------------------------------------------------
 * Structs for configuration: global and phase-specific
 *---------------------------------------------------------------------
 */

struct NKGlobalConfig {
    // global control based on step
    int numberOfStepsForSettingReferenceResiduals = 10;
    int freezeLimiterAtStep = -1;
    // stopping criterion
    int maxNewtonSteps = 1000;
    double stopOnRelativeResidual = 1.0e-99;
    double stopOnAbsoluteResidual = 1.0e-99;
    double stopOnMassBalance = -1.0;
    int maxConsecutiveBadSteps = 2;
    // CFL control
    double cflMax = 1.0e8;
    double cflMin = 0.001;
    Tuple!(int, "step", double, "cfl")[] cflSchedule;
    double cflReductionFactor = 0.5;
    // phase control
    int numberOfPhases = 1;
    int[] phaseChangesAtSteps;
    // Newton stepping and continuation
    bool useLocalTimestep = true;
    bool inviscidCFLOnly = true;
    bool useLineSearch = true;
    bool usePhysicalityCheck = true;
    double allowableRelativeMassChange = 0.2;
    double minRelaxationFactor = 0.1;
    double relaxationFactorReductionFactor = 0.7;
    // Linear solver and preconditioning
    int maxLinearSolverIterations = 10;
    int maxLinearSolverRestarts = 0;
    bool useScaling = true;
    double frechetDerivativePerturbation = 1.0e-30;
    bool usePreconditioner = true;
    double preconditionerPerturbation = 1.0e-30;
    PreconditionerType preconditioner = PreconditionerType.ilu;
    // ILU setting
    int iluFill = 0;
    // sub iterations for preconditioners
    int preconditionerSubIterations = 4;
    // output and diagnostics
    int totalSnapshots = 5;
    int stepsBetweenSnapshots = 10;
    int stepsBetweenDiagnostics = 10;
    int stepsBetweenLoadsUpdate = 20;

    void readValuesFromJSON(JSONValue jsonData)
    {
        numberOfStepsForSettingReferenceResiduals = getJSONint(jsonData, "number_of_steps_for_setting_reference_residuals", numberOfStepsForSettingReferenceResiduals);
        freezeLimiterAtStep = getJSONint(jsonData, "freeze_limiter_at_step", freezeLimiterAtStep);
        maxNewtonSteps = getJSONint(jsonData, "max_newton_steps", maxNewtonSteps);
        stopOnRelativeResidual = getJSONdouble(jsonData, "stop_on_relative_residual", stopOnRelativeResidual);
        stopOnAbsoluteResidual = getJSONdouble(jsonData, "stop_on_absolute_residual", stopOnAbsoluteResidual);
        stopOnMassBalance = getJSONdouble(jsonData, "stop_on_mass_balance", stopOnMassBalance);
        maxConsecutiveBadSteps = getJSONint(jsonData, "max_consecutive_bad_steps", maxConsecutiveBadSteps);
        cflMax = getJSONdouble(jsonData, "max_cfl", cflMax);
        cflMin = getJSONdouble(jsonData, "min_cfl", cflMin);
        auto jsonArray = jsonData["cfl_schedule"].array;
        foreach (entry; jsonArray) {
            auto values = entry.array;
            cflSchedule ~= tuple!("step", "cfl")(values[0].get!int, values[1].get!double);
        }
        cflReductionFactor = getJSONdouble(jsonData, "cfl_reduction_factor", cflReductionFactor);
        numberOfPhases = getJSONint(jsonData, "number_of_phases", numberOfPhases);
        phaseChangesAtSteps = getJSONintarray(jsonData, "phase_changes_at_steps", phaseChangesAtSteps);
        useLocalTimestep = getJSONbool(jsonData, "use_local_timestep", useLocalTimestep);
        inviscidCFLOnly = getJSONbool(jsonData, "inviscid_cfl_only", inviscidCFLOnly);
        useLineSearch = getJSONbool(jsonData, "use_line_search", useLineSearch);
        usePhysicalityCheck = getJSONbool(jsonData, "use_physicality_check", usePhysicalityCheck);
        allowableRelativeMassChange = getJSONdouble(jsonData, "allowable_relative_mass_change", allowableRelativeMassChange);
        minRelaxationFactor = getJSONdouble(jsonData, "min_relaxation_factor", minRelaxationFactor);
        relaxationFactorReductionFactor = getJSONdouble(jsonData, "relaxation_factor_reduction_factor", relaxationFactorReductionFactor);
        maxLinearSolverIterations = getJSONint(jsonData, "max_linear_solver_iterations", maxLinearSolverIterations);
        maxLinearSolverRestarts = getJSONint(jsonData, "max_linear_solver_restarts", maxLinearSolverRestarts);
        useScaling = getJSONbool(jsonData, "use_scaling", useScaling);
        frechetDerivativePerturbation = getJSONdouble(jsonData, "frechet_derivative_perturbation", frechetDerivativePerturbation);
        usePreconditioner = getJSONbool(jsonData, "use_preconditioner", usePreconditioner);
        preconditionerPerturbation = getJSONdouble(jsonData, "preconditioner_perturbation", preconditionerPerturbation);
        auto pString = getJSONstring(jsonData, "preconditioner", "NO_SELECTION_SUPPLIED");
        preconditioner = preconditionerTypeFromName(pString);
        iluFill = getJSONint(jsonData, "ilu_fill", iluFill);
        preconditionerSubIterations = getJSONint(jsonData, "preconditioner_sub_iterations", preconditionerSubIterations);
        totalSnapshots = getJSONint(jsonData, "total_snapshots", totalSnapshots);
        stepsBetweenSnapshots = getJSONint(jsonData, "steps_between_snapshots", stepsBetweenSnapshots);
        stepsBetweenDiagnostics = getJSONint(jsonData, "steps_between_diagnostics", stepsBetweenDiagnostics);
        stepsBetweenLoadsUpdate = getJSONint(jsonData, "steps_between_loads_update", stepsBetweenLoadsUpdate);
    }
}
NKGlobalConfig nkCfg;

struct NKPhaseConfig {
    int residualInterpolationOrder = 2;
    int jacobianInterpolationOrder = 2;
    bool frozenPreconditioner = true;
    int stepsBetweenPreconditionerUpdate = 10;
    bool useAdaptivePreconditioner = false;
    bool ignoreStoppingCriteria = true;
    bool frozenLimiterForJacobian = true;
    double linearSolveTolerance = 0.01;
    // Auto CFL control
    bool useAutoCFL = false;
    double thresholdRelativeResidualForCFLGrowth = 0.99;
    double startCFL = 1.0;
    double maxCFL = 1000.0;
    double autoCFLExponent = 0.75;

    void readValuesFromJSON(JSONValue jsonData)
    {
        residualInterpolationOrder = getJSONint(jsonData, "residual_interpolation_order", residualInterpolationOrder);
        jacobianInterpolationOrder = getJSONint(jsonData, "jacobian_interpolation_order", jacobianInterpolationOrder);
        frozenPreconditioner = getJSONbool(jsonData, "frozen_preconditioner", frozenPreconditioner);
        stepsBetweenPreconditionerUpdate = getJSONint(jsonData, "steps_between_preconditioner_update", stepsBetweenPreconditionerUpdate);
        useAdaptivePreconditioner = getJSONbool(jsonData, "use_adaptive_preconditioner", useAdaptivePreconditioner);
        ignoreStoppingCriteria = getJSONbool(jsonData, "ignore_stopping_criteria", ignoreStoppingCriteria);
        frozenLimiterForJacobian = getJSONbool(jsonData, "frozen_limiter_for_jacobian", frozenLimiterForJacobian);
        linearSolveTolerance = getJSONdouble(jsonData, "linear_solver_tolerance", linearSolveTolerance);
        useAutoCFL = getJSONbool(jsonData, "use_auto_cfl", useAutoCFL);
        thresholdRelativeResidualForCFLGrowth = getJSONdouble(jsonData, "threshold_relative_residual_for_cfl_growth", thresholdRelativeResidualForCFLGrowth);
        startCFL = getJSONdouble(jsonData, "start_cfl", startCFL);
        maxCFL = getJSONdouble(jsonData, "max_cfl", maxCFL);
        autoCFLExponent = getJSONdouble(jsonData, "auto_cfl_exponent", autoCFLExponent);
    }
}

NKPhaseConfig[] nkPhases;
NKPhaseConfig activePhase;

/*---------------------------------------------------------------------
 * Classes to handle CFL selection
 *---------------------------------------------------------------------
 */

/**
 * Interface to define the behaviour of a generic CFL selector.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-08
 */
interface CFLSelector {
    @nogc double nextCFL(double cfl, int step, double currResidual, double prevResidual, double relativeGlobalResidual);
}

/**
 * A CFL selctor that simply returns a linear interpolation between
 * start and end points (in terms of steps).
 *
 * The CFL is constant off the ends the step range, as shown below.
 *
 *                   endStep    
 *                       |            cflEnd
 *                       +------------------
 *                      /
 *                     /
 *                    /
 *                   /
 *                  /
 *                 /
 *                /
 *               /
 *   cflStart   /
 * ------------+
 *             |
 *          startStep
 *
 * Authors: RJG and KAD
 * Date: 2022-03-08
 */

class LinearRampCFL : CFLSelector {
    this(int startStep, int endStep, double startCFL, double endCFL)
    {
        mStartStep = startStep;
        mEndStep = endStep;
        mStartCFL = startCFL;
        mEndCFL = endCFL;
    }

    @nogc
    override double nextCFL(double cfl, int step, double currResidual, double prevResidual, double relativeGlobalResidual)
    {
        if (step <= mStartStep) return mStartCFL;
        if (step >= mEndStep) return mEndCFL;
        double frac = (step - mStartStep)/(cast(double)(mEndStep - mStartStep));
        return (1.0-frac)*mStartCFL + frac*mEndCFL;
    }

private:
    int mStartStep;
    int mEndStep;
    double mStartCFL;
    double mEndCFL;
}

/**
 * CFL selector based on residual drop.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-08
 */

class ResidualBasedAutoCFL : CFLSelector {
    this(double p, double cfl_max, double thresholdResidualDrop)
    {
        mP = p;
        mMaxCFL = cfl_max;
        mThresholdResidualDrop = thresholdResidualDrop;
    }

    @nogc
    override double nextCFL(double cfl, int step, double currResidual, double prevResidual, double relativeGlobalResidual)
    {
        if (relativeGlobalResidual > mThresholdResidualDrop) {
            // No adjustment required yet, return what we are given.
            return cfl;
        }
        auto residRatio = prevResidual/currResidual;
        double cflTrial = cfl*pow(residRatio, mP);
        // Apply some safeguards on the value.
        cflTrial = fmin(cflTrial, mLimitOnCFLIncreaseRatio*cfl);
        cflTrial = fmax(cflTrial, mLimitOnCFLDecreaseRatio*cfl);
        cflTrial = fmin(cflTrial, mMaxCFL);
        return cflTrial;
    }

private:
    immutable double mLimitOnCFLIncreaseRatio = 2.0;
    immutable double mLimitOnCFLDecreaseRatio = 0.1;
    double mP;
    double mMaxCFL;
    double mThresholdResidualDrop;
}



/*---------------------------------------------------------------------
 * Module-local globals
 *---------------------------------------------------------------------
 */

static int fnCount = 0;
immutable double dummySimTime = -1.0;
immutable double minScaleFactor = 1.0;
immutable string refResidFname = "config/reference-residuals.saved";
    
ConservedQuantities referenceResiduals, currentResiduals, scale;
double referenceGlobalResidual, globalResidual, prevGlobalResidual;
bool residualsUpToDate = false;
bool referenceResidualsAreSet = false;


// Module-local, global memory arrays and matrices
// TODO: Think about these, maybe they shouldn't be globals
number[] g0;
number[] g1;
number[] h;
number[] hR;
Matrix!number H0;
Matrix!number H1;
Matrix!number Gamma;
Matrix!number Q0;
Matrix!number Q1;

/*---------------------------------------------------------------------
 * Locally used data structures
 *---------------------------------------------------------------------
 */

struct RestartInfo {
    double pseudoSimTime;
    double dt;
    double cfl;
    int step;
    double globalResidual;
    ConservedQuantities residuals;

    this(size_t n)
    {
        residuals = new ConservedQuantities(n);
    }
}
RestartInfo[] snapshots;

/*
struct LinearSystemInput {
    int step;
    double dt;
    bool computePreconditioner;
};
LinearSystemInput lsi;
*/

struct GMRESInfo {
    int nRestarts;
    double initResidual;
    double finalResidual;
    int iterationCount;
};
GMRESInfo gmresInfo;

/*--------------------------------------------------------------------
 * Initialisation for NK solve.
 *--------------------------------------------------------------------
 */

void initNewtonKrylovSimulation(int snapshotStart, int maxCPUs, int threadsPerMPITask, string maxWallClock)
{
    alias cfg = GlobalConfig;
    if (cfg.verbosity_level > 0 && cfg.is_master_task) {
        writeln("Begin initNewtonKrylovSimulation()...");
    }
    // Start timer
    SimState.maxWallClockSeconds = timeStringToSeconds(maxWallClock);
    SimState.wall_clock_start = Clock.currTime();

    version(enable_fp_exceptions) {
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
    }

    // Initialise baseline configuration
    initConfiguration();
    if (cfg.nFluidBlocks == 0 && cfg.is_master_task) {
        throw new NewtonKrylovException("No FluidBlocks; no point in continuing with simulation initialisation.");
    }
    cfg.n_flow_time_levels = 2;
    
    initLocalFluidBlocks();
    initBlockIDs();

    initThreadPool(maxCPUs);

    initFluidBlocksBasic();
    initFluidBlocksGridsAndGeom();
    initFluidBlocksGlobalCellIDs();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    initFullFaceDataExchange();
    initMappedCellDataExchange();
    initGhostCellGeometry();
    initLeastSquaresStencils();

    if ((cfg.interpolation_order > 1) &&
	((cfg.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp) ||
	 (cfg.unstructured_limiter == UnstructuredLimiter.venkat_mlp))) {
        initMLPlimiter();
    }

    // [TODO] Think about whether re-ordering localBlocksBySize is needed.

    initMasterLuaState();
    initCornerCoordinates();
    if (cfg.turb_model.needs_dwall) initWallDistances();

    // [TODO] Add in electric field solver initialisation.

    // Do some memory clean-up and reporting.
    GC.collect();
    GC.minimize();
    debug {
        if (cfg.verbosity_level > 0) {
            auto myStats = GC.stats();
            auto heapUsed = to!double(myStats.usedSize)/(2^^20);
            auto heapFree = to!double(myStats.freeSize)/(2^^20);
            writefln("Heap memory used for task %d: %.2f  free: %.2f  total: %.1f MB",
                     cfg.mpi_rank_for_local_task, heapUsed, heapFree, heapUsed+heapFree);
            stdout.flush();
        }
    }

    if (cfg.verbosity_level > 0 && cfg.is_master_task) {
        // For reporting wall-clock time, convert to seconds with precision of milliseconds.
        double wall_clock_elapsed = to!double((Clock.currTime() - SimState.wall_clock_start).total!"msecs"())/1000.0;
        writefln("Done initNewtonKrylovSimulation() at wall-clock(WC)= %.1f sec", wall_clock_elapsed);
        stdout.flush();
    }
}

void initConfiguration()
{
    // Read in config file and set parameters
    auto cfgData = readJSONfile(lmrCfg.cfgFile);
    set_config_for_core(cfgData);
    set_config_for_blocks(cfgData);
}

void initLocalFluidBlocks()
{
    // [TODO] mpi version
    foreach (blk; globalBlocks) {
        auto myfblk = cast(FluidBlock) blk;
        if (myfblk) { localFluidBlocks ~= myfblk; }
        /+ [TODO] add in solid blocks 
        auto mysblk = cast(SSolidBlock) blk;
        if (mysblk) { localSolidBlocks ~= mysblk; }
        +/
    }
}

void initBlockIDs()
{
    alias cfg = GlobalConfig;
    foreach (blk; localFluidBlocks) cfg.localFluidBlockIds ~= blk.id;
    // [TODO] add solid blocks here
}

void initThreadPool(int maxCPUs)
{
    auto nBlocksInThreadParallel = localFluidBlocks.length; // [TODO] add solid blocks
    int extraThreadsInPool;
    // [TODO] handle MPI case
    extraThreadsInPool = min(maxCPUs-1, nBlocksInThreadParallel-1);
    defaultPoolThreads(extraThreadsInPool); // total = main thread + extra-threads-in-Pool
    if (GlobalConfig.verbosity_level > 0) {
        writeln("Single process running with ", extraThreadsInPool+1, " threads.");
    }
}

void initFluidBlocksBasic()
{
    foreach (myblk; localFluidBlocks) {
        myblk.myConfig.init_gas_model_bits();
        myblk.init_workspace();
        myblk.init_lua_globals();
        foreach (bci; myblk.bc) { bci.post_bc_construction(); }
        // NOTE: Removed userPad in NK solver.
        if (GlobalConfig.udf_source_terms) {
            luaL_dofile(myblk.myL, GlobalConfig.udf_source_terms_file.toStringz);
        }
        // After fully constructing the blocks and its boundary conditions,
        // we can optionally print their representation for checking.
        if (GlobalConfig.verbosity_level > 1) {
            writeln("  Block[", myblk.id, "]: ", myblk);
        }
    }
}

void initFluidBlocksGridsAndGeom()
{
    bool anyBlockFail = false;
    foreach (blk; parallel(localFluidBlocks,1)) {
        try {
            string gName = gridFilenameWithoutExt(blk.id);
            if (GlobalConfig.grid_format == "gziptext") {
                gName ~= "." ~ lmrCfg.gzipExt;
            }
            else if (GlobalConfig.grid_format == "rawbinary") {
                gName ~= "." ~ lmrCfg.rawBinExt;
            }
            else {
                throw new Error(format("Oops, invalid grid_format: %s", GlobalConfig.grid_format));
            }
            debug { writeln("Calling init_grid_and_flow_arrays for grid: ", gName); }
            blk.init_grid_and_flow_arrays(gName);
            blk.compute_primary_cell_geometric_data(0);
            blk.add_IO();
        }
        catch (Exception e) {
            writefln("Block[%d] failed to initialise geometry, msg=%s", blk.id, e.msg);
            anyBlockFail = true;
        }
    }
    if (anyBlockFail) {
        throw new NewtonKrylovException("Failed at initialisation stage during grid reading and geometry calculations.");
    }
}

void initFluidBlocksGlobalCellIDs()
{
    // Note that the global id is across all processes, not just the local collection of blocks.
    foreach (i, blk; globalBlocks) {
        auto fluidblk = cast(FluidBlock) blk;
        if (fluidblk) {
            if (i == 0) {
                fluidblk.globalCellIdStart = 0;
            } else {
                auto prev_fluidblk = cast(FluidBlock) globalBlocks[i-1];
                fluidblk.globalCellIdStart = prev_fluidblk.globalCellIdStart + prev_fluidblk.ncells_expected;
            }
        }
    }
}

void initFluidBlocksZones()
{
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.identify_reaction_zones(0);
        blk.identify_turbulent_zones(0);
        blk.identify_suppress_reconstruction_zones();
        blk.identify_suppress_viscous_stresses_zones();
    }
}

void initFluidBlocksFlowField(int snapshotStart)
{
    bool anyBlockFail = false;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.read_zip_solution(steadyFlowFilename(snapshotStart, blk.id));
        foreach (iface; blk.faces) iface.gvel.clear();
        foreach (cell; blk.cells) {
            cell.encode_conserved(0, 0, blk.omegaz);
            // Even though the following call appears redundant at this point,
            // fills in some gas properties such as Prandtl number that is
            // needed for both the cfl_check and the BaldwinLomax turbulence model.
            if (0 != cell.decode_conserved(0, 0, blk.omegaz)) {
                writefln("Block[%d] Bad cell decode_conserved while initialising flow.", blk.id);
                anyBlockFail = true;
            }
        }
        blk.set_cell_dt_chem(-1.0);
    }
    // [TODO] mpi reduce anyBlockFail
    if (anyBlockFail) {
        throw new NewtonKrylovException("Failed at initialisation stage during flow field initialisation.");
    }
}

void initFullFaceDataExchange()
{
    bool anyBlockFail = false;
    foreach (blk; localFluidBlocks) {
        foreach (j, bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto my_gce = cast(GhostCellFullFaceCopy)gce;
                if (my_gce) {
                    // The local block thinks that it has an exchange boundary with another block,
                    // so we need to check the ghost-cell effects of the other block's face to see
                    // that it points back to the local block face.
                    auto other_blk = my_gce.neighbourBlock;
                    bool ok = false;
                    auto other_blk_bc = other_blk.bc[my_gce.neighbourFace];
                    foreach (gce2; other_blk_bc.preReconAction) {
                        auto other_gce = cast(GhostCellFullFaceCopy)gce2;
                        if (other_gce &&
                            (other_gce.neighbourBlock.id == blk.id) &&
                            (other_gce.neighbourFace == j)) {
                            ok = true;
                        }
                    }
                    if (!ok) {
                        string msg = format("FullFaceCopy for local blk_id=%d face=%d", blk.id, j);
                        msg ~= format(" is not correctly paired with other block id=%d face=%d.",
                                      other_blk.id, my_gce.neighbourFace);
                        writeln(msg);
                        anyBlockFail = true;
                    }
                }
            }
        }
    }
    // [TODO] mpi version: reduce anyBlockFail flag
    if (anyBlockFail) {
        throw new NewtonKrylovException("Failed at initialisation stage during full-face boundary data exchange.");
    }
}

void initMappedCellDataExchange()
{
    // Serial loops follow because the cell-mapping function searches across
    // all blocks local to the process.
    // Also, there are several loops because the MPI communication,
    // if there is any, needs to be done in phases of posting of non-blocking reads,
    // followed by all of the sends and then waiting for all requests to be filled.
    //
    bool anyBlockFail = false;
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellMappedCellCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce && (mygce.check_cell_mapping() != 0)) { anyBlockFail = true; }
            }
        }
    }
    // [TODO] mpi version: reduce anyBlockFail flag
    if (anyBlockFail) {
        throw new NewtonKrylovException("Failed at initialisation stage during locating mapped-cell boundaries.");
    }
    
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase0(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase1(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase2(); }
            }
        }
    }
}

void initGhostCellGeometry()
{
    exchange_ghost_cell_geometry_data();
}

void initLeastSquaresStencils()
{
    foreach (blk; localFluidBlocks) blk.compute_least_squares_setup(0);
}

void initMLPlimiter()
{
    foreach (blk; localFluidBlocks) {
        auto ublock = cast(UFluidBlock) blk;
        if (ublock) { ublock.build_cloud_of_cell_references_at_each_vertex(); }
    }
}

void initMasterLuaState()
{
    auto L = GlobalConfig.master_lua_State;
    lua_pushboolean(L, GlobalConfig.in_mpi_context);
    lua_setglobal(L, "in_mpi_context");
    lua_pushnumber(L, GlobalConfig.mpi_size);
    lua_setglobal(L, "mpi_size");
    lua_pushnumber(L, GlobalConfig.mpi_rank_for_local_task);
    lua_setglobal(L, "mpi_rank_for_local_task");
    lua_pushboolean(L, GlobalConfig.is_master_task);
    lua_setglobal(L, "is_master_task");
    push_array_to_Lua(L, GlobalConfig.localFluidBlockIds, "localFluidBlockIds");
    // [TODO] think about user_pad -- does it have a use case in steady-state?
}

void initCornerCoordinates()
{
    synchronize_corner_coords_for_all_blocks();
}

void initWallDistances()
{
    compute_wall_distances();
}

/*---------------------------------------------------------------------
 * Main iteration algorithm
 *---------------------------------------------------------------------
 */

void performNewtonKrylovUpdates(int snapshotStart, int maxCPUs, int threadsPerMPITask)
{
    alias cfg = GlobalConfig;
    string jobName = cfg.base_file_name;
    cfg.print_count = 1;

    if (cfg.verbosity_level > 1) writeln("Read N-K config file.");
    JSONValue jsonData = readJSONfile(lmrCfg.nkCfgFile);
    nkCfg.readValuesFromJSON(jsonData);
    // Allocate space and configure phases
    nkPhases.length = nkCfg.numberOfPhases;
    foreach (i, ref phase; nkPhases) {
        string key = "NewtonKrylovPhase_" ~ to!string(i);
        phase.readValuesFromJSON(jsonData[key]);
    }

    int nWrittenSnapshots;
    int startStep = 1;
    bool finalStep = false;
    double cfl;
    double dt;
    bool updatePreconditionerThisStep = false;
    CFLSelector cflSelector;
    
    /*----------------------------------------------
     * Initialisation
     *----------------------------------------------
     */
    setAndReportThreads(maxCPUs, threadsPerMPITask);
    if (nkCfg.usePreconditioner) initPreconditioner();
    size_t nConserved = cfg.cqi.n;
    referenceResiduals = new ConservedQuantities(nConserved);
    currentResiduals = new ConservedQuantities(nConserved);
    scale = new ConservedQuantities(nConserved);
    if (cfg.is_master_task) {
        initialiseDiagnosticsFile();
    }
    allocateGlobalGMRESWorkspace();
    foreach (blk; localFluidBlocks) {
        blk.allocate_GMRES_workspace(nkCfg.maxLinearSolverIterations);
    }
    /* solid blocks don't work just yet.
    allocate_global_solid_workspace();
    foreach (sblk; localSolidBlocks) {
        sblk.allocate_GMRES_workspace();
    }
    */

    /*
    debug {
        writeln("DEBUG: performNewtoKrylovUpdate()");
        writeln("       Initialisation done.");
    }
    */
    
    // Look for global CFL schedule and use to set CFL
    if (nkCfg.cflSchedule.length > 0) {
        foreach (i, startRamp; nkCfg.cflSchedule[0 .. $-1]) {
            if (startStep >= startRamp.step) {
                auto endRamp = nkCfg.cflSchedule[i+1];
                cflSelector = new LinearRampCFL(startRamp.step, endRamp.step, startRamp.cfl, endRamp.cfl);
                break;
            }
        }
        // Or check we aren't at end of cfl schedule
        auto lastEntry = nkCfg.cflSchedule[$-1];
        if (startStep >= lastEntry.step) {
            // Set a flat CFL beyond limit of scheule.
            cflSelector = new LinearRampCFL(lastEntry.step, nkCfg.maxNewtonSteps, lastEntry.cfl, lastEntry.cfl);
        }
    }
    
    if (snapshotStart > 0) {
        /*----------------------------------------------
         * Special actions on restart
         *----------------------------------------------
         *  + extract restart information
         *  + read in reference residuals from file
         *  + determine how many snapshots have already been written
         *  + determine phase and set as active phase
         */
        extractRestartInfoFromTimesFile(jobName);
        referenceGlobalResidual = setReferenceResidualsFromFile();
        referenceResidualsAreSet = true;
        nWrittenSnapshots = determineNumberOfSnapshots();
        RestartInfo restart = snapshots[snapshotStart];
        startStep = restart.step + 1;
        globalResidual = restart.globalResidual;
        prevGlobalResidual = globalResidual;
        // Determine phase
        foreach (phase, phaseStep; nkCfg.phaseChangesAtSteps) {
            if (startStep < phaseStep) {
                setPhaseSettings(phase); 
                break;
            }
        }
        // end condition when step is past final phase
        if (startStep >= nkCfg.phaseChangesAtSteps[$-1]) {
            auto finalPhase = nkCfg.phaseChangesAtSteps.length;
            setPhaseSettings(finalPhase);
        }
        if (activePhase.useAutoCFL) {
            cflSelector = new ResidualBasedAutoCFL(activePhase.autoCFLExponent, activePhase.maxCFL, activePhase.thresholdRelativeResidualForCFLGrowth);
            cfl = restart.cfl;
        }
        else { // Assume we have a global (phase-independent) schedule
            cfl = cflSelector.nextCFL(-1.0, startStep, -1.0, -1.0, -1.0);
        }
    }
    else {
        // On fresh start, the phase setting must be at 0
        setPhaseSettings(0);
        if (activePhase.useAutoCFL) {
            cflSelector = new ResidualBasedAutoCFL(activePhase.autoCFLExponent, activePhase.maxCFL, activePhase.thresholdRelativeResidualForCFLGrowth);
            cfl = activePhase.startCFL;
        }
        else { // Assume we have a global (phase-independent) schedule
            cfl = cflSelector.nextCFL(-1.0, startStep, -1.0, -1.0, -1.0);
        }
    }

    // Start timer right at beginning of stepping.
    auto wallClockStart = Clock.currTime();
    double wallClockElapsed;
    int numberBadSteps = 0;
    bool startOfNewPhase = false;

    foreach (step; startStep .. nkCfg.maxNewtonSteps+1) {
        /*
        debug {
            writeln("DEBUG: performNewtoKrylovUpdate()");
            writefln("       STEP %2d ", step);
        }
        */
        
        /*---
         * 0. Check for any special actions based on step to perform at START of step
         *---
         *    a. change of phase
         *    b. limiter freezing
         *    c. set the timestep
         *    d. set flag on preconditioner
         */
        residualsUpToDate = false;
        // 0a. change of phase 
        size_t currentPhase = countUntil(nkCfg.phaseChangesAtSteps, step);
        startOfNewPhase = false;
        if (currentPhase != -1) { // start of new phase detected
            startOfNewPhase = true;
            setPhaseSettings(currentPhase);
            if (activePhase.useAutoCFL) {
                // When we change phase, we reset CFL to user's selection if using auto CFL.
                cfl = activePhase.startCFL;
            }
        }

        // 0b. Check if we need to freeze limiter
        if (step == nkCfg.freezeLimiterAtStep) {
            // We need to compute the limiter a final time before freezing it.
            // This is achieved via evalRHS
            if (cfg.frozen_limiter == false) {
                evalResidual(0);
            }
            cfg.frozen_limiter = true;
        }

        if (step < nkCfg.freezeLimiterAtStep) {
            // Make sure that limiter is on in case it has been set
            // to frozen during a Jacobian evaluation.
            cfg.frozen_limiter = true;
        }

        // 0c. Set the timestep for this step
        dt = setDtInCells(cfl);

        // 0d. determine if we need to update preconditioner

        if (step == startStep || startOfNewPhase || (step % activePhase.stepsBetweenPreconditionerUpdate) == 0) {
            updatePreconditionerThisStep = true;
        }
        else {
            updatePreconditionerThisStep = false;
        }
        /*---
         * 1. Perforn Newton update
         *---
         */
        /*
        debug {
            writeln("DEBUG: performNewtoKrylovUpdate()");
            writefln("       about to solve for Newton step.", step);
        }
        */

        globalResidual = solveNewtonStep(updatePreconditionerThisStep);

        /*
        debug {
            writeln("DEBUG: performNewtoKrylovUpdate()");
            writefln("       Newton step done.", step);
        }
        */


        double omega = nkCfg.usePhysicalityCheck ? determineRelaxationFactor() : 1.0;
        if (omega >= nkCfg.minRelaxationFactor) {
            // Things are good. Apply omega-scaled update and continue on.
            // We think??? If not, we bail at this point.
            try {
                applyNewtonUpdate(omega);
            }
            catch (NewtonKrylovException e) {
                // We need to bail out at this point.
                // User can probably restart with less aggressive CFL schedule.
                if (GlobalConfig.is_master_task) {
                    writeln("Update failure in Newton step.");
                    writefln("step= %d, CFL= %e, dt= %e, global-residual= %e ", step, cfl, dt, globalResidual);
                    writeln("Error message from failed update:");
                    writefln("%s", e.msg);
                    writeln("You might be able to try a smaller CFL.");
                    writeln("Bailing out!");
                    exit(1);
                }
            }
            // Since things are good, we can select a new CFL
            if (activePhase.useAutoCFL) {
                // We need to be careful in the early steps with the auto CFL.
                // On step 1, we have no previous residual, so we can't make an adjustment.
                // Also, we need to check on a residual drop, but this makes no sense
                // until the reference residuals are established.
                if (step > startStep && step > nkCfg.numberOfStepsForSettingReferenceResiduals) {
                    cfl = cflSelector.nextCFL(cfl, step, globalResidual, prevGlobalResidual, globalResidual/referenceGlobalResidual);
                }
            }
            else {
                cfl = cflSelector.nextCFL(-1.0, step, -1.0, -1.0, -1.0);
            }
            numberBadSteps = 0;
        }
        else if (omega < nkCfg.minRelaxationFactor && activePhase.useAutoCFL) {
            numberBadSteps++;
            if (numberBadSteps == nkCfg.maxConsecutiveBadSteps) {
                writeln("Too many consecutive bad steps while trying to update flow state.\n");
                writefln("Number of bad steps = %d", numberBadSteps);
                writefln("Last attempted CFL = %e", cfl);
                writeln("Bailing out!");
                exit(1);
            }

            // Do NOT apply update. Simply change CFL for next step.
            cfl *= nkCfg.cflReductionFactor;
            if (cfl <= nkCfg.cflMin) {
                writeln("The CFL has been reduced due a bad step, but now it has dropped below the minimum allowable CFL.");
                writefln("current cfl = %e  \t minimum allowable cfl = %e", cfl, nkCfg.cflMin);
                writeln("Bailing out!");
                exit(1);
            }
            
            // Return flow states to their original state for next attempt.
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (cell; blk.cells) {
                    cell.decode_conserved(0, 0, 0.0);
                }
            }
        }
        else {
            if (GlobalConfig.is_master_task) {
                writeln("WARNING: relaxation factor for Newton update is very small.");
                writefln("step= %d, relaxation factor= %f");
                writeln("Bailing out!");
                exit(1);
            }
        }
        
        /*---
         * 2. Post-update actions
         *---
         * Here we need to do some house-keeping and see if we continue with iterations.
         */
        /*
        debug {
            writeln("DEBUG: performNewtoKrylovUpdate()");
            writefln("       post-update actions on step=  %d", step);
        }
        */
        
        /*----
         * 2a. Search for reference residuals if needed
         *----
         */
        if (!referenceResidualsAreSet) {
            referenceGlobalResidual = fmax(referenceGlobalResidual, globalResidual);
            if (!residualsUpToDate) {
                computeResiduals(currentResiduals);
                residualsUpToDate = true;
            }
            foreach (ivar; 0 .. nConserved) {
                referenceResiduals[ivar] = fmax(referenceResiduals[ivar], currentResiduals[ivar]);
            }
            if (step == nkCfg.numberOfStepsForSettingReferenceResiduals) {
                referenceResidualsAreSet = true;
                if (GlobalConfig.is_master_task) {
                    writeln("*************************************************************************");
                    writeln("*");
                    writefln("*  After first %d steps, reference residuals have been set.", nkCfg.numberOfStepsForSettingReferenceResiduals);
                    writefln("*  Reference global residual: %.12e", referenceGlobalResidual);
                    writeln("*");
                    writeln("*  Reference residuals for each conservation equation:");
                    foreach (ivar; 0 .. nConserved) {
                        writefln("* %12s: %.12e", cfg.cqi.names[ivar], referenceResiduals[ivar].re);
                    }
                    writeln("*");
                    writeln("*************************************************************************\n");
                }
            }
        }

        // We can now set previous residual in preparation for next step.
        prevGlobalResidual = globalResidual;
        
        /*---
         * 2b. Stopping checks.
         *---
         */
        if (step == nkCfg.maxNewtonSteps) {
            finalStep = true;
            if (cfg.is_master_task) {
                writeln("STOPPING: Reached maximum number of steps.");
            }
        }
        if (!activePhase.ignoreStoppingCriteria) {
            if (globalResidual <= nkCfg.stopOnAbsoluteResidual) {
                finalStep = true;
                if (cfg.is_master_task) {
                    writeln("STOPPING: The absolute global residual is below target value.");
                    writefln("         current global residual= %.6e  target value= %.6e", globalResidual, nkCfg.stopOnAbsoluteResidual);
                }
            }
            if ((globalResidual/referenceGlobalResidual) <= nkCfg.stopOnRelativeResidual) {
                finalStep = true;
                if (cfg.is_master_task) {
                    writeln("STOPPING: The relative global residual is below target value.");
                    writefln("         current residual= %.6e  target value= %.6e", (globalResidual/referenceGlobalResidual), nkCfg.stopOnRelativeResidual);
                }
            }
        }
        // [TODO] Add in a halt_now condition.
        /*---
         * 2c. Reporting (to files and screen)
         *---
         */
        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();
        if (((step % nkCfg.stepsBetweenDiagnostics) == 0) || finalStep) {
            writeDiagnostics(step, dt, cfl, wallClockElapsed, omega, residualsUpToDate);
        }

        if (((step % nkCfg.stepsBetweenSnapshots) == 0) || finalStep) {
            writeSnapshot(step, nWrittenSnapshots);
        }
        
        
        // [TODO] Write loads. We only need one lot of loads.
        // Any intermediate loads before steady-state have no physical meaning.
        // They might have some diagnostic purpose?

        // Reporting to screen on progress.
        if ( ((step % cfg.print_count) == 0) || finalStep ) {
            printStatusToScreen(step, cfl, dt, wallClockElapsed, residualsUpToDate);
        }

        if (finalStep) break;
        
    }
}



/*---------------------------------------------------------------------
 * Auxiliary functions related to initialisation of stepping
 *---------------------------------------------------------------------
 */

void setAndReportThreads(int maxCPUs, int threadsPerMPITask)
{
    alias cfg = GlobalConfig;

    // 1. Check we aren't using more task threads than blocks
    int extraThreadsInPool;
    auto nBlocksInThreadParallel = localFluidBlocks.length;
    version(mpi_parallel) {
        extraThreadsInPool = min(threadsPerMPITask-1, nBlocksInThreadParallel-1);
    } else {
        extraThreadsInPool = min(maxCPUs-1, nBlocksInThreadParallel-1);
    }
    defaultPoolThreads(extraThreadsInPool);

    // 2. Report out the thread configuration for run-time
    version(mpi_parallel) {
        writefln("MPI-task %d : running with %d threads.", cfg.mpi_rank_for_local_task, extraThreadsInPool+1);
    }
    else {
        writefln("Single process running with %d threads.", extraThreadsInPool+1); // +1 for main thread.
    }
}

void initPreconditioner()
{
    if (nkCfg.usePreconditioner) {
        evalResidual(0);
        // initialize the flow Jacobians used as local precondition matrices for GMRES
        final switch (nkCfg.preconditioner) {
        case PreconditionerType.jacobi:
            foreach (blk; localFluidBlocks) { blk.initialize_jacobian(-1, nkCfg.preconditionerPerturbation); }
            break;
        case PreconditionerType.ilu:
            foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, nkCfg.preconditionerPerturbation, nkCfg.iluFill); }
            break;
        case PreconditionerType.sgs:
            foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, nkCfg.preconditionerPerturbation); }
            break;
        case PreconditionerType.lusgs:
            // do nothing
            break;
        case PreconditionerType.diagonal:
            // do nothing
            break; 
        } // end switch
    }
}

void extractRestartInfoFromTimesFile(string jobName)
{
    size_t nConserved = GlobalConfig.cqi.n;
    RestartInfo restartInfo = RestartInfo(nConserved);
    // Start reading the times file, looking for the snapshot index
    auto timesFile = File("./config/" ~ jobName ~ ".times");
    auto line = timesFile.readln().strip();
    while (line.length > 0) {
        if (line[0] != '#') {
            // Process a non-comment line.
            auto tokens = line.split();
            auto idx = to!int(tokens[0]);
            restartInfo.pseudoSimTime = to!double(tokens[1]);
            restartInfo.dt = to!double(tokens[2]);
            restartInfo.cfl = to!double(tokens[3]);
            restartInfo.step = to!int(tokens[4]);
            restartInfo.globalResidual = to!double(tokens[5]);
            size_t startIdx = 6;
            foreach (ivar; 0 .. nConserved) {
                restartInfo.residuals[ivar] = to!double(tokens[startIdx+ivar]);
            }
            snapshots ~= restartInfo;
        }
        line = timesFile.readln().strip();
    }
    timesFile.close();
    return;
}

double setReferenceResidualsFromFile()
{
    double refGlobalResidual;
    size_t nConserved = GlobalConfig.cqi.n;
    
    auto refResid = File(refResidFname, "r");
    auto line = refResid.readln().strip();
    auto tokens = line.split();
    refGlobalResidual = to!double(tokens[0]);
    size_t startIdx = 1;
    foreach (ivar; 0 .. nConserved) {
        referenceResiduals[ivar] = to!double(tokens[startIdx+ivar]);
    }
    refResid.close();
    return refGlobalResidual;
}

int determineNumberOfSnapshots()
{
    string jobName = GlobalConfig.base_file_name;
    int nWrittenSnapshots = 0;
    auto timesFile = File("./config/" ~ jobName ~ ".times");
    auto line = timesFile.readln().strip();
    while (line.length > 0) {
        if (line[0] != '#') {
            nWrittenSnapshots++;
        }
        line = timesFile.readln().strip();
    }
    timesFile.close();
    nWrittenSnapshots--; // We don't count the initial solution as a written snapshot
    return nWrittenSnapshots; 
}

void allocateGlobalGMRESWorkspace()
{
    size_t m = to!size_t(nkCfg.maxLinearSolverIterations);
    g0.length = m+1;
    g1.length = m+1;
    h.length = m+1;
    hR.length = m+1;
    H0 = new Matrix!number(m+1, m);
    H1 = new Matrix!number(m+1, m);
    Gamma = new Matrix!number(m+1, m+1);
    Q0 = new Matrix!number(m+1, m+1);
    Q1 = new Matrix!number(m+1, m+1);
}

/*---------------------------------------------------------------------
 * Auxiliary functions related to iteration algorithm
 *---------------------------------------------------------------------
 */

void setPhaseSettings(size_t phase)
{
    activePhase = nkPhases[phase];
    foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(activePhase.residualInterpolationOrder);
}

void computeResiduals(ref ConservedQuantities residuals)
{
    size_t nConserved = GlobalConfig.cqi.n;

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.residuals[] = blk.cells[0].dUdt[0][];
        foreach (ivar; 0 .. nConserved) blk.residuals[ivar] = fabs(blk.residuals[ivar]);

        foreach (cell; blk.cells) {
            foreach (ivar; 0 .. nConserved) blk.residuals[ivar] = fmax(blk.residuals[ivar], fabs(cell.dUdt[0][ivar]));
        }
    }
    // Do next bit in serial to reduce information to current thread
    residuals[] = localFluidBlocks[0].residuals[];
    foreach (blk; localFluidBlocks) {
        foreach (ivar; 0 .. nConserved) residuals[ivar] = fmax(residuals[ivar], blk.residuals[ivar]);
    }
    // and for MPI a reduce onto master rank
    version(mpi_parallel) {
        foreach (ivar; 0 .. nConserved) {
            if (GlobalConfig.is_master_task) {
                MPI_Reduce(MPI_IN_PLACE, &(residuals[ivar].re), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            }
            else {
                MPI_Reduce(&(residuals.vec[ivar].re), &(residuals[ivar].re), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            }
        }
    }
}

/**
 * This function sets a dt in all cells and finds the minimum across the domain.
 *
 * If we are using local timestepping, the local value will (likely) differ in each cell.
 * If we are using global timestepping, we explicitly the local values to a consistent
 * global minimum.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-08
 */
double setDtInCells(double cfl)
{
    double signal;
    double dt;
    bool inviscidCFLOnly = nkCfg.inviscidCFLOnly;
    bool firstCell;
    foreach (blk; parallel(localFluidBlocks,1)) {
        firstCell = true;
        foreach (i, cell; blk.cells) {
            if (inviscidCFLOnly) {
                // calculate the signal using the maximum inviscid wave speed
                // ref. Computational Fluid Dynamics: Principles and Applications, J. Blazek, 2015, pg. 175
                // Note: this approximates cell width dx as the cell volume divided by face area
                signal = 0.0;
                foreach (f; cell.iface) {
                    number un = fabs(f.fs.vel.dot(f.n));
                    number signal_f = (un+f.fs.gas.a)*f.area[0];
                    signal += signal_f.re;
                }
                signal *= (1.0/cell.volume[0].re);
            }
            else {
                // use the default signal frequency routine from the time-accurate code path
                signal = cell.signal_frequency();
            }
            cell.dt_local = cfl / signal;
            if (firstCell) {
                blk.dtMin = cell.dt_local;
                firstCell = false;
            }
            else {
                blk.dtMin = fmin(blk.dtMin, cell.dt_local);
            }
        }
    }
    // Find smallest dt across these local blocks in serial search
    dt = localFluidBlocks[0].dtMin;
    foreach (blk; localFluidBlocks) {
        dt = fmin(dt, blk.dtMin);
    }

    version(mpi_parallel) {
        // Find smallest dt globally if using distributed memory.
        MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }

    // When using global timestepping, make all dt_local consistent.
    if (!nkCfg.useLocalTimestep) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (cell; blk.cells) cell.dt_local = dt;
        }
    }
    
    return dt;
}

/**
 * Set dt_local in cells based on CFL.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-08
 */
void setDtLocalInCells(double cfl)
{
    double signal;
    bool inviscidCFLOnly = nkCfg.inviscidCFLOnly;
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (cell; blk.cells) {
            signal = cell.signal_frequency();
            if (inviscidCFLOnly) signal = cell.signal_hyp.re;
            cell.dt_local = cfl / signal;
        }
    }
}


/*---------------------------------------------------------------------
 * Mixins to handle shared-memory dot product and norm
 *---------------------------------------------------------------------
 */

string dotOverBlocks(string dot, string A, string B)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.dotAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.dotAcc += blk.`~A~`[k].re*blk.`~B~`[k].re;
   }
}
`~dot~` = 0.0;
foreach (blk; localFluidBlocks) `~dot~` += blk.dotAcc;`;

}

string norm2OverBlocks(string norm2, string blkMember)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.normAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.normAcc += blk.`~blkMember~`[k].re*blk.`~blkMember~`[k].re;
   }
}
`~norm2~` = 0.0;
foreach (blk; localFluidBlocks) `~norm2~` += blk.normAcc;
`~norm2~` = sqrt(`~norm2~`);`;

}


/**
 * This function solves a linear system to provide a Newton step for the flow field.
 *
 * The particular linear solver used here is GMRES. It uses right-preconditioning,
 * scaling and possibly restarts.
 *
 * The linear solver is baked in-place here because the "sytem" is so closely
 * tied to the the flow field. In other words, we require several flow field
 * residual evaluations. This coupling between linear solver and fluid updates
 * makes it difficult to extract as a stand-alone routine.
 *
 * NOTE: GMRES algorithm does not need to compute the global residual.
 * However, it's convenient and efficient to do this computation here as a side effect.
 * The reason is because we've done all the work of packing the residual vector, R.
 * While the data is current, we can do the parallel computation of norm on that vector.
 *
 * Algorithm overview:
 *    0. Preparation for iterations.
 *       0a) Evaluate dUdt and store in R.
 *       0b) Determine scale factors.
 *       0c) Compute global residual. (Efficiency measure)
 *       0d) Compute scaled r0.
 *       0e) Compute preconditioner matrix
 *    1. Start loop on restarts
 *       1a) Perform inner iterations.
 *       1b) Check tolerance for convergence
 *    2. Post-iteration actions
 *       2a) Unscale values
 *     
 * Authors: RJG and KAD
 * Date: 2022-03-02
 * History:
 *    2022-03-02  Major clean-up as part of refactoring work.
 */
double solveNewtonStep(bool updatePreconditionerThisStep)
{
    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" entered function.");
    }
    */
            
    alias cfg = GlobalConfig;
    size_t nConserved = cfg.cqi.n;

    bool isConverged = false;
    /*---
     * 0. Preparation for iterations.
     *---
     */
    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" calling evalResidual.");
    }
    */

    evalResidual(0);

    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" evalResidual done.");
    }
    */
    
    setResiduals();

    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" calling computeGlobalResidual.");
    }
    */
    
    double globalResidual = computeGlobalResidual();

    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" computeGlobalResidual done.");
    }
    */

    
    determineScaleFactors(scale);
    // r0 = A*x0 - b
    compute_r0(scale);

    // beta = ||r0||
    number beta = computeLinearSystemResidual();
    number beta0 = beta; // Store a copy as initial residual
                         // because we will look for relative drop in residual
                         // compared to this.
    // v = r0/beta
    prepareKrylovSpace(beta);

    auto targetResidual = activePhase.linearSolveTolerance * beta;
    if (nkCfg.usePreconditioner && updatePreconditionerThisStep) {
        computePreconditioner();
    }

    /*---
     * 1. Outer loop of restarted GMRES
     *---
     */

    size_t r_break;
    int maxIterations = nkCfg.maxLinearSolverIterations;
    // We add one here because input is to do with number of *restarts*.
    // We need at least one attempt (+1) plus the number of restarts chosen.
    int nAttempts = nkCfg.maxLinearSolverRestarts + 1;
    int iterationCount;

    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" starting restarted GMRES iterations.");
    }
    */


    foreach (r; 0 .. nAttempts) {
        // Initialise some working arrays and matrices for this step.
        g0[] = to!number(0.0);
        g1[] = to!number(0.0);
        H0.zeros();
        H1.zeros();
        // Set first residual entry.
        g0[0] = beta;

        // Delegate inner iterations
        /*
        debug {
            writeln("DEBUG: solveNewtonStep()");
            writeln(" calling performIterations.");
        }
        */
        
        isConverged = performIterations(maxIterations, beta0, targetResidual, scale, iterationCount);
        int m = iterationCount;

        /*
        debug {
            writeln("DEBUG: solveNewtonStep()");
            writeln(" done performIterations.");
        }
        */
        
        // At end H := R up to row m
        //        g := gm up to row m
        upperSolve!number(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot!number(blk.V, blk.nvars, m, blk.g1, blk.zed);
            unscaleVector(blk.zed, nConserved);
        }

        // Prepare dU values (for Newton update)
        if (nkCfg.usePreconditioner) {
            // Remove preconditioner effect from values.
            removePreconditioning();
        }
        else {
            foreach(blk; parallel(localFluidBlocks,1)) {
                blk.dU[] = blk.zed[];
            }
        }

        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) blk.dU[k] += blk.x0[k];
        }

        if (isConverged || (r == nkCfg.maxLinearSolverRestarts)) {
            // We are either converged, or
            // we've run out of restart attempts.
            // In either case, we can leave now.
            r_break = r;
            break;
        }

        // If we get here, we need to prepare for next restart.
        // This requires setting x0[] and r0[].
        // We'll compute r0[] using the approach of Fraysee et al. (2005)
        foreach (blk; parallel(localFluidBlocks, 1)) {
            blk.x0[] = blk.dU[];
        }

        foreach (blk; localFluidBlocks) copy(Q1, blk.Q1);
        // Set all values in g0 to 0.0 except for final (m+1) value
        foreach (i; 0 .. m) g0[i] = 0.0;
        foreach (blk; localFluidBlocks) blk.g0[] = g0[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.Q1, m, m+1, blk.g0, blk.g1);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.V, blk.nvars, m+1, blk.g1, blk.r0);
        }

        mixin(dotOverBlocks("beta", "r0", "r0"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
        }
        beta = sqrt(beta);

        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.r0[k]/beta;
                blk.V[k,0] = blk.v[k];
            }
        }
    }

    /*
    debug {
        writeln("DEBUG: solveNewtonStep()");
        writeln(" done restarted GMRES iterations.");
    }
    */
    
    // Set some information before leaving. This might be used in diagnostics file.
    gmresInfo.nRestarts = to!int(r_break);
    gmresInfo.initResidual = beta0.re;
    gmresInfo.finalResidual = beta.re;
    gmresInfo.iterationCount = iterationCount;
    
    // NOTE: global residual is value at START, before applying Newton update.
    // The caller applies the update, so we can't compute the new global
    // residual in here.
    return globalResidual;
}

/**
 * Copy values from dUdt into R.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
void setResiduals()
{
    size_t nConserved = GlobalConfig.cqi.n;
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t startIdx = 0;
        foreach (i, cell; blk.cells) {
            blk.R[startIdx .. startIdx+nConserved] = cell.dUdt[0][0 .. nConserved];
            startIdx += nConserved;
        }
    }
}

/**
 * Determine scale factors per conserved quantity.
 *
 * The scale factors are typically the maximum rate of change found globally
 * (over all cells) per each conserved quantity. However, we place some gaurds
 * on determining those scales when rates of change are small.
 * 
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
void determineScaleFactors(ref ConservedQuantities scale)
{
    scale[] = to!number(0.0);
    // First do this for each block.
    size_t nConserved = GlobalConfig.cqi.n;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.maxRate[] = to!number(0.0);
        foreach (i, cell; blk.cells) {
            foreach (ivar; 0 .. nConserved) {
                blk.maxRate[ivar] = fmax(blk.maxRate[ivar], fabs(cell.dUdt[0][ivar]));
            }
        }
    }
    // Next, reduce that maxR information across all blocks and processes
    foreach (blk; localFluidBlocks) {
        foreach (ivar; 0 .. nConserved) {
            scale[ivar] = fmax(scale[ivar], blk.maxRate[ivar]);
        }
    }
    // In distributed memory, reduce max values and make sure everyone has a copy.
    version(mpi_parallel) {
        foreach (ivar; 0 .. nConserved) {
            MPI_Allreduce(MPI_IN_PLACE, &(scale[ivar].re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
    }

    // Use a guard on scale values if they get small
    foreach (ivar; 0 .. nConserved) {
        scale[ivar] = fmax(scale[ivar], minScaleFactor);
    }
    
    // Value is presently maxRate. Store as scale = 1/maxRate.
    foreach (ivar; 0 .. nConserved) {
        scale[ivar] = 1./scale[ivar];
    }
    
}


/**
 * Compute the global residual based on vector R.
 *
 * For certain turbulence models, we scale the contribution in the global residual
 * so that certain quantities do not dominate this norm. For example, the turbulent
 * kinetic energy in the k-omega model has its contribution scaled down.
 *
 * Authors: KAD and RJG
 * Date: 2022-03-02
 */
double computeGlobalResidual()
{
    double globalResidual;
    mixin(dotOverBlocks("globalResidual", "R", "R"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    globalResidual = sqrt(globalResidual);

    return globalResidual;
}


/**
 * Compute the initial scaled residual vector for GMRES method.
 *
 * r0 = b - A*x0
 *
 * With x0 = [0] (as is common), r0 = b = R
 * However, we apply scaling also at this point.
 *
 * Authors: KAD and RJG
 * Date: 2022-03-02
 */
void compute_r0(ConservedQuantities scale)
{
    size_t nConserved = GlobalConfig.cqi.n;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.x0[] = to!number(0.0);
        foreach (i; 0 .. blk.r0.length) {
            size_t ivar = i % nConserved;
            blk.r0[i] = scale[ivar]*blk.R[i];
        }
    }
}

/**
 * Compute the residual of the linear system from r0.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
number computeLinearSystemResidual()
{
    number beta;
    mixin(dotOverBlocks("beta", "r0", "r0"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
    }
    beta = sqrt(beta);
    return beta;
}

/**
 * Prepare the Krylov space for iterating.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
void prepareKrylovSpace(number beta)
{
    g0[0] = beta;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.v[] = (1./beta)*blk.r0[];
        foreach (k; 0 .. blk.nvars) {
            blk.V[k,0] = blk.v[k];
        }
    }
}

/**
 * Computes the preconditioner for use in GMRES solver.
 *
 * Authors: KAD and RJG
 * Date: 2022-03-02
 */
void computePreconditioner()
{
    size_t nConserved = GlobalConfig.cqi.n;
    
    final switch (nkCfg.preconditioner) {
    case PreconditionerType.lusgs:
        // do nothing
        break;
    case PreconditionerType.diagonal:
        goto case PreconditionerType.sgs;
    case PreconditionerType.jacobi:
        goto case PreconditionerType.sgs;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.evaluate_jacobian();
            blk.flowJacobian.augment_with_dt(blk.cells, nConserved);
            nm.smla.invert_block_diagonal(blk.flowJacobian.local, blk.flowJacobian.D, blk.flowJacobian.Dinv, blk.cells.length, nConserved);
        }
        break;
    case PreconditionerType.ilu:
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.evaluate_jacobian();
            blk.flowJacobian.augment_with_dt(blk.cells, nConserved);
            nm.smla.decompILU0(blk.flowJacobian.local);
        }
        break;
    }
}

/*---------------------------------------------------------------------
 * Mixins to scale and unscale vectors
 *---------------------------------------------------------------------
 */

void scaleVector(ref number[] vec, size_t nConserved)
{
    foreach (k, ref v; vec) {
        size_t ivar = k % nConserved;
        v *= scale[ivar];
    }
}

void unscaleVector(ref number[] vec, size_t nConserved)
{
    foreach (k, ref v; vec) {
        size_t ivar = k % nConserved;
        v /= scale[ivar];
    }
}


/**
 * Perform GMRES iterations to fill Krylov subspace and get solution estimate.
 *
 * Returns boolean indicating converged or not.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
bool performIterations(int maxIterations, number beta0, number targetResidual,
                       ref ConservedQuantities scale, ref int iterationCount)
{

    /*
    debug {
        writeln("DEBUG: performIterations()");
        writeln(" entered function.");
    }
    */
    
    alias cfg = GlobalConfig;

    bool isConverged = false;
    size_t nConserved = cfg.cqi.n;

    /*
    debug {
        writeln("DEBUG: performIterations()");
        writeln(" starting Krylov iterations.");
    }
    */
        
    foreach (j; 0 .. maxIterations) {
        iterationCount = j+1;

        /*
        debug {
            writeln("DEBUG: performIterations()");
            writeln(" scaling vector.");
        }
        */

        
        // 1. Unscale v
        // v is scaled earlier when r0 copied in.
        // However, to compute Jv via Frechet, we will need
        // unscaled values.
        foreach (blk; parallel(localFluidBlocks,1)) {
            unscaleVector(blk.v, nConserved);
        }

        /*
        debug {
            writeln("DEBUG: performIterations()");
            writeln(" apply preconditioning.");
        }
        */
        
        // 2. Apply preconditioning (if requested)
        if (nkCfg.usePreconditioner) {
            applyPreconditioning();
        }
        else {
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.zed[] = blk.v[];
            }
        }

        /*
        debug {
            writeln("DEBUG: performIterations()");
            writeln(" include 1/dt term.");
        }
        */
        
        // 3. Jacobian-vector product
        // 3a. Prepare w vector with 1/dt term.
        //      (I/dt)(P^-1)v
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t startIdx = 0;
            foreach (cell; blk.cells) {
                number dtInv = 1.0/cell.dt_local;
                blk.w[startIdx .. startIdx + nConserved] = dtInv*blk.zed[startIdx .. startIdx + nConserved];
                startIdx += nConserved;
            }
        }

        // 3b. Determine perturbation size, sigma
        // Kyle's experiments show one needs to recompute on every step
        // if using the method to estimate a perturbation size based on
        // vector 'v'.
        double sigma;
        version (complex_numbers) {
            // For complex-valued Frechet derivative, a very small perturbation
            // works well (almost) all the time.
            sigma = nkCfg.frechetDerivativePerturbation;
        }
        else {
            // For real-valued Frechet derivative, we may need to attempt to compute
            // a perturbation size.
            sigma =  (nkCfg.frechetDerivativePerturbation < 0.0) ? computePerturbationSize() : nkCfg.frechetDerivativePerturbation;
        }

        /*
        debug {
            writeln("DEBUG: performIterations()");
            writeln(" Jz calc.");
        }
        */

        // 3b. Evaluate Jz and place result in z
        evalJacobianVectorProduct(sigma);

        // 3c. Complete the calculation of w
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars)  blk.w[k] = blk.w[k] - blk.zed[k];
            scaleVector(blk.w, nConserved);
        }

        /*
        debug {
            writeln("DEBUG: performIterations()");
            writeln(" remainder of GMRES.");
        }
        */

        
        // 4. The remainder of the algorithm looks a lot like any standard
        // GMRES implementation (for example, see smla.d)
        foreach (i; 0 .. j+1) {
            foreach (blk; parallel(localFluidBlocks,1)) {
                // Extract column 'i'
                foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.V[k,i];
            }
            number H0_ij;
            mixin(dotOverBlocks("H0_ij", "w", "v"));
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
            }
            H0[i,j] = H0_ij;
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k];
            }
        }
        number H0_jp1j;
        mixin(dotOverBlocks("H0_jp1j", "w", "w"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(H0_jp1j.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(H0_jp1j.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
        }
        H0_jp1j = sqrt(H0_jp1j);
        H0[j+1,j] = H0_jp1j;

        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.w[k]/H0_jp1j;
                blk.V[k,j+1] = blk.v[k];
            }
        }

        // Build rotated Hessenberg progressively
        if (j != 0) {
            // Extract final column in H
            foreach (i; 0 .. j+1) h[i] = H0[i,j];
            // Rotate column by previous rotations (stored in Q0)
            nm.bbla.dot(Q0, j+1, j+1, h, hR);
            // Place column back in H
            foreach (i; 0 .. j+1) H0[i,j] = hR[i];
        }
        // Now form new Gamma
        Gamma.eye();
        auto denom = sqrt(H0[j,j]*H0[j,j] + H0[j+1,j]*H0[j+1,j]);
        auto s_j = H0[j+1,j]/denom;
        auto c_j = H0[j,j]/denom;
        Gamma[j,j] = c_j; Gamma[j,j+1] = s_j;
        Gamma[j+1,j] = -s_j; Gamma[j+1,j+1] = c_j;
        // Apply rotations
        nm.bbla.dot(Gamma, j+2, j+2, H0, j+1, H1);
        nm.bbla.dot(Gamma, j+2, j+2, g0, g1);
        // Accumulate Gamma rotations in Q.
        if (j == 0) {
            copy(Gamma, Q1);
        }
        else {
            nm.bbla.dot!number(Gamma, j+2, j+2, Q0, j+2, Q1);
        }

        // Prepare for next step
        copy(H1, H0);
        g0[] = g1[];
        copy(Q1, Q0);

        // Get residual
        auto resid = fabs(g1[j+1]);
        auto linSolResid = (resid/beta0).re;
        if (resid <= targetResidual) {
            isConverged = true;
            break;
        }
    }

    /*
    debug {
        writeln("DEBUG: performIterations()");
        writeln(" done Krylov iterations.");
    }
    */
    
    
    return isConverged;
}


/**
 * Apply preconditioning to GMRES iterations.
 *
 * NOTE: Although you don't see "dt" in the code here, it is passed in because
 *       it is accessed in the mixin code.
 * Authors: KAD and RJG
 * Date: 2022-07-09
 */
void applyPreconditioning()
{
    auto nConserved = GlobalConfig.cqi.n;

    final switch (nkCfg.preconditioner) {
    case PreconditionerType.lusgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = to!number(0.0); }
        mixin(lusgs_solve("zed", "v"));
        break;
    case PreconditionerType.diagonal:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = to!number(0.0); }
        mixin(diagonal_solve("zed", "v"));
        break;
    case PreconditionerType.jacobi:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = to!number(0.0); }
        mixin(jacobi_solve("zed", "v"));
        break;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = to!number(0.0); }
        mixin(sgs_solve("zed", "v"));
        break;
    case PreconditionerType.ilu:
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.zed[] = blk.v[];
            nm.smla.solve(blk.flowJacobian.local, blk.zed);
        }
        break;
    } // end switch
}


/**
 * Remove preconditioning on values and place in dU.
 *
 * NOTE: Although you don't see "dt" in the code here, it is passed in because
 *       it is accessed in the mixin code.
 *
 * Authors: KAD and RJG
 * Date: 2022-07-09
 */

void removePreconditioning()
{
    auto nConserved = GlobalConfig.cqi.n;
    
    final switch (nkCfg.preconditioner) {
    case PreconditionerType.lusgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = to!number(0.0); }
        mixin(lusgs_solve("dU", "zed"));
        break;
    case PreconditionerType.diagonal:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = to!number(0.0); }
        mixin(diagonal_solve("dU", "zed"));
        break;
    case PreconditionerType.jacobi:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = to!number(0.0); }
        mixin(jacobi_solve("dU", "zed"));
        break;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = to!number(0.0); }
        mixin(sgs_solve("dU", "zed"));
        break;
    case PreconditionerType.ilu:
        foreach(blk; parallel(localFluidBlocks,1)) {
            blk.dU[] = blk.zed[];
            nm.smla.solve(blk.flowJacobian.local, blk.dU);
        }
        break;
    } // end switch
}


/**
 * Compute perturbation size estimate for real-valued Frechet derivative.
 *
 * REFERENCE: [ASK KYLE]
 *
 * Authors: KAD and RJG
 * Date: 2022-03-02
 */
number computePerturbationSize()
{
    // calculate sigma without scaling
    number sumv = 0.0;
    mixin(dotOverBlocks("sumv", "zed", "zed"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(sumv.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    size_t nConserved = GlobalConfig.cqi.n;
    auto eps0 = 1.0e-6; // FIX ME: Check with Kyle about what to call this
                        //         and how to set the value.
    number N = 0.0;
    number sume = 0.0;
    foreach (blk; parallel(localFluidBlocks,1)) {
        int cellCount = 0;
        foreach (cell; blk.cells) {
            foreach (val; cell.U[0]) {
                sume += eps0*abs(val) + eps0;
                N += 1;
            }
            cellCount += nConserved;
        }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(sume.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(N.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    number sigma = (sume/(N*sqrt(sumv))).re;
    return sigma;
}

/**
 * Evaluate J*v via a Frechet derivative.
 *
 * This is just a delegator function through to a complex-valued perturbation evaulation
 * or a real-valued perturbation evalulation.
 * Authors: KAD and RJG
 * Date: 2022-03-02
 */
void evalJacobianVectorProduct(double sigma)
{
    version (complex_numbers) {
        evalComplexMatVecProd(sigma);
    }
    else {
        evalRealMatVecProd(sigma);
    }
}



void evalResidual(int ftl)
{
    fnCount++;
    int gtl = 0;
    double dummySimTime = -1.0;

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
    }
    exchange_ghost_cell_boundary_data(dummySimTime, 0, ftl);
    foreach (blk; localFluidBlocks) {
        blk.applyPreReconAction(dummySimTime, 0, ftl);
    }

    // We don't want to switch between flux calculator application while
    // doing the Frechet derivative, so we'll only search for shock points
    // at ftl = 0, which is when the F(U) evaluation is made.
    if (ftl == 0 && GlobalConfig.do_shock_detect) { detect_shocks(0, ftl); }

    // We need to apply the copy_cell_data BIE at this point to allow propagation of
    // "shocked" cell information (fs.S) to the boundary interface BEFORE the convective
    // fluxes are evaluated. This is important for a real-valued Frechet derivative
    // with adaptive fluxes, to ensure that each interface along the boundary uses a consistent flux
    // calculator for both the baseline residual R(U) and perturbed residual R(U+dU) evaluations.
    // [TODO] KD 2021-11-30 This is a temporary fix until a more formal solution has been decided upon.
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach(boundary; blk.bc) {
            if (boundary.preSpatialDerivActionAtBndryFaces[0].desc == "CopyCellData") {
                boundary.preSpatialDerivActionAtBndryFaces[0].apply(dummySimTime, gtl, ftl);
            }
        }
    }

    bool allow_high_order_interpolation = true;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase0(allow_high_order_interpolation, gtl);
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(dummySimTime, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase1(allow_high_order_interpolation, gtl);
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(dummySimTime, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase2(allow_high_order_interpolation, gtl);
    }
    
    foreach (blk; localFluidBlocks) {
        blk.applyPostConvFluxAction(dummySimTime, gtl, ftl);
    }

    if (GlobalConfig.viscous) {
        foreach (blk; localFluidBlocks) {
            blk.applyPreSpatialDerivActionAtBndryFaces(dummySimTime, gtl, ftl);
            blk.applyPreSpatialDerivActionAtBndryCells(dummySimTime, gtl, ftl);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.flow_property_spatial_derivatives(0);
        }
        // for unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
        // we need to transfer the viscous gradients before the flux calc
        exchange_ghost_cell_boundary_viscous_gradient_data(dummySimTime, gtl, ftl);
        foreach (blk; parallel(localFluidBlocks,1)) {
            // we need to average cell-centered spatial (/viscous) gradients to get approximations of the gradients
            // at the cell interfaces before the viscous flux calculation.
            if (blk.myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                foreach(f; blk.faces) {
                    f.average_cell_deriv_values(0);
                }
            }
            blk.estimate_turbulence_viscosity();
        }
        // we exchange boundary data at this point to ensure the
        // ghost cells along block-block boundaries have the most
        // recent mu_t and k_t values.
        exchange_ghost_cell_turbulent_viscosity();
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.viscous_flux();
        }
        foreach (blk; localFluidBlocks) {
            blk.applyPostDiffFluxAction(dummySimTime, 0, ftl);
        }
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        // the limit_factor is used to slowly increase the magnitude of the
        // thermochemical source terms from 0 to 1 for problematic reacting flows
        double limit_factor = 1.0;
        if (blk.myConfig.nsteps_of_chemistry_ramp > 0) {
            double S = SimState.step/to!double(blk.myConfig.nsteps_of_chemistry_ramp);
            limit_factor = min(1.0, S);
        }
        foreach (i, cell; blk.cells) {
            cell.add_inviscid_source_vector(0, 0.0);
            if (blk.myConfig.viscous) {
                cell.add_viscous_source_vector();
            }
            if (blk.myConfig.reacting) {
                cell.add_thermochemical_source_vector(blk.thermochem_source, limit_factor);
            }
            if (blk.myConfig.udf_source_terms) {
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                if (blk.grid_type == Grid_t.structured_grid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                    i_cell = ijk_indices[0];
                    j_cell = ijk_indices[1];
                    k_cell = ijk_indices[2];
                }
                getUDFSourceTermsForCell(blk.myL, cell, 0, dummySimTime, blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                cell.add_udf_source_vector();
            }
            cell.time_derivatives(0, ftl);
        }
    }
}


/**
 * Evaluate J*v via a Frechet derivative with perturbation in imaginary plane.
 *
 * J*v = Im( R(U + sigma*j)/sigma )
 *
 * Authors: KAD and RJG
 * Date: 2022-03-03
 */
void evalComplexMatVecProd(double sigma)
{
    version(complex_numbers) {
        alias cfg = GlobalConfig;
        
        foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(activePhase.jacobianInterpolationOrder); }
        if (activePhase.frozenLimiterForJacobian) {
            foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = true; }
        }
        // Make a stack-local copy of conserved quantities info
        size_t nConserved = GlobalConfig.cqi.n;

        //writeln("in evalComplexMatVecProd");
        //writefln("sigma= %e", sigma);
        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.clear_fluxes_of_conserved_quantities();
            foreach (i, cell; blk.cells) cell.clear_source_vector();
            size_t startIdx = 0;
            foreach (i, cell; blk.cells) {
                cell.U[1][] = cell.U[0][];
                foreach (ivar; 0 .. nConserved) {
                    cell.U[1][ivar] += complex(0.0, sigma * blk.zed[startIdx+ivar].re);
                }
                cell.decode_conserved(0, 1, 0.0);
                startIdx += nConserved;
            }
        }
        evalResidual(1);
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t startIdx = 0;
            foreach (cell; blk.cells) {
                foreach (ivar; 0 .. nConserved) {
                    blk.zed[startIdx+ivar] = cell.dUdt[1][ivar].im/(sigma);
                }
                startIdx += nConserved;
            }
            // we must explicitly remove the imaginary components from the cell and interface flowstates
            foreach(cell; blk.cells) { cell.fs.clear_imaginary_components(); }
            foreach(bc; blk.bc) {
                foreach(ghostcell; bc.ghostcells) { ghostcell.fs.clear_imaginary_components(); }
            }
            foreach(face; blk.faces) { face.fs.clear_imaginary_components(); }
        }
        foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(activePhase.residualInterpolationOrder); }
    } else {
        throw new Error("Oops. Steady-State Solver setting: useComplexMatVecEval is not compatible with real-number version of the code.");
    }
}

/**
 * Evaluate J*v via a Frechet derivative
 *
 * J*v = (R(U + sigma) - R(U)) / sigma
 *
 * Authors: KAD and RJG
 * Date: 2022-03-03
 */
void evalRealMatVecProd(double sigma)
{
    foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(activePhase.jacobianInterpolationOrder); }
    if (activePhase.frozenLimiterForJacobian) {
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = true; }
    }
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;

    // We perform a Frechet derivative to evaluate J*D^(-1)v
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
        size_t startIdx;
        foreach (cell; blk.cells) {
            cell.U[1][] = cell.U[0][];
            foreach (ivar; 0 .. nConserved) {
                cell.U[1][ivar] += sigma*blk.zed[startIdx+ivar];
            }
            cell.decode_conserved(0, 1, 0.0);
            startIdx += nConserved;
        }
    }
    evalResidual(1);
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t startIdx = 0;
        foreach (cell; blk.cells) {
            foreach (ivar; 0 .. nConserved) {
                blk.zed[startIdx+ivar] = (cell.dUdt[1][ivar] - blk.R[startIdx+ivar])/(sigma);
            }
            cell.decode_conserved(0, 0, 0.0);
            startIdx += nConserved;
        }
    }
    foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(activePhase.residualInterpolationOrder); }
}

/**
 * Determine a relaxation factor.
 *
 * In this algorithm, the relaxation factor keeps falling as we search across cells in order.
 * This is efficient since we're searching for a worst case. If at any point, the relaxation
 * factor gets smaller than what we're prepared to accept then we just break the search
 * in that block of cells.
 *
 * Authors: KAD and RJG
 * Date: 2022-03-05
 */
double determineRelaxationFactor()
{
    alias GlobalConfig cfg;
    double theta = nkCfg.allowableRelativeMassChange;
    double minOmega = nkCfg.minRelaxationFactor;
    double omegaReductionFactor = nkCfg.relaxationFactorReductionFactor;
    
    size_t nConserved = cfg.cqi.n;
    size_t massIdx = cfg.cqi.mass;
    
    double omega = 1.0;

    //----
    // 1. First determine a relaxation factor based on an allowable amount of mass change
    //----
    
    foreach (blk; parallel(localFluidBlocks,1)) {
        auto cqi = blk.myConfig.cqi;
        int startIdx = 0;
        blk.omegaLocal = 1.0;
        number U, dU, relDiffLimit;
        foreach (cell; blk.cells) {
            if (cqi.n_species == 1) {
                U = cell.U[0][cqi.mass];
                dU = blk.dU[startIdx+cqi.mass];
            }
            else {
                U = 0.0;
                dU = 0.0;
                foreach (isp; 0 .. cqi.n_species) {
                    U += cell.U[0][cqi.species+isp];
                    dU += blk.dU[startIdx+cqi.species+isp];
                }
            }
            relDiffLimit = fabs(dU/(theta*U));
            blk.omegaLocal = 1.0/(fmax(relDiffLimit.re, 1.0/blk.omegaLocal));
            startIdx += nConserved;
        }
    }
    // In serial, find minimum omega across all blocks.
    foreach (blk; localFluidBlocks) omega = fmin(omega, blk.omegaLocal);
    version (mpi_parallel) {
        // In parallel, find minimum and communicate to all processes
        MPI_Allreduce(MPI_IN_PLACE, &(omega), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }

    //----
    // 2. Now check if primitives are valid. If not adjust relaxation factor until they are.
    //----

    foreach (blk; parallel(localFluidBlocks, 1)) {
        int startIdx = 0;
        foreach (cell; blk.cells) {
            bool failedDecode = false;
            while (blk.omegaLocal >= minOmega) {
                cell.U[1][] = cell.U[0][];
                foreach (i; 0 .. nConserved) {
                    cell.U[1][i] = cell.U[0][i] + blk.omegaLocal * blk.dU[startIdx+i];
                }
                try {
                    cell.decode_conserved(0, 1, 0.0);
                }
                catch (FlowSolverException e) {
                    failedDecode = true;
                }
                // return cell to original state
                cell.decode_conserved(0, 0, 0.0);

                if (failedDecode) {
                    blk.omegaLocal *= omegaReductionFactor;
                    failedDecode = false;
                }
                else {
                    // Update was good, so omega is ok for this cell.
                    break;
                }
            }
            startIdx += nConserved;
        }
    }
    // In serial, find minimum omega across all blocks.
    foreach (blk; localFluidBlocks) omega = fmin(omega, blk.omegaLocal);
    version (mpi_parallel) {
        // In parallel, find minimum and communicate to all processes
        MPI_Allreduce(MPI_IN_PLACE, &(omega), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    
    return omega;
}


/**
 * Apply Newton update, scaled by relaxation factor, omega.
 *
 * Authors: KAD and RJG
 * Date: 2022-03-05
 */
void applyNewtonUpdate(double relaxationFactor)
{
    size_t nConserved = GlobalConfig.cqi.n;
    double omega = relaxationFactor;
    
    foreach (blk; parallel(localFluidBlocks, 1)) {
        size_t startIdx = 0;
        foreach (cell; blk.cells) {
            foreach (ivar; 0 .. nConserved) {
                cell.U[0][ivar] += omega * blk.dU[startIdx + ivar];
            }
            try {
                cell.decode_conserved(0, 0, 0.0);
            }
            catch (FlowSolverException e) {
                // If the physicality check + relaxation factor is in play,
                // then we should never fail here since we've pre-tested
                // the update before applying.
                // However, we may not be using that check.
                // In which case, it's probably best to bail out here.
                string errMsg;
                errMsg ~= "Failed update for cell when applying Newton step.\n";
                errMsg ~= format("blk-id: %d", blk.id);
                errMsg ~= format("cell-id: %d", cell.id);
                errMsg ~= format("cell-pos: %s", cell.pos[0]);
                errMsg ~= "Error message from decode_conserved:\n";
                errMsg ~= e.msg;
                throw new NewtonKrylovException(errMsg);
            }
            startIdx += nConserved;
        }
    }
}

/**
 * Initialise diagnostics file.
 *
 * Authors: KAD and RJG
 * Date: 2022-07-24
 */
void initialiseDiagnosticsFile()
{
    ensure_directory_is_present(diagnosticsDir);
    auto fname = diagnosticsDir ~ "/" ~ diagnosticsFilename;
    diagnostics = File(fname, "w");
    
    alias cfg = GlobalConfig;

    diagnostics.writeln("#  1: step");
    diagnostics.writeln("#  2: dt");
    diagnostics.writeln("#  3: CFL");
    diagnostics.writeln("#  4: linear-solve-residual-target");
    diagnostics.writeln("#  5: nRestarts");
    diagnostics.writeln("#  6: nIters");
    diagnostics.writeln("#  7: nFnCalls");
    diagnostics.writeln("#  8: wall-clock, s");
    diagnostics.writeln("#  9: global-residual-abs");
    diagnostics.writeln("# 10: global-residual-rel");
    diagnostics.writeln("# 11: mass-balance");
    diagnostics.writeln("# 12: linear-solve-residual");
    diagnostics.writeln("# 13: omega");
    int n_entry = 14;
    foreach (ivar; 0 .. cfg.cqi.n) {
        diagnostics.writefln("# %02d: %s-abs", n_entry, cfg.cqi.names[ivar]);
        n_entry++;
        diagnostics.writefln("# %02d: %s-rel", n_entry, cfg.cqi.names[ivar]);
        n_entry++;
    }
    diagnostics.close();
}
        

/**
 * Update diagnostics file with current status.
 *
 * Authors: KAD and RJG
 * Date: 2022-07-24
 */
void writeDiagnostics(int step, double dt, double cfl, double wallClockElapsed, double omega, ref bool residualsUpToDate)
{
    alias cfg = GlobalConfig;

    // [TODO] compute mass balance
    // double massBalance = computeMassBalance();
    double massBalance = 100.0;
    if (!residualsUpToDate) {
        computeResiduals(currentResiduals);
        residualsUpToDate = true;
    }

    // We don't need to proceed on ranks other than master.
    if (!cfg.is_master_task) return;
    
    auto fname = diagnosticsDir ~ "/" ~ diagnosticsFilename;
    diagnostics = File(fname, "a");
    diagnostics.writef("%8d %20.16e %20.16e %20.16e %.8f %2d %3d %8d %20.16e %20.16e %20.16e %20.16e ",
                       step, dt, cfl, activePhase.linearSolveTolerance, wallClockElapsed, 
                       gmresInfo.nRestarts, gmresInfo.iterationCount, fnCount,
                       globalResidual, globalResidual/referenceGlobalResidual,
                       massBalance, gmresInfo.finalResidual, omega);
    foreach (ivar; 0 .. cfg.cqi.n) {
        diagnostics.writef("%20.16e %20.16e ", currentResiduals[ivar].re, currentResiduals[ivar].re/referenceResiduals[ivar].re);
    }
    diagnostics.writef("\n");
    diagnostics.close();
}

void writeSnapshot(int step, ref int nWrittenSnapshots)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writeln();
        writefln("+++++++++++++++++++++++++++++++++++++++");
        writefln("+   Writing snapshot at step = %4d   +", step);
        writefln("+++++++++++++++++++++++++++++++++++++++");
        writeln();
    }

    nWrittenSnapshots++;

    if (nWrittenSnapshots <= nkCfg.totalSnapshots) {
        auto dirName = steadyFlowDirectory(nWrittenSnapshots);        
        if (cfg.is_master_task) {
            ensure_directory_is_present(dirName);
        }

        foreach (blk; localFluidBlocks) {
            foreach (io; blk.block_io) {
                auto fileName = steadyFlowFilename(nWrittenSnapshots, blk.id);
                if (io.do_save()) io.save_to_file(fileName, dummySimTime);
            }
        }
    }
    else {
        // We need to shuffle all of the snapshots
        foreach (iSnap; 2 .. nkCfg.totalSnapshots+1) {
            foreach (blk; localFluidBlocks) {
                foreach (io; blk.block_io) {
                    auto fromName = steadyFlowFilename(iSnap, blk.id);
                    auto toName = steadyFlowFilename(iSnap-1, blk.id);
                    rename(fromName, toName);
                }
            }
        }
        foreach (blk; localFluidBlocks) {
            foreach (io; blk.block_io) {
                auto fileName = steadyFlowFilename(nkCfg.totalSnapshots, blk.id);
                if (io.do_save()) io.save_to_file(fileName, dummySimTime);
            }
        }
    }
  
}


void printStatusToScreen(int step, double cfl, double dt, double wallClockElapsed, ref bool residualsUpToDate)
{
    alias cfg = GlobalConfig;
    
    if (!residualsUpToDate) {
        computeResiduals(currentResiduals);
        residualsUpToDate = true;
    }

    // We don't need to proceed on ranks other than master.
    if (!cfg.is_master_task) return;

    auto cqi = GlobalConfig.cqi;
    auto writer = appender!string();
    string hrule = "--------------------------------------------------------------\n";
    formattedWrite(writer, hrule);
    formattedWrite(writer, "step= %6d\tcfl=%10.3e\tdt=%10.3e\tWC=%.1f\n",
                   step, cfl, dt, wallClockElapsed);
    formattedWrite(writer, hrule);
    formattedWrite(writer, "%12s\t\tRELATIVE\t\tABSOLUTE\n", "RESIDUALS");
    formattedWrite(writer, "%12s\t\t%10.6e\t\t%10.6e\n", "global", globalResidual.re/referenceGlobalResidual.re, globalResidual.re);
    foreach (ivar; 0 .. cqi.n) {
        formattedWrite(writer, "%12s\t\t%10.6e\t\t%10.6e\n", cqi.names[ivar], currentResiduals[ivar].re/referenceResiduals[ivar].re, currentResiduals[ivar].re);
    }        
    writeln(writer.data);
}


/*---------------------------------------------------------------------
 * Mixins for performing preconditioner actions
 *---------------------------------------------------------------------
 */

string diagonal_solve(string lhs_vec, string rhs_vec)
{
    string code = "{

    foreach (blk; parallel(localFluidBlocks,1)) {
        nm.smla.multiply(blk.flowJacobian.local, blk."~rhs_vec~", blk."~lhs_vec~"[]);
    }

    }";
    return code;
}

string jacobi_solve(string lhs_vec, string rhs_vec)
{
    string code = "{

    int kmax = nkCfg.preconditionerSubIterations;
    foreach (k; 0 .. kmax) {
         foreach (blk; parallel(localFluidBlocks,1)) {
               blk.rhs[] = to!number(0.0);
               nm.smla.multiply_block_upper_triangular(blk.flowJacobian.local, blk."~lhs_vec~"[], blk.rhs[], blk.cells.length, nConserved);
               nm.smla.multiply_block_lower_triangular(blk.flowJacobian.local, blk."~lhs_vec~"[], blk.rhs[], blk.cells.length, nConserved);
               blk.rhs[] = blk."~rhs_vec~"[] - blk.rhs[];
               blk."~lhs_vec~"[] = to!number(0.0);
               nm.smla.multiply_block_diagonal(blk.flowJacobian.local, blk.rhs[], blk."~lhs_vec~"[], blk.cells.length, nConserved);
         }
    }

    }";
    return code;
}

string sgs_solve(string lhs_vec, string rhs_vec)
{
    string code = "{

    int kmax = nkCfg.preconditionerSubIterations;
    foreach (k; 0 .. kmax) {
         foreach (blk; parallel(localFluidBlocks,1)) {
             // forward sweep
             blk.rhs[] = to!number(0.0);
             nm.smla.multiply_block_upper_triangular(blk.flowJacobian.local, blk."~lhs_vec~", blk.rhs, blk.cells.length, nConserved);
             blk.rhs[] = blk."~rhs_vec~"[] - blk.rhs[];
             nm.smla.block_lower_triangular_solve(blk.flowJacobian.local, blk.rhs[], blk."~lhs_vec~"[], blk.cells.length, nConserved);

             // backward sweep
             blk.rhs[] = to!number(0.0);
             nm.smla.multiply_block_lower_triangular(blk.flowJacobian.local, blk."~lhs_vec~"[], blk.rhs[], blk.cells.length, nConserved);
             blk.rhs[] = blk."~rhs_vec~"[] - blk.rhs[];
             nm.smla.block_upper_triangular_solve(blk.flowJacobian.local, blk.rhs, blk."~lhs_vec~"[], blk.cells.length, nConserved);
         }
    }

    }";
    return code;
}

string lusgs_solve(string lhs_vec, string rhs_vec)
{
    string code = "{

    int kmax = nkCfg.preconditionerSubIterations;
    bool matrix_based = false;
    double omega = 2.0;
    double dummySimTime = -1.0;

    // 1. initial subiteration
    foreach (blk; parallel(localFluidBlocks,1)) {
        int startIdx = 0;
        foreach (cell; blk.cells) {
            number dtInv = 1.0/cell.dt_local;
            auto z_local = blk."~lhs_vec~"[startIdx .. startIdx+nConserved]; // this is actually a reference not a copy
            cell.lusgs_startup_iteration(dtInv, omega, z_local, blk."~rhs_vec~"[startIdx .. startIdx+nConserved]);
            startIdx += nConserved;
        }
    }

    // 2. kmax subiterations
    foreach (k; 0 .. kmax) {
         // shuffle dU values
         foreach (blk; parallel(localFluidBlocks,1)) {
             int startIdx = 0;
             foreach (cell; blk.cells) {
                 cell.dUk[0 .. nConserved] = blk."~lhs_vec~"[startIdx .. startIdx+nConserved];
                 startIdx += nConserved;
             }
         }

         // exchange boundary dU values
         exchange_ghost_cell_boundary_data(dummySimTime, 0, 0);

         // perform subiteraion
         foreach (blk; parallel(localFluidBlocks,1)) {
             int startIdx = 0;
             foreach (cell; blk.cells) {
                  auto z_local = blk."~lhs_vec~"[startIdx .. startIdx+nConserved]; // this is actually a reference not a copy
                  cell.lusgs_relaxation_iteration(omega, matrix_based, z_local, blk."~rhs_vec~"[startIdx .. startIdx+nConserved]);
                  startIdx += nConserved;
             }
         }
    }

    }";
    return code;
}
