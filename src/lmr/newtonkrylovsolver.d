/**
 * Core set of functions used in the Newton-Krylov updates for steady-state convergence.
 *
 * Authors: RJG and KAD
 * Date: 2022-02-28
 * History:
 *   2022-02-28 Fairly aggressive refactor of code that was in: steadystate_core.d
 *   2022-10-22 Reworked as part of lmr5
 */

module newtonkrylovsolver;

import core.stdc.stdlib : exit;
import core.memory : GC;
import std.algorithm : min;
import std.algorithm.searching : countUntil;
import std.datetime : DateTime, Clock;
import std.parallelism : parallel, defaultPoolThreads;
import std.stdio : File, writeln, writefln, stdout;
import std.file : rename, readText;
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

import lmrexceptions;
import lmrconfig;
import conservedquantities : ConservedQuantities, copy_values_from;
import fileutil : ensure_directory_is_present;

import globalconfig;
import globaldata;
import init;
import simcore_gasdynamic_step : detect_shocks;
import simcore_exchange;
import bc;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;
import user_defined_source_terms : getUDFSourceTermsForCell;
import blockio;
import fvcellio;

version(mpi_parallel) {
    import mpi;
}

/*---------------------------------------------------------------------
 * Module globals
 * Typically for diagnostics and/or debugging.
 *---------------------------------------------------------------------
 */

File diagnostics;
string diagDir;
string diagFile;
static this()
{
    diagDir = lmrCfg.simDir ~ "/" ~ lmrJSONCfg["diagnostics-directory"].str;
    diagFile = diagDir ~ "/" ~ lmrJSONCfg["nk-diagnostics-filename"].str;
}

FVCellIO limCIO;
BlockIO limBlkIO;

/*---------------------------------------------------------------------
 * Enums for preconditioners
 *---------------------------------------------------------------------
 */
enum PreconditionerType { diagonal, jacobi, sgs, ilu }

string preconditionerTypeName(PreconditionerType i)
{
    final switch (i) {
    case PreconditionerType.diagonal: return "diagonal";
    case PreconditionerType.jacobi: return "jacobi";
    case PreconditionerType.sgs: return "sgs";
    case PreconditionerType.ilu: return "ilu";
    }
}

PreconditionerType preconditionerTypeFromName(string name)
{
    switch (name) {
    case "diagonal": return PreconditionerType.diagonal;
    case "jacobi": return PreconditionerType.jacobi;
    case "sgs": return PreconditionerType.sgs;
    case "ilu": return PreconditionerType.ilu;
    default:
        string errMsg = "The selected 'preconditioner' is unavailable.\n";
        errMsg ~= format("You selected: '%s'\n", name);
        errMsg ~= "The available strategies are: \n";
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
    int numberOfStepsForSettingReferenceResiduals = 0;
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
    bool inviscidCFLOnly = true;
    bool useLineSearch = true;
    bool usePhysicalityCheck = true;
    double allowableRelativeMassChange = 0.2;
    double minRelaxationFactor = 0.1;
    double relaxationFactorReductionFactor = 0.7;
    bool useResidualSmoothing = false;
    // Linear solver and preconditioning
    int maxLinearSolverIterations = 10;
    int maxLinearSolverRestarts = 0;
    bool useScaling = true;
    bool useRealValuedFrechetDerivative = false;
    double frechetDerivativePerturbation = 1.0e-30;
    bool usePreconditioner = true;
    double preconditionerPerturbation = 1.0e-30;
    PreconditionerType preconditioner = PreconditionerType.ilu;
    // ILU setting
    int iluFill = 0;
    // sub iterations for preconditioners
    int preconditionerSubIterations = 4;
    // output and diagnostics
    int stepsBetweenStatus = 10;
    int totalSnapshots = 5;
    int stepsBetweenSnapshots = 10;
    int stepsBetweenDiagnostics = 10;
    int stepsBetweenLoadsUpdate = 20;
    bool writeSnapshotOnLastStep = true;
    bool writeDiagnosticsOnLastStep = true;
    bool writeLimiterValues = false;

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
        inviscidCFLOnly = getJSONbool(jsonData, "inviscid_cfl_only", inviscidCFLOnly);
        useLineSearch = getJSONbool(jsonData, "use_line_search", useLineSearch);
        usePhysicalityCheck = getJSONbool(jsonData, "use_physicality_check", usePhysicalityCheck);
        allowableRelativeMassChange = getJSONdouble(jsonData, "allowable_relative_mass_change", allowableRelativeMassChange);
        minRelaxationFactor = getJSONdouble(jsonData, "min_relaxation_factor", minRelaxationFactor);
        relaxationFactorReductionFactor = getJSONdouble(jsonData, "relaxation_factor_reduction_factor", relaxationFactorReductionFactor);
	useResidualSmoothing = getJSONbool(jsonData, "use_residual_smoothing", useResidualSmoothing);
	maxLinearSolverIterations = getJSONint(jsonData, "max_linear_solver_iterations", maxLinearSolverIterations);
        maxLinearSolverRestarts = getJSONint(jsonData, "max_linear_solver_restarts", maxLinearSolverRestarts);
        useScaling = getJSONbool(jsonData, "use_scaling", useScaling);
        useRealValuedFrechetDerivative = getJSONbool(jsonData, "use_real_valued_frechet_derivative", useRealValuedFrechetDerivative);
        frechetDerivativePerturbation = getJSONdouble(jsonData, "frechet_derivative_perturbation", frechetDerivativePerturbation);
        usePreconditioner = getJSONbool(jsonData, "use_preconditioner", usePreconditioner);
        preconditionerPerturbation = getJSONdouble(jsonData, "preconditioner_perturbation", preconditionerPerturbation);
        auto pString = getJSONstring(jsonData, "preconditioner", "NO_SELECTION_SUPPLIED");
        preconditioner = preconditionerTypeFromName(pString);
        iluFill = getJSONint(jsonData, "ilu_fill", iluFill);
        preconditionerSubIterations = getJSONint(jsonData, "preconditioner_sub_iterations", preconditionerSubIterations);
        stepsBetweenStatus = getJSONint(jsonData, "steps_between_status", stepsBetweenStatus);
        totalSnapshots = getJSONint(jsonData, "total_snapshots", totalSnapshots);
        stepsBetweenSnapshots = getJSONint(jsonData, "steps_between_snapshots", stepsBetweenSnapshots);
        stepsBetweenDiagnostics = getJSONint(jsonData, "steps_between_diagnostics", stepsBetweenDiagnostics);
        stepsBetweenLoadsUpdate = getJSONint(jsonData, "steps_between_loads_update", stepsBetweenLoadsUpdate);
        writeSnapshotOnLastStep = getJSONbool(jsonData, "write_snapshot_on_last_step", writeSnapshotOnLastStep);
        writeDiagnosticsOnLastStep = getJSONbool(jsonData, "write_diagnostics_on_last_step", writeDiagnosticsOnLastStep);
        writeLimiterValues = getJSONbool(jsonData, "write_limiter_values", writeLimiterValues);
    }
}
NKGlobalConfig nkCfg;

struct NKPhaseConfig {
    bool useLocalTimestep = true;
    int residualInterpolationOrder = 2;
    int jacobianInterpolationOrder = 2;
    bool frozenPreconditioner = true;
    int stepsBetweenPreconditionerUpdate = 10;
    bool useAdaptivePreconditioner = false;
    bool ignoreStoppingCriteria = true;
    bool frozenLimiterForJacobian = false;
    double linearSolveTolerance = 0.01;
    // Auto CFL control
    bool useAutoCFL = false;
    double thresholdRelativeResidualForCFLGrowth = 0.99;
    double startCFL = 1.0;
    double maxCFL = 1000.0;
    double autoCFLExponent = 0.75;

    void readValuesFromJSON(JSONValue jsonData)
    {
        useLocalTimestep = getJSONbool(jsonData, "use_local_timestep", useLocalTimestep);
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
    double dt;
    double cfl;
    int step;
    double globalResidual;
    double prevGlobalResidual;
    ConservedQuantities residuals;

    this(size_t n)
    {
        residuals = new ConservedQuantities(n);
    }

    this(this)
    {
	residuals = residuals.dup;
    }

}
RestartInfo[] snapshots;

// For MPI, it's more robust on large clusters to broadcast
// residuals and restart info than to rely on reading from file.
// This was reported by Nick Gibbons. It seems the fact that these
// files are very small might be the issue. Reading from them
// from many processes on a networked filesystem does not seem
// to work reliably.

version(mpi_parallel){

/**
 * Helper function for sending around a conserved quantities object.
 *
 * Authors: NNG and RJG
 * Date: 2023-04-11
 */
void broadcastResiduals(ref ConservedQuantities residuals)
{
    double[] buffer;

    int size = to!int(residuals.length);
    version(complex_numbers) { size *= 2; }
    buffer.length = size;

    size_t i = 0;
    foreach (f; residuals){
        buffer[i] = f.re;
        i++;
        version(complex_numbers){
            buffer[i] = f.im;
            i++;
        }
    }

    MPI_Bcast(buffer.ptr, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    i = 0;
    foreach (j; 0 .. residuals.length){
        residuals[j].re = buffer[i];
        i++;
        version(complex_numbers) {
            residuals[j].im = buffer[i];
            i++;
        }
    }
    return;
}

/**
 * Helper function for sending around a collection of RestartInfo objects
 *
 * Notes: For some reason, this function does need a ref in its argument,
 * even though dynamic arrays are supposed to be reference types...
 *
 * Authors: NNG and RJG
 * Date: 2023-04-11
 */
void broadcastRestartInfo()
{
    // The main issue here is that the non-master processes do not actually know
    // how many of these things there are. First we need sort that out.
    // FIXME: Note that this could cause problems if this function is reused.
    int n_snapshots = to!int(snapshots.length);
    MPI_Bcast(&n_snapshots, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!GlobalConfig.is_master_task) {
        foreach (i; 0 .. n_snapshots){
            snapshots ~= RestartInfo(GlobalConfig.cqi.n);
        }
    }

    // Now that we have them allocated, we can fill them out
    foreach (i; 0 .. n_snapshots){
        MPI_Bcast(&(snapshots[i].dt), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(snapshots[i].cfl), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(snapshots[i].step), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(snapshots[i].globalResidual), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(snapshots[i].prevGlobalResidual), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        broadcastResiduals(snapshots[i].residuals);
    }
}

} // end version(mpi_parallel)



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
        writeln("lmr run-steady: Begin initNewtonKrylovSimulation()...");
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

    initThreadPool(maxCPUs, threadsPerMPITask);

    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowFieldSteadyMode(snapshotStart);

    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

    initFullFaceDataExchange();
    initMappedCellDataExchange();
    initGhostCellGeometry();
    initLeastSquaresStencils();

    if ((cfg.interpolation_order > 1) &&
	((cfg.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp) ||
	 (cfg.unstructured_limiter == UnstructuredLimiter.venkat_mlp))) {
        initMLPlimiter();
    }

    orderBlocksBySize();
    initMasterLuaState();
    initCornerCoordinates();
    if (cfg.turb_model.needs_dwall) initWallDistances();

    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

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
        writefln("lmr run-steady: Done initNewtonKrylovSimulation() at wall-clock(WC)= %.1f sec", wall_clock_elapsed);
        stdout.flush();
    }
}

void readNewtonKrylovConfig()
{
    alias cfg = GlobalConfig;

    if (cfg.verbosity_level > 1) writeln("Read N-K config file.");
    JSONValue jsonData = readJSONfile(lmrCfg.nkCfgFile);
    nkCfg.readValuesFromJSON(jsonData);
    // Allocate space and configure phases
    nkPhases.length = nkCfg.numberOfPhases;
    foreach (i, ref phase; nkPhases) {
        string key = "NewtonKrylovPhase_" ~ to!string(i);
        phase.readValuesFromJSON(jsonData[key]);
    }

    // Perform some consistency checks
    if (nkCfg.phaseChangesAtSteps.length < nkCfg.numberOfPhases-1) {
        string errMsg;
        errMsg ~= format("ERROR: number of phases = %d but number of entries in phase changes list is: %d\n",
                nkCfg.numberOfPhases, nkCfg.phaseChangesAtSteps.length);
        errMsg ~= "       We expect at least (number of phases) - 1 entries in phase changes list.\n";
        errMsg ~= "       These entries are the step number at which to change from one phase to the next.\n";
        throw new Error(errMsg);
    }
    foreach (i, phase; nkPhases) {
        // Check that the interpolation order within any phases does not exceed the globally requested order.
        if (phase.residualInterpolationOrder > cfg.interpolation_order) {
            string errMsg = format("ERROR: The residual interpolation order in phase %d exceeds the globally selected interpolation order.\n", i);
            errMsg ~= format("       phase interpolation order= %d  globally-requested interpolation order= %d\n",
                    phase.residualInterpolationOrder, cfg.interpolation_order);
            errMsg ~= "       This is not allowed because memory is allocated based on the globally selected interpolation order.\n";
            throw new Error(errMsg);
        }            
        if (phase.jacobianInterpolationOrder > cfg.interpolation_order) {
            string errMsg = format("ERROR: The Jacobian interpolation order in phase %d exceeds the globally selected interpolation order.\n", i);
            errMsg ~= format("       phase interpolation order= %d  globally-requested interpolation order= %d\n",
                    phase.jacobianInterpolationOrder, cfg.interpolation_order);
            errMsg ~= "       This is not allowed because memory is allocated based on the globally selected interpolation order.\n";
            throw new Error(errMsg);
        }            
    }
}

/*---------------------------------------------------------------------
 * Main iteration algorithm
 *---------------------------------------------------------------------
 */

void performNewtonKrylovUpdates(int snapshotStart, double startCFL, int maxCPUs, int threadsPerMPITask)
{
    alias cfg = GlobalConfig;

    readNewtonKrylovConfig();

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
    if (nkCfg.writeLimiterValues) {
        limCIO = new FVCellLimiterIO(buildLimiterVariables());
        if (cfg.flow_format == "rawbinary")
            limBlkIO = new BinaryBlockIO(limCIO);
        else
            limBlkIO = new GzipBlockIO(limCIO);
        limBlkIO.writeMetadataToFile(lmrCfg.limiterMetadataFile);
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
	// Only let master read from disk.
	if (GlobalConfig.is_master_task) {
	    readRestartMetadata();
	    readReferenceResidualsFromFile();
	    // Now ditch snapshots BEYOND what was requested.
	    snapshots.length = snapshotStart;
	}
	// And, in MPI, broadcast information.
	version(mpi_parallel) {
	    broadcastRestartInfo();
	    MPI_Bcast(&referenceGlobalResidual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	    broadcastResiduals(referenceResiduals);
	}
        referenceResidualsAreSet = true;
        nWrittenSnapshots = to!int(snapshots.length);
        RestartInfo restart = snapshots[$-1];
        startStep = restart.step + 1;
        globalResidual = restart.globalResidual;
        prevGlobalResidual = restart.prevGlobalResidual;
        // Determine phase
	if (nkCfg.phaseChangesAtSteps.length == 0) {
	    // Only a single phase
	    setPhaseSettings(0);
	}
	else {
	    foreach (phase, phaseStep; nkCfg.phaseChangesAtSteps) {
		if (startStep < phaseStep) {
		    setPhaseSettings(phase);
		    break;
		}
	    }
	    // end condition when step is past final phase
	    if (startStep >= nkCfg.phaseChangesAtSteps[$-1]) {
		auto finalPhase = nkCfg.phaseChangesAtSteps.length;
		setPhaseSettings(finalPhase-1);
	    }
	}
        if (activePhase.useAutoCFL) {
            cflSelector = new ResidualBasedAutoCFL(activePhase.autoCFLExponent, activePhase.maxCFL,
						   activePhase.thresholdRelativeResidualForCFLGrowth);
            cfl = cflSelector.nextCFL(restart.cfl, startStep, globalResidual, prevGlobalResidual, globalResidual/referenceGlobalResidual);
        }
        else { // Assume we have a global (phase-independent) schedule
            cfl = cflSelector.nextCFL(-1.0, startStep, -1.0, -1.0, -1.0);
        }
    }
    else {
        // On fresh start, the phase setting must be at 0
        setPhaseSettings(0);
        if (activePhase.useAutoCFL) {
            cflSelector = new ResidualBasedAutoCFL(activePhase.autoCFLExponent, activePhase.maxCFL,
						   activePhase.thresholdRelativeResidualForCFLGrowth);
            cfl = activePhase.startCFL;
        }
        else { // Assume we have a global (phase-independent) schedule
            cfl = cflSelector.nextCFL(-1.0, startStep, -1.0, -1.0, -1.0);
        }

        // On fresh start, may need to set reference residuals based on initial condition.
        evalResidual(0);
        setResiduals();
        computeGlobalResidual();
        referenceGlobalResidual = globalResidual;
        computeResiduals(referenceResiduals);
        if (nkCfg.numberOfStepsForSettingReferenceResiduals == 0) {
            referenceResidualsAreSet = true;
            if (GlobalConfig.is_master_task) {
                writeln("*************************************************************************");
                writeln("*");
                writeln("*  At step 0, reference residuals have been set.");
                writefln("*  Reference global residual: %.12e", referenceGlobalResidual);
                writeln("*");
                writeln("*  Reference residuals for each conservation equation:");
                foreach (ivar; 0 .. nConserved) {
                    writefln("* %12s: %.12e", cfg.cqi.names[ivar], referenceResiduals[ivar].re);
                }
                writeln("*");
                writeln("*************************************************************************\n");
                writeReferenceResidualsToFile();
            }
        }   
    } // end actions for fresh start

    // Override CFL if supplied at command line.
    if (startCFL > 0.0) {
	cfl = startCFL;
	if (cfg.verbosity_level > 0 && cfg.is_master_task) {
	    writefln("lmr run-steady: CFL overridden with command-line choice, cfl=%f", cfl);
	}
    }

    // Start timer right at beginning of stepping.
    auto wallClockStart = Clock.currTime();
    double wallClockElapsed;
    int numberBadSteps = 0;
    bool startOfNewPhase = false;

    foreach (step; startStep .. nkCfg.maxNewtonSteps+1) {

        /*---
         * 0. Check for any special actions based on step to perform at START of step
         *---
         *    a. change of phase
         *    b. limiter freezing
	 *    c. compute CFl for this step
         *    d. set the timestep
         *    e. set flag on preconditioner
         */
        residualsUpToDate = false;
        // 0a. change of phase
        auto currentPhase = countUntil(nkCfg.phaseChangesAtSteps, step);
        startOfNewPhase = false;
        if (currentPhase != -1) { // start of new phase detected
	    currentPhase++; // increment because countUntil counts from 0
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
            cfg.frozen_limiter = false;
        }

	// 0c. Compute the CFL for this step
	if (step > startStep) { // because starting step is taken care of specially BEFORE this loop
	    if (numberBadSteps == 0) {
		// then previous update was fine, so proceed to get a new CFL value
		if (activePhase.useAutoCFL) {
		    // We need to be careful in the early steps with the auto CFL.
		    // On step 1, we have no previous residual, so we can't make an adjustment.
		    // Also, we need to check on a residual drop, but this makes no sense
		    // until the reference residuals are established.
		    if (step > nkCfg.numberOfStepsForSettingReferenceResiduals) {
			cfl = cflSelector.nextCFL(cfl, step, globalResidual, prevGlobalResidual, globalResidual/referenceGlobalResidual);
		    }
		}
		else {
		    cfl = cflSelector.nextCFL(-1.0, step, -1.0, -1.0, -1.0);
		}
            }
	    else {
		// Presume last step was bad, so we reduce the CFL (if doing auto CFL)
		cfl *= nkCfg.cflReductionFactor;
		if (cfl <= nkCfg.cflMin) {
		    writeln("The CFL has been reduced due a bad step, but now it has dropped below the minimum allowable CFL.");
		    writefln("current cfl = %e  \t minimum allowable cfl = %e", cfl, nkCfg.cflMin);
		    writeln("Bailing out!");
		    exit(1);
		}
	    }
	}

        // 0d. Set the timestep for this step
        dt = setDtInCells(cfl, activePhase.useLocalTimestep);

        // 0e. determine if we need to update preconditioner

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
        prevGlobalResidual = globalResidual;
        solveNewtonStep(updatePreconditionerThisStep);
	computeGlobalResidual();

	/* 1a. perform a physicality check if required */
        double omega = nkCfg.usePhysicalityCheck ? determineRelaxationFactor() : 1.0;

	/* 1b. do a line search if required */
	if ( (omega > nkCfg.minRelaxationFactor) && nkCfg.useLineSearch ) {
	    omega = applyLineSearch(omega);
	}

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
                writefln("step= %d, relaxation factor= %f", step, omega);
                writeln("Bailing out!");
                exit(1);
            }
        }

        /*---
         * 2. Post-update actions
         *---
         * Here we need to do some house-keeping and see if we continue with iterations.
         */

        /*----
         * 2a. Set reference residuals if needed
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

		    writeReferenceResidualsToFile();
                }
            }
        }

        // [TODO] Add in a halt_now condition.

        /*---
         * 2c. Reporting (to files and screen)
         *---
         */
        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();
        if (((step % nkCfg.stepsBetweenDiagnostics) == 0) || (finalStep && nkCfg.writeDiagnosticsOnLastStep)) {
            writeDiagnostics(step, dt, cfl, wallClockElapsed, omega, residualsUpToDate);
        }

        if (((step % nkCfg.stepsBetweenSnapshots) == 0) || (finalStep && nkCfg.writeSnapshotOnLastStep)) {
            writeSnapshot(step, dt, cfl, nWrittenSnapshots);
            if (nkCfg.writeLimiterValues) {
                writeLimiterValues(step, nWrittenSnapshots);
            }
        }

	version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }


        // [TODO] Write loads. We only need one lot of loads.
        // Any intermediate loads before steady-state have no physical meaning.
        // They might have some diagnostic purpose?


	/*---
         * 2b. Stopping checks.
         *---
         */
        string reasonForStop;
        if (step == nkCfg.maxNewtonSteps) {
            finalStep = true;
            if (cfg.is_master_task) {
                writeln("*** STOPPING: Reached maximum number of steps.");
                reasonForStop = "STOP-REASON: maximum-steps";
            }
        }
        if (!activePhase.ignoreStoppingCriteria) {
            if (globalResidual <= nkCfg.stopOnAbsoluteResidual) {
                finalStep = true;
                if (cfg.is_master_task) {
                    writeln("*** STOPPING: The absolute global residual is below target value.");
                    writefln("         current global residual= %.6e  target value= %.6e", globalResidual, nkCfg.stopOnAbsoluteResidual);
                    reasonForStop = "STOP-REASON: absolute-global-residual-target";
                }
            }
            if ((globalResidual/referenceGlobalResidual) <= nkCfg.stopOnRelativeResidual) {
                finalStep = true;
                if (cfg.is_master_task) {
                    writeln("*** STOPPING: The relative global residual is below target value.");
                    writefln("              current residual= %.6e  target value= %.6e", (globalResidual/referenceGlobalResidual), nkCfg.stopOnRelativeResidual);
                    reasonForStop = "STOP-REASON: relative-global-residual-target";
                }
            }
        }

        // Reporting to screen on progress.
        if (((step % nkCfg.stepsBetweenStatus) == 0) || finalStep) {
            printStatusToScreen(step, cfl, dt, wallClockElapsed, residualsUpToDate);
        }

        if (finalStep) {
	    if (cfg.is_master_task) {
                writeln(reasonForStop);
		writefln("FINAL-STEP: %d", step);
		writefln("FINAL-CFL: %.3e", cfl);
	    }
	    break;
	}
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
        case PreconditionerType.diagonal:
            // do nothing
            break;
        } // end switch
    }
}

void writeReferenceResidualsToFile()
{
    auto f = File(lmrCfg.referenceResidualsFile, "w");
    f.writef("%.18e", referenceGlobalResidual);
    foreach (r; referenceResiduals) {
	f.writef(" %.18e", r.re);
    }
    f.writef("\n");
    f.close();
}

void readReferenceResidualsFromFile()
{
    size_t nConserved = GlobalConfig.cqi.n;
    auto f = File(lmrCfg.referenceResidualsFile, "r");
    auto line = f.readln().strip();
    auto tokens = line.split();
    referenceGlobalResidual = to!double(tokens[0]);
    size_t startIdx = 1;
    foreach (ivar; 0 .. nConserved) {
        referenceResiduals[ivar] = to!double(tokens[startIdx+ivar]);
    }
    f.close();
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
                MPI_Reduce(&(residuals[ivar].re), &(residuals[ivar].re), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
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
double setDtInCells(double cfl, bool useLocalTimestep)
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
    if (!useLocalTimestep) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (cell; blk.cells) cell.dt_local = dt;
        }
    }

    return dt;
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
void solveNewtonStep(bool updatePreconditionerThisStep)
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
    evalResidual(0);

    setResiduals();

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

    if (nkCfg.useResidualSmoothing) {
	applyResidualSmoothing();
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
void computeGlobalResidual()
{
    mixin(dotOverBlocks("globalResidual", "R", "R"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    globalResidual = sqrt(globalResidual);
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
 * Apply residual smoothing.
 *
 * Reference:
 * Mavriplis (2021)
 * A residual smoothing strategy for accelerating Newton method continuation.
 * Computers & Fluids
 *
 * Authors: KAD and RJG
 * Date: 2023-03-13
 */
void applyResidualSmoothing()
{
    size_t nConserved = GlobalConfig.cqi.n;

    //----
    // 1. Compute approximate solution via dU = D^{-1} * R(U) where D = preconditioner matrix
    //----

    final switch (nkCfg.preconditioner) {
    case PreconditionerType.diagonal:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = to!number(0.0); }
        mixin(diagonal_solve("DinvR", "R"));
        break;
    case PreconditionerType.jacobi:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = to!number(0.0); }
        mixin(jacobi_solve("DinvR", "R"));
        break;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = to!number(0.0); }
        mixin(sgs_solve("DinvR", "R"));
        break;
    case PreconditionerType.ilu:
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.DinvR[] = blk.R[];
            nm.smla.solve(blk.flowJacobian.local, blk.DinvR);
        }
        break;
    } // end switch

    //----
    // 2. Add smoothing source term to RHS
    //----
    foreach (blk; parallel(localFluidBlocks,1)) {
	size_t startIdx = 0;
	foreach (cell; blk.cells) {
	    double dtInv = 1.0/cell.dt_local;
	    blk.R[startIdx .. startIdx+nConserved] += dtInv * blk.DinvR[startIdx .. startIdx+nConserved];
	    startIdx += nConserved;
	}
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
 * Authors: KAD and RJG
 * Date: 2022-07-09
 */
void applyPreconditioning()
{
    auto nConserved = GlobalConfig.cqi.n;

    final switch (nkCfg.preconditioner) {
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
        if (nkCfg.useRealValuedFrechetDerivative) {
            evalRealMatVecProd(sigma);
        } else {
            evalComplexMatVecProd(sigma);
        }
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
        size_t startIdx = 0;
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
    foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(activePhase.residualInterpolationOrder);
}

/**
 * Determine a relaxation factor based on a physicality check.
 *
 * In this algorithm, the relaxation factor keeps falling as we search across cells in order.
 * This is efficient since we're searching for a worst case. If at any point, the relaxation
 * factor gets smaller than what we're prepared to accept then we just break the search
 * in that block of cells.
 *
 * [TODO:KAD] Add reference please
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
        size_t startIdx = 0;
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
 * Apply a line search to modify relaxation of factor.
 *
 * The line search is used *after* the physicality check is performed
 * to ensure that the unsteady residual is reduced in this step.
 *
 * [TODO:KAD] Add reference please
 *
 * Authors: KAD and RJG
 * Date: 2023-03-13
 */
double applyLineSearch(double omega) {

    double minOmega = nkCfg.minRelaxationFactor;
    double omegaReductionFactor = nkCfg.relaxationFactorReductionFactor;

    size_t nConserved = GlobalConfig.cqi.n;

    double RU0 = globalResidual;
    double RUn;

    bool reduceOmega = true;
    while (reduceOmega) {

	//----
	// 1. Compute unsteady term
	//----
	foreach (blk; parallel(localFluidBlocks,1)) {
	    size_t startIdx = 0;
	    foreach (cell; blk.cells) {
		blk.R[startIdx .. startIdx+nConserved] = -cell.dt_local * omega * blk.dU[startIdx .. startIdx+nConserved];
		startIdx += nConserved;
	    }
	}

	//----
	// 2. Compute residual at updated state
	//----
	foreach (blk; parallel(localFluidBlocks,1)) {
	    size_t startIdx = 0;
	    foreach (cell; blk.cells) {
		cell.U[1].copy_values_from(cell.U[0]);
		foreach (ivar; 0 .. nConserved) {
		    cell.U[1][ivar] = cell.U[0][ivar] + omega * blk.dU[startIdx+ivar];
		}
		cell.decode_conserved(0, 1, 0.0);
		startIdx += nConserved;
	    }
	}
	evalResidual(1);
	foreach (blk; parallel(localFluidBlocks,1)) {
	    size_t startIdx = 0;
	    foreach (cell; blk.cells) {
		blk.R[startIdx .. startIdx+nConserved] += cell.dUdt[1][0 .. nConserved];
		// return cell to original state
		cell.decode_conserved(0, 0, 0.0);
		startIdx += nConserved;
	    }
	}

	//----
	// 3. Add smoothing source term
	//----
	if (nkCfg.useResidualSmoothing) {
	    foreach (blk; parallel(localFluidBlocks,1)) {
		size_t startIdx = 0;
		foreach (cell; blk.cells) {
		    auto dtInv = 1.0/cell.dt_local;
		    blk.R[startIdx .. startIdx+nConserved] += dtInv * blk.DinvR[startIdx .. startIdx+nConserved];
		}
		startIdx += nConserved;
	    }
	}

	//----
	// 4. Compute norm of unsteady residual
	//----
	mixin(dotOverBlocks("RUn", "R", "R"));
	version(mpi_parallel) {
	    MPI_Allreduce(MPI_IN_PLACE, &RUn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
	RUn = sqrt(RUn);

	//----
	// 5. Check if unsteady residual is reduced
	//----
	if ( (RUn < RU0) || (omega < minOmega) ) {
	    // we are done, don't attempt to reduce omega any further
	    reduceOmega = false;
	}
	else {
	    omega *= omegaReductionFactor;
	}
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
    ensure_directory_is_present(diagDir);
    diagnostics = File(diagFile, "w");

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

    diagnostics = File(diagFile, "a");
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

void writeSnapshot(int step, double dt, double cfl, ref int nWrittenSnapshots)
{
    alias cfg = GlobalConfig;
    size_t nConserved = cfg.cqi.n;
    if (cfg.is_master_task) {
        writeln();
        writefln("+++++++++++++++++++++++++++++++++++++++");
        writefln("+   Writing snapshot at step = %4d   +", step);
        writefln("+++++++++++++++++++++++++++++++++++++++");
    }

    if (!residualsUpToDate) {
	computeResiduals(currentResiduals);
	residualsUpToDate = true;
    }

    nWrittenSnapshots++;

    if (nWrittenSnapshots <= nkCfg.totalSnapshots) {
	// Write snapshot
        auto dirName = snapshotDirectory(nWrittenSnapshots);
        if (cfg.is_master_task) {
            ensure_directory_is_present(dirName);
        }

        // Wait for master to complete building the directory for the next snapshot
        version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

        foreach (blk; localFluidBlocks) {
	    auto fileName = flowFilename(nWrittenSnapshots, blk.id);
	    blkIO.writeVariablesToFile(fileName, blk.cells);
        }

	// Add restart info
	snapshots ~= RestartInfo(nConserved);
	snapshots[$-1].step = step;
	snapshots[$-1].dt = dt;
	snapshots[$-1].cfl = cfl;
	snapshots[$-1].globalResidual = globalResidual;
	snapshots[$-1].prevGlobalResidual = prevGlobalResidual;
	snapshots[$-1].residuals = currentResiduals.dup;

    }
    else {
        // We need to shuffle all of the snapshots
        foreach (iSnap; 2 .. nkCfg.totalSnapshots+1) {
            foreach (blk; localFluidBlocks) {
		auto fromName = flowFilename(iSnap, blk.id);
                auto toName = flowFilename(iSnap-1, blk.id);
                rename(fromName, toName);
            }
        }
        foreach (blk; localFluidBlocks) {
	    auto fileName = flowFilename(nkCfg.totalSnapshots, blk.id);
	    blkIO.writeVariablesToFile(fileName, blk.cells);
        }

	// Shuffle the restart info
	snapshots[0 .. nkCfg.totalSnapshots-1] = snapshots[1 .. nkCfg.totalSnapshots];
	// Add new info at end
	snapshots[$-1].step = step;
	snapshots[$-1].dt = dt;
	snapshots[$-1].cfl = cfl;
	snapshots[$-1].globalResidual = globalResidual;
	snapshots[$-1].prevGlobalResidual = prevGlobalResidual;
	snapshots[$-1].residuals = currentResiduals.dup;
    }

    if (GlobalConfig.is_master_task) {
	writeRestartMetadata(snapshots);
    }
}

void writeRestartMetadata(RestartInfo[] snapshots)
{
    auto f = File(lmrCfg.restartFile, "w");
    auto timeNow =  cast(DateTime)(Clock.currTime());
    f.writefln("# Restart metadata written at: %s", timeNow.toSimpleString());
    f.writefln("# step,  dt,  cfl, global-residual, pre-global-residual, n-conserved quantities residuals");
    foreach (snap; snapshots) {
	f.writef("%04d %.18e %.18e %.18e %.18e", snap.step, snap.dt, snap.cfl, snap.globalResidual, snap.prevGlobalResidual);
	foreach (r; snap.residuals) f.writef(" %.18e", r.re);
	f.writef("\n");
    }
    f.close();
}

void readRestartMetadata()
{
    auto nConserved = GlobalConfig.cqi.n;
    auto f = File(lmrCfg.restartFile, "r");
    auto line = f.readln().strip();
    while (line.length > 0) {
	if (line[0] != '#') {
	    // Process non-comment line
	    auto tks = line.split();
	    snapshots ~= RestartInfo(nConserved);
	    snapshots[$-1].step = to!int(tks[0]);
	    snapshots[$-1].dt = to!double(tks[1]);
	    snapshots[$-1].cfl = to!double(tks[2]);
	    snapshots[$-1].globalResidual = to!double(tks[3]);
	    snapshots[$-1].prevGlobalResidual = to!double(tks[4]);
	    size_t start_idx = 5;
	    foreach (i, ref r; snapshots[$-1].residuals) r = to!double(tks[start_idx+i]);
	}
	line = f.readln().strip();
    }
    f.close();
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

/**
 * Write limiter values to disk.
 *
 * Presently, the only use case for this is on the final step.
 * We write the limiter values so that the adjoint solver can
 * read them in at initialisation time.
 *
 * Authors: RJG
 * Date: 2023-08-13
 */
void writeLimiterValues(int step, int nWrittenSnapshots)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("    |");
        writefln("    |-->  Writing limiter values at step = %4d  ", step);
    }

    int iSnap = (nWrittenSnapshots <= nkCfg.totalSnapshots) ? nWrittenSnapshots : nkCfg.totalSnapshots;
    foreach (blk; localFluidBlocks) {
        auto fileName = limiterFilename(iSnap, blk.id);
        limBlkIO.writeVariablesToFile(fileName, blk.cells);
    }
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

