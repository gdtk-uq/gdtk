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
import std.range : walkLength;
import std.datetime : DateTime, Clock;
import std.parallelism : parallel, defaultPoolThreads;
import std.stdio : File, writeln, writefln, stdout;
import std.file: copy, rename, dirEntries, SpanMode, DirEntry, readText, exists;
import std.array : appender;
import std.format : formattedWrite;
import std.json : JSONValue;
import std.math;
import std.string;
import std.conv : to;
import std.typecons : Tuple, tuple;

import ntypes.complex;
import nm.number : number;
import nm.smla;
import nm.bbla;
import util.time_utils : timeStringToSeconds;
import util.lua;
import util.lua_service;
import lua_helper;
import util.json_helper;
import geom;

import lmrexceptions;
import lmrconfig;
import conservedquantities : ConservedQuantities, copy_values_from;
import fileutil : ensure_directory_is_present;

import globalconfig;
import globaldata;
import init;
import simcore: compute_mass_balance;
import simcore_gasdynamic_step : detect_shocks;
import simcore_exchange;
import bc;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;
import user_defined_source_terms : getUDFSourceTermsForCell;
import blockio;
import fvcellio;
import lmr.fvcell : FVCell;
import lmr.loads : writeLoadsToFile;
import lmr.loads : init_current_loads_indx_dir,
                   wait_for_current_indx_dir,
                   writeLoadsToFile,
                   count_written_loads,
                   update_loads_metadata_file;
import lmr.lmrwarnings;

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

FVCellIO residCIO;
BlockIO residBlkIO;

FVCellIO gradCIO;
BlockIO gradBlkIO;

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
    int[] maxStepsInInitialPhases;
    double[] phaseChangesAtRelativeResidual;
    // Newton stepping and continuation
    bool inviscidCFLOnly = true;
    bool useLineSearch = true;
    int lineSearchOrder = 1;
    bool usePhysicalityCheck = true;
    double allowableRelativeMassChange = 0.2;
    double minRelaxationFactorForUpdate = 0.01;
    double minRelaxationFactorForCFLGrowth = 0.1;
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
    bool writeLoadsOnLastStep = true;
    bool writeLimiterValues = false;
    bool writeResidualValues = false;
    bool writeGradientValues = false;
    bool writeLoads = false;

    void readValuesFromJSON(JSONValue jsonData)
    {
        numberOfStepsForSettingReferenceResiduals = getJSONint(jsonData, "number_of_steps_for_setting_reference_residuals", numberOfStepsForSettingReferenceResiduals);
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
        maxStepsInInitialPhases = getJSONintarray(jsonData, "max_steps_in_initial_phases", maxStepsInInitialPhases);
        phaseChangesAtRelativeResidual = getJSONdoublearray(jsonData, "phase_changes_at_relative_residual", phaseChangesAtRelativeResidual);
        inviscidCFLOnly = getJSONbool(jsonData, "inviscid_cfl_only", inviscidCFLOnly);
        useLineSearch = getJSONbool(jsonData, "use_line_search", useLineSearch);
        lineSearchOrder = getJSONint(jsonData, "line_search_order", lineSearchOrder);
        usePhysicalityCheck = getJSONbool(jsonData, "use_physicality_check", usePhysicalityCheck);
        allowableRelativeMassChange = getJSONdouble(jsonData, "allowable_relative_mass_change", allowableRelativeMassChange);
        minRelaxationFactorForUpdate = getJSONdouble(jsonData, "min_relaxation_factor_for_update", minRelaxationFactorForUpdate);
        minRelaxationFactorForCFLGrowth = getJSONdouble(jsonData, "min_relaxation_factor_for_cfl_growth", minRelaxationFactorForCFLGrowth);
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
        writeLoadsOnLastStep = getJSONbool(jsonData, "write_loads_on_last_step", writeLoadsOnLastStep);
        writeLimiterValues = getJSONbool(jsonData, "write_limiter_values", writeLimiterValues);
        writeResidualValues = getJSONbool(jsonData, "write_residual_values", writeResidualValues);
        writeGradientValues = getJSONbool(jsonData, "write_gradient_values", writeGradientValues);
        writeLoads = getJSONbool(jsonData, "write_loads", writeLoads);
    }
}
NKGlobalConfig nkCfg;

struct NKPhaseConfig {
    bool useLocalTimestep = true;
    int residualInterpolationOrder = 2;
    int jacobianInterpolationOrder = 2;
    bool frozenPreconditioner = true;
    int stepsBetweenPreconditionerUpdate = 10;
    bool enforceLinearSolverTolerance = false;
    bool useAdaptivePreconditioner = false;
    bool ignoreStoppingCriteria = true;
    bool frozenShockDetector = false;
    bool frozenLimiterForResidual = false;
    bool frozenLimiterForJacobian = false;
    double linearSolveTolerance = 0.01;
    // Auto CFL control
    bool useAutoCFL = false;
    double thresholdRelativeResidualForCFLGrowth = 0.99;
    double startCFL = 1.0;
    double maxCFL = 1000.0;
    double autoCFLExponent = 0.75;
    double limitOnCFLIncreaseRatio = 2.0;
    double limitOnCFLDecreaseRatio = 0.1;

    void readValuesFromJSON(JSONValue jsonData)
    {
        useLocalTimestep = getJSONbool(jsonData, "use_local_timestep", useLocalTimestep);
        residualInterpolationOrder = getJSONint(jsonData, "residual_interpolation_order", residualInterpolationOrder);
        jacobianInterpolationOrder = getJSONint(jsonData, "jacobian_interpolation_order", jacobianInterpolationOrder);
        frozenPreconditioner = getJSONbool(jsonData, "frozen_preconditioner", frozenPreconditioner);
        stepsBetweenPreconditionerUpdate = getJSONint(jsonData, "steps_between_preconditioner_update", stepsBetweenPreconditionerUpdate);
        enforceLinearSolverTolerance = getJSONbool(jsonData, "enforce_linear_solver_tolerance", enforceLinearSolverTolerance);
        useAdaptivePreconditioner = getJSONbool(jsonData, "use_adaptive_preconditioner", useAdaptivePreconditioner);
        ignoreStoppingCriteria = getJSONbool(jsonData, "ignore_stopping_criteria", ignoreStoppingCriteria);
        frozenShockDetector = getJSONbool(jsonData, "frozen_shock_detector", frozenShockDetector);
        frozenLimiterForResidual = getJSONbool(jsonData, "frozen_limiter_for_residual", frozenLimiterForResidual);
        frozenLimiterForJacobian = getJSONbool(jsonData, "frozen_limiter_for_jacobian", frozenLimiterForJacobian);
        linearSolveTolerance = getJSONdouble(jsonData, "linear_solve_tolerance", linearSolveTolerance);
        useAutoCFL = getJSONbool(jsonData, "use_auto_cfl", useAutoCFL);
        thresholdRelativeResidualForCFLGrowth = getJSONdouble(jsonData, "threshold_relative_residual_for_cfl_growth", thresholdRelativeResidualForCFLGrowth);
        startCFL = getJSONdouble(jsonData, "start_cfl", startCFL);
        maxCFL = getJSONdouble(jsonData, "max_cfl", maxCFL);
        autoCFLExponent = getJSONdouble(jsonData, "auto_cfl_exponent", autoCFLExponent);
        limitOnCFLIncreaseRatio = getJSONdouble(jsonData, "limit_on_cfl_increase_ratio", limitOnCFLIncreaseRatio);
        limitOnCFLDecreaseRatio = getJSONdouble(jsonData, "limit_on_cfl_decrease_ratio", limitOnCFLDecreaseRatio);
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
    this(double p, double cfl_max, double thresholdResidualDrop, double limit_on_cfl_increase_ratio, double limit_on_cfl_decrease_ratio)
    {
        mP = p;
        mMaxCFL = cfl_max;
        mThresholdResidualDrop = thresholdResidualDrop;
        mLimitOnCFLIncreaseRatio = limit_on_cfl_increase_ratio;
        mLimitOnCFLDecreaseRatio = limit_on_cfl_decrease_ratio;
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
    double mLimitOnCFLIncreaseRatio;
    double mLimitOnCFLDecreaseRatio;
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

ConservedQuantities referenceResiduals, currentResiduals, rowScale, colScale;
double referenceGlobalResidual, globalResidual, prevGlobalResidual;
bool residualsUpToDate = false;
bool referenceResidualsAreSet = false;


// Module-local, global memory arrays and matrices
// TODO: Think about these, maybe they shouldn't be globals
double[] g0;
double[] g1;
double[] h;
double[] hR;
Matrix!double H0;
Matrix!double H1;
Matrix!double Gamma;
Matrix!double Q0;
Matrix!double Q1;

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
    int phase;
    int stepsIntoPhase;

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
        MPI_Bcast(&(snapshots[i].phase), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(snapshots[i].stepsIntoPhase), 1, MPI_INT, 0, MPI_COMM_WORLD);
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

    // Initialise baseline configuration
    initConfiguration();
    // After loading configuration, we can check if there are any settings
    // that are not particularly appropriate for the steady-state solver.
    if (cfg.extrema_clipping) {
        writeln(warningMessage(LmrWarning.useOfExtremaClipping));
    }
    if (cfg.nFluidBlocks == 0 && cfg.is_master_task) {
        throw new NewtonKrylovException("No FluidBlocks; no point in continuing with simulation initialisation.");
    }
    cfg.n_flow_time_levels = 2;

    initLocalBlocks();

    /* RJG, 2024-03-05
     * The shared memory model is not playing nicely with how the Newton-Krylov module
     * is set up. For the present, we limit the threads to 1.
     */
    initThreadPool(1, 1);

    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

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

    if ((cfg.interpolation_order > 1) &&
        ((cfg.unstructured_limiter == UnstructuredLimiter.hvenkat) ||
         (cfg.unstructured_limiter == UnstructuredLimiter.venkat) ||
         (cfg.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp) ||
         (cfg.unstructured_limiter == UnstructuredLimiter.venkat_mlp))) {
        initUSGlimiters();
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

    // we enforce that the user always has the correct number of entries in the max_steps_in_initial_phases list
    if (cfg.is_master_task && nkCfg.maxStepsInInitialPhases.length != nkCfg.numberOfPhases-1) {
        string errMsg;
        errMsg ~= "\n";
        errMsg ~= format("ERROR: number of phases = %d but number of entries in max_steps_in_initial_phases list is: %d\n",
                         nkCfg.numberOfPhases, nkCfg.maxStepsInInitialPhases.length);
        errMsg ~= "       We expect  (number of phases) - 1 entries in max_steps_in_initial_phases list.\n";
        errMsg ~= "       These entries are the maximum number of steps that will be taken in each initial phase,\n";
        errMsg ~= "       note that we do not need an entry for the terminal phase.\n";
        throw new Error(errMsg);
    }

    // we will allow the phase_changes_at_relative_residual list to be empty, if it is empty, we will only key the phase changes off of the steps
    if (nkCfg.phaseChangesAtRelativeResidual.length == 0) {
        size_t nentries = nkCfg.maxStepsInInitialPhases.length;
        nkCfg.phaseChangesAtRelativeResidual.length = nentries;
        foreach (ref val; nkCfg.phaseChangesAtRelativeResidual) { val = 0.0; }
    }

    // if the user has defined a nonempty phase_changes_at_relative_residual list, we enforce that the user has the correct number of entries
    if (cfg.is_master_task && nkCfg.phaseChangesAtRelativeResidual.length != nkCfg.maxStepsInInitialPhases.length ) {
        string errMsg;
        errMsg ~= "\n";
        errMsg ~= format("ERROR: number of phases = %d but number of entries in phase_changes_at_relative_residual list is: %d\n",
                         nkCfg.numberOfPhases, nkCfg.phaseChangesAtRelativeResidual.length);
        errMsg ~= "       We expect at least (number of phases) - 1 entries in phase_changes_at_relative_residual list.\n";
        errMsg ~= "       These entries are the relative residual at which to change from one phase to the next.\n";
        throw new Error(errMsg);
    }

    foreach (i, phase; nkPhases) {
        // Check that the interpolation order within any phases does not exceed the globally requested order.
        if (cfg.is_master_task && phase.residualInterpolationOrder > cfg.interpolation_order) {
            string errMsg;
            errMsg ~= "\n";
            errMsg ~= format("ERROR: The residual interpolation order in phase %d exceeds the globally selected interpolation order.\n", i+1);
            errMsg ~= format("       phase interpolation order= %d  globally-requested interpolation order= %d\n",
                    phase.residualInterpolationOrder, cfg.interpolation_order);
            errMsg ~= "       This is not allowed because memory is allocated based on the globally selected interpolation order.\n";
            throw new Error(errMsg);
        }
        if (cfg.is_master_task && phase.jacobianInterpolationOrder > cfg.interpolation_order) {
            string errMsg;
            errMsg ~= "\n";
            errMsg ~= format("ERROR: The Jacobian interpolation order in phase %d exceeds the globally selected interpolation order.\n", i+1);
            errMsg ~= format("       phase interpolation order= %d  globally-requested interpolation order= %d\n",
                    phase.jacobianInterpolationOrder, cfg.interpolation_order);
            errMsg ~= "       This is not allowed because memory is allocated based on the globally selected interpolation order.\n";
            throw new Error(errMsg);
        }

        // Check for a nonsensical combination of residual and jacobian interpolation orders.
        if (cfg.is_master_task && phase.jacobianInterpolationOrder > phase.residualInterpolationOrder) {
            string errMsg;
            errMsg ~= "\n";
            errMsg ~= format("ERROR: The Jacobian interpolation order in phase %d exceeds the Residual interpolation order.\n", i+1);
            errMsg ~= format("       Jacobian interpolation order= %d  Residual interpolation order= %d\n",
                    phase.jacobianInterpolationOrder,  phase.residualInterpolationOrder);
            errMsg ~= "       This is not allowed, the Jacobian interpolation order should be <= the Residual interpolation order.\n";
            throw new Error(errMsg);
        }

        // Check for nonsensical frozen limiter settings
        if (cfg.is_master_task && phase.frozenLimiterForResidual && phase.frozenLimiterForJacobian != phase.frozenLimiterForResidual) {
            string errMsg;
            errMsg ~= "\n";
            errMsg ~= format("ERROR: incompatible settings in phase %d,\n", i+1);
            errMsg ~= "       frozen_limiter_for_jacobian cannot be false when frozen_limiter_for_residual is true.\n";
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
    int nWrittenLoads = count_written_loads();
    int startStep = 1;
    bool finalStep = false;
    double cfl;
    double dt;
    int currentPhase = 0;
    int stepsIntoCurrentPhase = 0;
    bool updatePreconditionerThisStep = false;
    CFLSelector cflSelector;

    /*----------------------------------------------
     * Initialisation
     *----------------------------------------------
     */
    if (nkCfg.usePreconditioner) initPreconditioner();
    size_t nConserved = cfg.cqi.n;
    referenceResiduals = new ConservedQuantities(nConserved);
    currentResiduals = new ConservedQuantities(nConserved);
    rowScale = new ConservedQuantities(nConserved);
    colScale = new ConservedQuantities(nConserved);
    if (nkCfg.writeLoads && (nWrittenLoads == 0)) {
        if (cfg.is_master_task) {
            initLoadsFiles();
        }
    }
    if (nkCfg.writeLimiterValues) {
        limCIO = new FluidFVCellLimiterIO(buildLimiterVariables(localFluidBlocks[0].grid_type));
        if (cfg.field_format == "rawbinary")
            limBlkIO = new BinaryBlockIO(limCIO);
        else
            limBlkIO = new GzipBlockIO(limCIO);
        limBlkIO.writeMetadataToFile(lmrCfg.limiterMetadataFile);
    }
    if (nkCfg.writeResidualValues) {
        residCIO = new FluidFVCellResidualIO(buildResidualVariables());
        if (cfg.field_format == "rawbinary")
            residBlkIO = new BinaryBlockIO(residCIO);
        else
            residBlkIO = new GzipBlockIO(residCIO);
        residBlkIO.writeMetadataToFile(lmrCfg.residualMetadataFile);
    }
    if (nkCfg.writeGradientValues) {
        gradCIO = new FluidFVCellGradientIO(buildGradientVariables(localFluidBlocks[0].grid_type));
        if (cfg.field_format == "rawbinary")
            gradBlkIO = new BinaryBlockIO(gradCIO);
        else
            gradBlkIO = new GzipBlockIO(gradCIO);
        gradBlkIO.writeMetadataToFile(lmrCfg.gradientMetadataFile);
    }
    allocateGlobalGMRESWorkspace();
    foreach (blk; localFluidBlocks) {
        blk.allocate_GMRES_workspace(nkCfg.maxLinearSolverIterations, nkCfg.useRealValuedFrechetDerivative);
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

    // At this point we should have all our large data structures initialized, so compute some memory usage for reporting
    auto myStats = GC.stats();
    double heapUsed = to!double(myStats.usedSize)/(2^^20);
    double heapFree = to!double(myStats.freeSize)/(2^^20);
    double minTotal = heapUsed+heapFree;
    double maxTotal = heapUsed+heapFree;
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE,&heapUsed,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&heapFree,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&minTotal,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&maxTotal,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }
    if (cfg.is_master_task) {
        writefln("Heap memory used: %.0f MB, unused: %.0f MB, total: %.0f MB (%.0f-%.0f MB per task)",
                 heapUsed, heapFree, heapUsed+heapFree, minTotal, maxTotal);
        stdout.flush();
    }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

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
        currentPhase = restart.phase;
        stepsIntoCurrentPhase = restart.stepsIntoPhase;
        setPhaseSettings(currentPhase);
        if (activePhase.useAutoCFL) {
            cflSelector = new ResidualBasedAutoCFL(activePhase.autoCFLExponent, activePhase.maxCFL,
                                                   activePhase.thresholdRelativeResidualForCFLGrowth,
                                                   activePhase.limitOnCFLIncreaseRatio, activePhase.limitOnCFLDecreaseRatio);
            cfl = cflSelector.nextCFL(restart.cfl, startStep, globalResidual, prevGlobalResidual, globalResidual/referenceGlobalResidual);
        }
        else { // Assume we have a global (phase-independent) schedule
            cfl = cflSelector.nextCFL(-1.0, startStep, -1.0, -1.0, -1.0);
        }
        // On restart, we need to do some diagonstics file housekeeping
        if (cfg.is_master_task) {
            // we first ensure the diagnostics directory and file exist before proceeding
            ensure_directory_is_present(diagDir);
            if (!exists(diagFile)) {
                writeln("Bailing out!...unable to find the Newton-Krylov solver diagnostics file: ", diagFile);
                exit(1);
            }
            // we next make a copy of the current diagnostics file,
            // being careful not to overwrite any previously saved files
            auto numberOfExistingDiagFiles = walkLength(dirEntries(diagDir, "*", SpanMode.shallow));
            string diagFileCopy = diagFile;
            foreach (i; 0..numberOfExistingDiagFiles) {
                diagFileCopy ~= ".save";
            }
            copy(diagFile, diagFileCopy);
            // we now wind the base diagnostics file back to the restart step
            auto outputFile = File(diagFile, "w");
            auto inputFile = File(diagFileCopy, "r");
            foreach (line; inputFile.byLine) {
                if (line.split[0] == to!string(startStep)) break;
                outputFile.writeln(line);
            }
            outputFile.close();
        }
        if (cfg.is_master_task) {
            writeln("*** RESTARTING SIMULATION ***");
            writefln("RESTART-SNAPSHOT: %d", snapshots.length);
            writefln("RESTART-STEP: %d", startStep);
        }
    }
    else {
        // On fresh start, we need to create the diagnostics file
        if (cfg.is_master_task) {
            initialiseDiagnosticsFile();
        }
        // On fresh start, the phase setting must be at 0
        setPhaseSettings(0);
        if (activePhase.useAutoCFL) {
            cflSelector = new ResidualBasedAutoCFL(activePhase.autoCFLExponent, activePhase.maxCFL,
                                                   activePhase.thresholdRelativeResidualForCFLGrowth,
                                                   activePhase.limitOnCFLIncreaseRatio, activePhase.limitOnCFLDecreaseRatio);
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
        // Add value of 1.0 to each residaul.
        // If values are very large, 1.0 makes no difference.
        // If values are zero, the 1.0 should mean the reference residual
        // asymptotes to an absolute residual.
        referenceResiduals[] = referenceResiduals[] + to!number(1.0);

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
                    writefln("* %15s: %.12e", cfg.cqi.names[ivar], referenceResiduals[ivar].re);
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
    double omega = 1.0;
    bool terminalPhase = false;
    if (currentPhase == nkCfg.numberOfPhases-1) terminalPhase = true; // we are beginning in the terminal phase
    if (cfg.is_master_task) {
        writefln("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        writefln("+   Starting in Phase %d (out of %d) at step = %6d; global relative residual = %10.6e   +",
                 currentPhase+1, nkCfg.numberOfPhases, startStep, globalResidual/referenceGlobalResidual);
        writefln("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        writeln();
    }

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
        stepsIntoCurrentPhase++;
        startOfNewPhase = false;
        if (!terminalPhase &&
            (stepsIntoCurrentPhase >= nkCfg.maxStepsInInitialPhases[currentPhase] ||
             globalResidual/referenceGlobalResidual <= nkCfg.phaseChangesAtRelativeResidual[currentPhase])) {
            // start of new phase detected
            startOfNewPhase = true;
            currentPhase++;
            stepsIntoCurrentPhase = 0;
            setPhaseSettings(currentPhase);
            if (currentPhase == nkCfg.numberOfPhases-1) terminalPhase = true;
            if (activePhase.useAutoCFL) {
                // If the user gives us a positive startCFL, use that.
                // Otherwise we continue with the CFL we have. In other words, do nothing special here.
                if (activePhase.startCFL > 0.0) cfl = activePhase.startCFL;
            }
            if (cfg.is_master_task) {
                writefln("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                writefln("+   Start of Phase %d (out of %d) at step = %6d; global relative residual = %10.6e   +",
                         currentPhase+1, nkCfg.numberOfPhases, step, globalResidual/referenceGlobalResidual);
                if (activePhase.startCFL > 0.0) {
                    writefln("+   ---> CFL reset for the start of this phase: cfl=%10.3e                              +", cfl);
                }
                else {
                    writefln("+   ---> CFL continues from previous phase end: cfl=%10.3e                              +", cfl);
                }
                writefln("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
                writeln();
            }
        }

	// 0b. Compute the CFL for this step
	if (step > startStep) { // because starting step is taken care of specially BEFORE this loop
	    if (numberBadSteps == 0) {
		// then previous update was fine, so proceed to get a new CFL value
		if (activePhase.useAutoCFL) {
		    // We need to be careful in the early steps with the auto CFL.
		    // On step 1, we have no previous residual, so we can't make an adjustment.
		    // Also, we need to check on a residual drop, but this makes no sense
		    // until the reference residuals are established.
		    if (step > nkCfg.numberOfStepsForSettingReferenceResiduals &&
                        omega >= nkCfg.minRelaxationFactorForCFLGrowth) { // TODO: consider only limiting CFL growth based on the relaxation factor for the residual-based cflSelector
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

        // 0c. Set the timestep for this step
        dt = setDtInCells(cfl, activePhase.useLocalTimestep);

        // 0d. determine if we need to update preconditioner

        bool maxLinearSolverIterationsUsed = (nkCfg.maxLinearSolverRestarts == gmresInfo.nRestarts &&
                                              nkCfg.maxLinearSolverIterations == gmresInfo.iterationCount);
        if (step == startStep || startOfNewPhase || numberBadSteps > 0 ||
            (step % activePhase.stepsBetweenPreconditionerUpdate) == 0 ||
            (activePhase.useAdaptivePreconditioner && maxLinearSolverIterationsUsed)) {
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

	/* 1a. perform a physicality check if required */
        omega = nkCfg.usePhysicalityCheck ? determineRelaxationFactor() : 1.0;

	/* 1b. do a line search if required */
	if ( (omega > nkCfg.minRelaxationFactorForUpdate) && nkCfg.useLineSearch ) {
	    omega = applyLineSearch(omega);
	}

        /* 1c. check if we achived the allowable linear solver tolerance */
        bool failedToAchieveAllowableLinearSolverTolerance =
            (gmresInfo.finalResidual/gmresInfo.initResidual) > activePhase.linearSolveTolerance &&
            activePhase.enforceLinearSolverTolerance;

        if (omega >= nkCfg.minRelaxationFactorForUpdate && !failedToAchieveAllowableLinearSolverTolerance) {
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
        else if ((omega < nkCfg.minRelaxationFactorForUpdate || failedToAchieveAllowableLinearSolverTolerance) && activePhase.useAutoCFL) {
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
         * 2b. Stopping checks.
         *---
         */
        string reasonForStop;
        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();
        if (step == nkCfg.maxNewtonSteps) {
            finalStep = true;
            if (cfg.is_master_task) {
                writeln("*** STOPPING: Reached maximum number of steps.");
                reasonForStop = "STOP-REASON: maximum-steps";
            }
        }
        if ((SimState.maxWallClockSeconds > 0) && (wallClockElapsed > SimState.maxWallClockSeconds)) {
            finalStep = true;
            if (cfg.is_master_task) {
                writefln("*** STOPPING: Reached maximum wall-clock-time with elapsed time %s", to!string(wallClockElapsed));
                reasonForStop = "STOP-REASON: maximum-wall-clock";
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
        /*---
         * 2c. Reporting (to files and screen)
         *---
         */
        if (((step % nkCfg.stepsBetweenDiagnostics) == 0) || (finalStep && nkCfg.writeDiagnosticsOnLastStep)) {
            writeDiagnostics(step, dt, cfl, wallClockElapsed, omega, currentPhase, residualsUpToDate);
        }

        if (((step % nkCfg.stepsBetweenSnapshots) == 0) || (finalStep && nkCfg.writeSnapshotOnLastStep)) {
            writeSnapshot(step, dt, cfl, currentPhase, stepsIntoCurrentPhase, nWrittenSnapshots);
        }

        if (((step % nkCfg.stepsBetweenLoadsUpdate) == 0) || (finalStep && nkCfg.writeLoadsOnLastStep)) {
            if (nkCfg.writeLoads) {
                writeLoads(step, nWrittenLoads);
            }
        }
        version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

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
    H0 = new Matrix!double(m+1, m);
    H1 = new Matrix!double(m+1, m);
    Gamma = new Matrix!double(m+1, m+1);
    Q0 = new Matrix!double(m+1, m+1);
    Q1 = new Matrix!double(m+1, m+1);
}

/*---------------------------------------------------------------------
 * Auxiliary functions related to iteration algorithm
 *---------------------------------------------------------------------
 */

void setPhaseSettings(size_t phase)
{
    alias cfg = GlobalConfig;
    activePhase = nkPhases[phase];
    foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(activePhase.residualInterpolationOrder);

    // We need to compute the limiter and/or shock detector a final time before freezing it
    if ((activePhase.frozenLimiterForResidual && !cfg.frozen_limiter) ||
        (activePhase.frozenShockDetector && !cfg.frozen_shock_detector)) {
        /*
        debug {
            writeln("DEBUG: setPhaseSettings()");
            writeln(" Evaluating Residual.");
        }
        */
        evalResidual(0);
    }

    if (activePhase.frozenLimiterForResidual) {
        cfg.frozen_limiter = true;
    }
    if (activePhase.frozenShockDetector) {
        cfg.frozen_shock_detector = true;
    }
}

void computeResiduals(ref ConservedQuantities residuals)
{
    size_t nConserved = GlobalConfig.cqi.n;
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (ivar; 0 .. nConserved) blk.residuals[ivar] = fabs(blk.cells[0].dUdt[0][ivar]);

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
        bool gridFlag = blk.grid_type == Grid_t.structured_grid;
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
                signal = cell.signal_frequency(gridFlag);
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
 * In this implementation we opt to precondition first and scale second, the scaled
 * preconditioned linear system is then:
 *
 *     (Dr * J * P^(-1) * Dc) (Dc^(-1) * P * dU) = (Dr * R)
 *
 * where Dr and Dc are the row and column scaling (diagonal) matrices, respectively.
 * J is the Jacobian matrix, P is the precondition matrix, dU is the conserved
 * quantity update vector, and R is the residual vector.
 *
 * We never explicitly store J. We instead apply the scaling and preconditioning in a
 * Jacobian-free manner by following an approach similar to that presented as the
 * SP-GMRES algorithm in:
 *
 *     P. N. Brown and Y. Saad,
 *     Hybrid Krylov methods for nonlinear systems of equations,
 *     SIAM Journal of Scientific and Statistical Computing,
 *     Vol. 11, Iss. 3 (1990)
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

    double linearSolveResidual;
    alias cfg = GlobalConfig;
    size_t nConserved = cfg.cqi.n;

    bool isConverged = false;
    /*---
     * 0. Preparation for iterations.
     *---
     */

    evalResidual(0);

    setResiduals();

    if (nkCfg.useRealValuedFrechetDerivative) { setR0(); }

    computeGlobalResidual();

    determineScaleFactors(rowScale, colScale, nkCfg.useScaling);

    if (nkCfg.usePreconditioner && updatePreconditionerThisStep) {
        computePreconditioner();
    }

    if (nkCfg.useResidualSmoothing) {
        evalResidualSmoothing();
	applyResidualSmoothing();
    }

    setInitialGuess();

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

    double beta, beta0, targetResidual;
    foreach (r; 0 .. nAttempts) {
        // Initialise some working arrays and matrices for this step.
        g0[] = 0.0;
        g1[] = 0.0;
        H0.zeros();
        H1.zeros();
        Gamma.eye();

        // r0 = b - A*x0
        compute_r0();

        // scale the r0 vector
        foreach (blk; parallel(localFluidBlocks,1)) {
            scaleVector(rowScale, blk.r0, nConserved);
        }

        // beta = ||r0||
        beta = computeLinearSystemResidual();

        // Set first residual entry.
        g0[0] = beta;

        // v = r0/beta
        prepareKrylovSpace(beta);

        // Compute target residual on first restart
        if (r == 0) {
            targetResidual = activePhase.linearSolveTolerance * beta;
            beta0 = beta; // store a copy of the initial residual for the diagnostics output
        }

        // Delegate inner iterations
        /*
        debug {
            writeln("DEBUG: solveNewtonStep()");
            writeln(" calling performIterations.");
        }
        */

        isConverged = performIterations(maxIterations, targetResidual, rowScale, colScale, iterationCount, linearSolveResidual);
        int m = iterationCount;

        /*
        debug {
            writeln("DEBUG: solveNewtonStep()");
            writeln(" done performIterations.");
        }
        */

        // At end H := R up to row m
        //        g := gm up to row m
        upperSolve!double(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.transpose_and_dot!double(blk.VT, blk.nvars, m, blk.nvars, blk.g1, blk.zed);
        }

        // At this point we have zed = Dc^(-1) * P * dU
        // we need to prepare the dU values (for Newton update)

        // We first multiply zed by the column scaling to recover zed *= Dc = P * dU
        foreach (blk; parallel(localFluidBlocks,1)) {
            scaleVector(colScale, blk.zed, nConserved);
        }

        // Next we remove preconditioner effect to finally recover dU
        if (nkCfg.usePreconditioner) {
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

        // If we get here, we need to prepare for next restart by setting
        // the initial guess to the current estimate of the solution
        foreach (blk; parallel(localFluidBlocks, 1)) {
            blk.x0[] = blk.dU[];
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
    gmresInfo.finalResidual = linearSolveResidual.re;
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
            foreach (ivar; 0 .. nConserved) {
                blk.R[startIdx+ivar] = cell.dUdt[0][ivar].re;
            }
            startIdx += nConserved;
        }
    }
}

/**
 * Set the unperturbed state R0 used in the real-valued Frechet derivative.
 *
 * Note that this routine assumes that the R array attached to each FluidBlock has been set prior to being called.
 *
 * Authors: RJG and KAD
 * Date: 2024-03-28
 */
void setR0()
{
    size_t nConserved = GlobalConfig.cqi.n;
    bool defectCorrectionUpdate = activePhase.frozenLimiterForJacobian || activePhase.jacobianInterpolationOrder != activePhase.residualInterpolationOrder;
    if (defectCorrectionUpdate) {
        // When performing a defect-correction update with a real-valued Frechet derivative, we cannot use R as our unperturbed state.
        // We must evaluate the Residual vector with the defects applied and store it as R0 for later use in the Fechet derivative.
        foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(activePhase.jacobianInterpolationOrder); }
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = activePhase.frozenLimiterForJacobian; }
        evalResidual(1);
        foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(activePhase.residualInterpolationOrder);
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = activePhase.frozenLimiterForResidual; }
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t startIdx = 0;
            foreach (cell; blk.cells) {
                foreach (ivar; 0 .. nConserved) {
                    blk.R0[startIdx+ivar] = cell.dUdt[1][ivar].re;
                }
                startIdx += nConserved;
            }
        }
    } else {
        // If we aren't performing a defect-correction update, then just set R0 to the already evaluated R vector.
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.R0[] = blk.R[];
        }
    }
}

/**
 * Determine scale factors per conserved quantity.
 *
 * The row scale factors are typically the inverse of the maximum rate of change found
 * globally (over all cells) per each conserved quantity. However, we place some gaurds
 * on determining those scales when rates of change are small.
 *
 * The column scale factors are taken to be the reciprocal of the row scale factors
 * as suggested by Chisholm and Zingg (2009).
 *
 * REFERENCE: T. T. Chisholm and D. W. Zingg,
 *            A Jacobian-Free Newton-Krylov algorithm for compressible turbulent fluid flows,
 *            Journal of Computational Physics,
 *            Vol. 228, Iss. 9 (2009)
 *
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
void determineScaleFactors(ref ConservedQuantities rowScale, ref ConservedQuantities colScale, bool useScaling)
{
    if (!useScaling) {
        rowScale[] = to!number(1.0);
        colScale[] = to!number(1.0);
        return;
    }
    // else we go ahead and determine the scale factors
    rowScale[] = to!number(0.0);
    colScale[] = to!number(0.0);
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
            colScale[ivar] = fmax(colScale[ivar], blk.maxRate[ivar]);
        }
    }
    // In distributed memory, reduce max values and make sure everyone has a copy.
    version(mpi_parallel) {
        foreach (ivar; 0 .. nConserved) {
            MPI_Allreduce(MPI_IN_PLACE, &(colScale[ivar].re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
    }

    // Use a guard on scale values if they get small
    foreach (ivar; 0 .. nConserved) {
        colScale[ivar] = fmax(colScale[ivar], minScaleFactor);
    }

    // The row scaling is the reciprocal of the column scaling
    foreach (ivar; 0 .. nConserved) {
        rowScale[ivar] = 1./colScale[ivar];
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
 * Set initial guess vector, currently we set this to 0 (as is common)
 *
 * x0 = [0]
 *
 * Authors: KAD and RJG
 * Date: 2022-03-02
 */
void setInitialGuess()
{
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.x0[] = 0.0;
    }
}

/**
 * Compute the initial residual vector for GMRES method for any arbitrary x0,
 *
 * r0 = b - A*x0
 *
 * where b is the residual vector R, and A is the augmented Jacobian (I/dt - J).
 *
 * Note that this routine expects the x0 vector to be set prior to calling.
 *
 * Authors: KAD and RJG
 * Date: 2024-02-27
 */
void compute_r0()
{
    size_t nConserved = GlobalConfig.cqi.n;

    // place x0 in zed
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.zed[] = blk.x0[];
    }

    // evaluate A*zed
    evalAugmentedJacobianVectorProduct();

    // set r0 := R - A*zed
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) blk.r0[k] = blk.R[k] - blk.w[k];
    }
}

/**
 * Compute the residual of the linear system from r0.
 *
 * Authors: RJG and KAD
 * Date: 2022-03-02
 */
double computeLinearSystemResidual()
{
    double beta;
    mixin(dotOverBlocks("beta", "r0", "r0"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
void prepareKrylovSpace(double beta)
{
    g0[0] = beta;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.v[] = (1./beta)*blk.r0[];
        foreach (k; 0 .. blk.nvars) {
            blk.VT[k] = blk.v[k];
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

/*
 * Routine to scale a vector
 *
 * Authors: KAD and RJG
 * Date: 2023-03-13
 */

void scaleVector(ConservedQuantities scale, ref double[] vec, size_t nConserved)
{
    foreach (k, ref v; vec) {
        size_t ivar = k % nConserved;
        v *= scale[ivar].re;
    }
}

/**
 * Routines for evaluating and applying residual smoothing.
 *
 * REFERENCE:    D. J. Mavriplis
 *               A residual smoothing strategy for accelerating Newton method continuation.
 *               Computers & Fluids, Vol. 220 (2021)
 *
 * Authors: KAD and RJG
 * Date: 2023-03-13
 */
void evalResidualSmoothing()
{
    // Compute an approximate solution: dU = D^{-1} * R(U),
    // where D is some approximation of the Jacobian, here we use the
    // same approximation as the precondition matrix out of convenience
    size_t nConserved = GlobalConfig.cqi.n;
    final switch (nkCfg.preconditioner) {
    case PreconditionerType.diagonal:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = 0.0; }
        mixin(diagonal_solve("DinvR", "R"));
        break;
    case PreconditionerType.jacobi:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = 0.0; }
        mixin(jacobi_solve("DinvR", "R"));
        break;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = 0.0; }
        mixin(sgs_solve("DinvR", "R"));
        break;
    case PreconditionerType.ilu:
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.DinvR[] = blk.R[];
            nm.smla.solve(blk.flowJacobian.local, blk.DinvR);
        }
        break;
    } // end switch
}

void applyResidualSmoothing()
{
    // Add the smoothing source term to the Residual vector
    size_t nConserved = GlobalConfig.cqi.n;
    foreach (blk; parallel(localFluidBlocks,1)) {
	size_t startIdx = 0;
	foreach (cell; blk.cells) {
	    blk.R[startIdx .. startIdx+nConserved] += (1.0/cell.dt_local) * blk.DinvR[startIdx .. startIdx+nConserved];
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
bool performIterations(int maxIterations, double targetResidual,
                       ref ConservedQuantities rowScale, ref ConservedQuantities colScale,
                       ref int iterationCount,
                       ref double resid)
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


        // 1. apply column scaling to v
        foreach (blk; parallel(localFluidBlocks,1)) {
            scaleVector(colScale, blk.v, nConserved);
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
              writeln(" Jz calc.");
          }
        */

        // 3. Jacobian-vector product
        evalAugmentedJacobianVectorProduct();

        // apply row scaling to w
        foreach (blk; parallel(localFluidBlocks,1)) {
            scaleVector(rowScale, blk.w, nConserved);
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
                size_t offset = i*blk.nvars;
                foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.VT[offset + k];
            }
            double H0_ij;
            mixin(dotOverBlocks("H0_ij", "w", "v"));
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            }
            H0[i,j] = H0_ij;
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k];
            }
        }
        double H0_jp1j;
        mixin(dotOverBlocks("H0_jp1j", "w", "w"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(H0_jp1j.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        H0_jp1j = sqrt(H0_jp1j);
        H0[j+1,j] = H0_jp1j;

        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t offset = (j+1)*blk.nvars;
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.w[k]/H0_jp1j;
                blk.VT[offset + k] = blk.v[k];
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
            nm.bbla.copy(Gamma, Q1);
        }
        else {
            nm.bbla.dot!double(Gamma, j+2, j+2, Q0, j+2, Q1);
        }

        // Prepare for next step
        nm.bbla.copy(H1, H0);
        g0[] = g1[];
        nm.bbla.copy(Q1, Q0);

        // Get residual
        resid = fabs(g1[j+1]);
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
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = 0.0; }
        mixin(diagonal_solve("zed", "v"));
        break;
    case PreconditionerType.jacobi:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = 0.0; }
        mixin(jacobi_solve("zed", "v"));
        break;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = 0.0; }
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
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }
        mixin(diagonal_solve("dU", "zed"));
        break;
    case PreconditionerType.jacobi:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }
        mixin(jacobi_solve("dU", "zed"));
        break;
    case PreconditionerType.sgs:
        foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }
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
 * REFERENCE: Knoll, D.A. and McHugh, P.R.,
 *            Newton-Krylov Methods Applied to a System of Convection-Diffusion-Reaction Equations,
 *            Computer Physics Communications,
 *            1994
 *
 * Authors: KAD and RJG
 * Date: 2024-02-29
 */
double computePerturbationSize()
{
    // perturbation constant, Knoll and McHugh suggest a value on the order of the square
    // root of machine roundoff, however, we have observed that this value can be somewhat
    // problem dependent.
    auto a = nkCfg.frechetDerivativePerturbation;

    // calculate L2 norm of zed vector
    number z_L2 = 0.0;
    mixin(dotOverBlocks("z_L2", "zed", "zed"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(z_L2.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    z_L2 = sqrt(z_L2);
    // When the L2 norm of the zed vector is small, the sigma calculated from this routine
    // has been observed to result in noisy convergence as observed in the global residual.
    // We enforce a minimum value to guard against this behaviour.
    // From some experimentation a minimum value of sqrt(1.0e-07) appears sufficient.
    double min_z_L2 = sqrt(1.0e-07);
    if (z_L2.re <= min_z_L2) { z_L2.re = min_z_L2; }

    // compute sum of individual component perturbations, eps_i,
    // where eps_i is calculated using Eq. 12, i.e. eps_i = a*|var_i| + a
    double eps_sum = 0.0;
    double nvars = 0.0;
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            foreach (var; cell.U[0]) {
                eps_sum += a*fabs(var.re) + a;
                nvars += 1;
            }
        }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(eps_sum.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(nvars.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    // compute perturbation parameter via Eq. 11
    double sigma = eps_sum/(nvars*z_L2.re);

    return sigma;
}

/**
 * Evaluate the augmented Jacobian-vector product w = A*z where A = I/dt - J,
 * here J is the Jacobian matrix and I/dt is the time-step augmentation term.
 *
 * Note that this routine expects the z (or zed) vector to be set prior to calling,
 * and returns the product A*z stored in the w vector
 *
 * Authors: KAD and RJG
 * Date: 2023-02-27
 */
void evalAugmentedJacobianVectorProduct()
{
    size_t nConserved = GlobalConfig.cqi.n;

    // Prepare w vector with 1/dt term: (I/dt)z
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t startIdx = 0;
        foreach (cell; blk.cells) {
            double dtInv = 1.0/cell.dt_local;
            blk.w[startIdx .. startIdx + nConserved] = dtInv*blk.zed[startIdx .. startIdx + nConserved];
            startIdx += nConserved;
        }
    }

    // Determine perturbation size, sigma
    double sigma;
    version (complex_numbers) {
        // For complex-valued Frechet derivative, a very small perturbation
        // works well (almost) all the time.
        sigma = nkCfg.frechetDerivativePerturbation;
    }
    else {
        // For real-valued Frechet derivative, we compute a perturbation size
        // based on the zed vector.
        sigma = computePerturbationSize();
    }

    // Evaluate Jz and place result in zed vector
    evalJacobianVectorProduct(sigma);

    // Complete the calculation of w: (I/dt)z - Jz
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars)  blk.w[k] = blk.w[k] - blk.zed[k];
    }
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
    if (ftl == 0 && GlobalConfig.do_shock_detect && !GlobalConfig.frozen_shock_detector) { detect_shocks(0, ftl); }

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
            blk.average_turbulent_transprops_to_faces();
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
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = activePhase.frozenLimiterForJacobian; }

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
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = activePhase.frozenLimiterForResidual; }
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
    foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = activePhase.frozenLimiterForJacobian; }

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
                blk.zed[startIdx+ivar] = (cell.dUdt[1][ivar].re - blk.R0[startIdx+ivar])/(sigma);
            }
            cell.decode_conserved(0, 0, 0.0);
            startIdx += nConserved;
        }
    }
    foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(activePhase.residualInterpolationOrder);
    foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = activePhase.frozenLimiterForResidual; }
}

/**
 * Determine a relaxation factor based on a physicality check.
 *
 * In this algorithm, the relaxation factor keeps falling as we search across cells in order.
 * This is efficient since we're searching for a worst case. If at any point, the relaxation
 * factor gets smaller than what we're prepared to accept then we just break the search
 * in that block of cells.
 *
 * REFERENCE: A. Yildirim, G. K. W. Kenway, C. A. Mader, and J. R. R. A. Martins,
 *            A Jacobian-free approximate Newton-Krylov startup strategy for RANS simulations,
 *            Journal of Computational Physics,
 *            2019
 *
 * Authors: KAD and RJG
 * Date: 2022-03-05
 */
double determineRelaxationFactor()
{
    alias GlobalConfig cfg;
    double theta = nkCfg.allowableRelativeMassChange;
    double minOmega = nkCfg.minRelaxationFactorForUpdate;
    double omegaReductionFactor = nkCfg.relaxationFactorReductionFactor;

    size_t nConserved = cfg.cqi.n;
    size_t massIdx = cfg.cqi.mass;

    double omega = 1.0;

    //----
    // 1. First determine a relaxation factor based on an allowable amount of mass change (based on Algorithm 1 from ref.)
    //----

    foreach (blk; parallel(localFluidBlocks,1)) {
        auto cqi = blk.myConfig.cqi;
        int startIdx = 0;
        blk.omegaLocal = 1.0;
        double U, dU, relDiffLimit;
        foreach (cell; blk.cells) {
            if (cqi.n_species == 1) {
                U = cell.U[0][cqi.mass].re;
                dU = blk.dU[startIdx+cqi.mass];
            }
            else {
                U = 0.0;
                dU = 0.0;
                foreach (isp; 0 .. cqi.n_species) {
                    U += cell.U[0][cqi.species+isp].re;
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
        blk.omegaLocal = omega;
        size_t startIdx = 0;
        foreach (cell; blk.cells) {
            // make a copy of the original cell flow state
            blk.fs_save.copy_values_from(*(cell.fs));
            bool failedDecode = false;
            while (blk.omegaLocal >= minOmega) {
                foreach (i; 0 .. nConserved) {
                    cell.U[1][i] = cell.U[0][i] + blk.omegaLocal * blk.dU[startIdx+i];
                }
                try {
                    cell.decode_conserved(0, 1, 0.0);
                }
                catch (FlowSolverException e) {
                    failedDecode = true;
                }
                // return cell to original flow state
                cell.fs.copy_values_from(*(blk.fs_save));

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
 * Apply a backtracking line search to determine a relaxation factor that reduces the unsteady residual.
 *
 * The line search is applied *after* the physicality check is performed, so we assume
 * that the maximum allowable step (omega*dU) recovers a physically realizable state.
 *
 * This implementation is based on Algorithm A6.3.1 from pg. 325 of Dennis and Schnabel.
 *
 * The algorithm frames the problem in terms of finding a vector x* which
 * minimizes the objective function f defined as f(x) = 1/2 F(x)^T F(x)
 * where in our fluid solver context F = R (the residual function) and
 * x = U (the conserved quantities vector). See Section 6.5 from Dennis and
 * Schnabel for more details.
 *
 * Note that this particular algorithm is based on satisfying the Armijo
 * condition of sufficient decrease and does not consider the sufficient
 * curvature condition (together these constitute the Wolfe conditions).
 * This is typically acceptable when using backtracking to determine an
 * appropriate step length, see Nocedal (pg. 37).
 *
 *
 * REFERENCES: J.E. Dennis, Jr., and R. B. Schnabel,
 *             Numerical Methods for Unconstrained Optimization and Nonlinear Equations,
 *             Classics Appl. Math. 16, SIAM, Philadelphia (1996)
 *
 *             J. Nocedal, and S. J. Wright,
 *             Numerical Optimization (2nd edition),
 *             Springer Science & Business Media (2006)
 *
 *
 * Authors: KAD and RJG
 * Date: 2024-03-04
 */
double applyLineSearch(double omega)
{
    size_t nConserved = GlobalConfig.cqi.n;
    int lineSearchOrder = nkCfg.lineSearchOrder;
    if (lineSearchOrder < 1 || lineSearchOrder > 3) {
        throw new NewtonKrylovException("Invalid line_search_order; valid options: 1 (linear), 2 (quadratic), 3 (cubic).");
    }
    double minOmega = nkCfg.minRelaxationFactorForUpdate;
    double lambdaReductionFactor = nkCfg.relaxationFactorReductionFactor;
    double alpha = 1.0e-04;                // a constant used in the step-acceptance test
    double lambda = 1.0, lambdaPrev = 1.0; // step length
    double f = 0.0;                        // objective function evaluated at starting point
    double fPlus, fPlusPrev;               // objective function evaluated at new point
    double initSlope = 0.0;                // intial slope used in the step-acceptance test


    // evaluate objective function at starting point
    mixin(dotOverBlocks("f", "R", "R"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    f = 0.5*f;

    // evaluate the initial slope := - F(x)^T J(x) p
    // where we take p = omega*dU as the search direction, see Dennis and Schnabel (pg. 148)
    //
    // 1. evaluate J(x).p via a Frechet derivative
    // 1.a. place omega*dU in zed
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) blk.zed[k] = omega*blk.dU[k];
    }
    // 1.b. evaluate Jacobian-vector product w := A . zed
    evalAugmentedJacobianVectorProduct();
    //
    // 2. evaluate slope := - F(x)^T w(x,p)
    mixin(dotOverBlocks("initSlope", "w", "R"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &initSlope, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    initSlope = -initSlope;

    bool reduceLambda = true;
    while (reduceLambda) {

	//----
	// 1. Compute unsteady term
	//----
	foreach (blk; parallel(localFluidBlocks,1)) {
	    size_t startIdx = 0;
	    foreach (cell; blk.cells) {
		blk.R[startIdx .. startIdx+nConserved] = -(1.0/cell.dt_local) * lambda * omega * blk.dU[startIdx .. startIdx+nConserved];
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
		    cell.U[1][ivar] = cell.U[0][ivar] + lambda * omega * blk.dU[startIdx+ivar];
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
                    blk.R[startIdx+ivar] += cell.dUdt[1][ivar].re;
                }
                // return cell to original state
		cell.decode_conserved(0, 0, 0.0);
		startIdx += nConserved;
	    }
	}

	//----
	// 3. Add smoothing source term
	//----
	if (nkCfg.useResidualSmoothing) {
            applyResidualSmoothing();
	}

	//----
	// 4. Compute objective function at new point
	//----
	mixin(dotOverBlocks("fPlus", "R", "R"));
	version(mpi_parallel) {
	    MPI_Allreduce(MPI_IN_PLACE, &fPlus, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
	fPlus = 0.5*fPlus;

        //debug {
        //    writefln("lambda: %.8e    RUn: %.8e    RU0: %.8e    test: %.8e", lambda, fPlus, f, f + alpha*lambda*initSlope);
        //}

	//----
	// 5. Check if objective function is reduced
	//----
	if ( (fPlus <= f + alpha*lambda*initSlope) || (lambda*omega < minOmega) ) {
	    // we are done, don't attempt to reduce lambda any further
	    reduceLambda = false;
	}
	else {
            if (lineSearchOrder == 1) {
                // just backtrack
                lambda *= lambdaReductionFactor;
            } else {
                double lambdaNext;
                if (lineSearchOrder == 2 || (lineSearchOrder == 3 && lambda == 1.0)) {
                    // use a quadratic fit
                    lambdaNext = -initSlope / (2*(fPlus-f-initSlope));
                } else {
                    // use a cubic fit
                    double v1 = fPlus-f-lambda*initSlope;
                    double v2 = fPlusPrev-f-lambdaPrev*initSlope;
                    double a  = (v1/(lambda*lambda) - v2/(lambdaPrev*lambdaPrev))/(lambda-lambdaPrev);
                    double b  = (-lambdaPrev*v1/(lambda*lambda) + lambda*v2/(lambdaPrev*lambdaPrev))/(lambda-lambdaPrev);
                    double d  = b*b - 3*a*initSlope;
                    if (d < 0.0) d = 0.0; // prevent taking sqrt of -ve number
                    if (a == 0.0) {
                        // the cubic is a quadratic
                        lambdaNext = -initSlope/(2.0*b);
                    } else {
                        // legitimate cubic
                        lambdaNext = (-b + sqrt(d))/(3.0*a);
                    }
                }
                // setup for the next step
                lambdaPrev = lambda;
                fPlusPrev = fPlus;
                // place some bounds on lambda
                if (lambdaNext > 0.5*lambda) {
                    lambdaNext = 0.5*lambda;
                }
                if (lambdaNext <= 0.1*lambda) {
                    lambda = 0.1*lambda;
                } else {
                    lambda = lambdaNext;
                }
            }
	}
    }

    return lambda*omega;
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
 * History:
 *   2024-03-26 changed format to match other output files (e.g. loads)
 */
void initialiseDiagnosticsFile()
{
    ensure_directory_is_present(diagDir);
    diagnostics = File(diagFile, "w");

    alias cfg = GlobalConfig;

    string header;
    header ~= "step wall-clock phase dt CFL relaxation-factor linear-solve-residual-rel-target linear-solve-residual-rel-achieved ";
    header ~= "n-restarts n-iters n-fn-calls mass-balance global-residual-abs global-residual-rel ";
    foreach (ivar; 0 .. cfg.cqi.n) {
        header ~= format("%s-abs %s-rel ", cfg.cqi.names[ivar], cfg.cqi.names[ivar]);
    }
    diagnostics.writeln(header);
    diagnostics.close();
}


/**
 * Update diagnostics file with current status.
 *
 * Authors: KAD and RJG
 * Date: 2022-07-24
 * History:
 *   2024-03-26 added some additional entries and reordered contents
 */
void writeDiagnostics(int step, double dt, double cfl, double wallClockElapsed, double omega, int phase, ref bool residualsUpToDate)
{
    alias cfg = GlobalConfig;

    double massBalance = compute_mass_balance();
    if (!residualsUpToDate) {
        computeResiduals(currentResiduals);
        residualsUpToDate = true;
    }

    // We don't need to proceed on ranks other than master.
    if (!cfg.is_master_task) return;

    diagnostics = File(diagFile, "a");
    diagnostics.writef("%8d %.8f %2d %20.16e %20.16e %20.16e %20.16e %20.16e %2d %4d %8d %20.16e %20.16e %20.16e ",
                       step, wallClockElapsed, phase, dt, cfl, omega, activePhase.linearSolveTolerance, gmresInfo.finalResidual/gmresInfo.initResidual,
                       gmresInfo.nRestarts, nkCfg.maxLinearSolverIterations*gmresInfo.nRestarts+gmresInfo.iterationCount, fnCount,
                       massBalance, globalResidual, globalResidual/referenceGlobalResidual);
    foreach (ivar; 0 .. cfg.cqi.n) {
        diagnostics.writef("%20.16e %20.16e ", currentResiduals[ivar].re, currentResiduals[ivar].re/referenceResiduals[ivar].re);
    }
    diagnostics.writef("\n");
    diagnostics.close();
}

void writeSnapshot(int step, double dt, double cfl, int currentPhase, int stepsIntoCurrentPhase, ref int nWrittenSnapshots)
{
    alias cfg = GlobalConfig;
    size_t nConserved = cfg.cqi.n;
    if (cfg.is_master_task) {
        writefln("+++++++++++++++++++++++++++++++++++++++");
        writefln("+   Writing snapshot at step = %6d +", step);
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

        writeFluidValues(nWrittenSnapshots);
        if (nkCfg.writeLimiterValues) {
            writeLimiterValues(nWrittenSnapshots);
        }
        if (nkCfg.writeResidualValues) {
            writeResidualValues(nWrittenSnapshots);
        }
        if (nkCfg.writeGradientValues) {
            writeGradientValues(nWrittenSnapshots);
        }

	// Add restart info
	snapshots ~= RestartInfo(nConserved);
	snapshots[$-1].step = step;
	snapshots[$-1].dt = dt;
	snapshots[$-1].cfl = cfl;
	snapshots[$-1].phase = currentPhase;
	snapshots[$-1].stepsIntoPhase = stepsIntoCurrentPhase;
	snapshots[$-1].globalResidual = globalResidual;
	snapshots[$-1].prevGlobalResidual = prevGlobalResidual;
	snapshots[$-1].residuals = currentResiduals.dup;
    }
    else {
        // We need to shuffle all of the snapshots
        foreach (iSnap; 2 .. nkCfg.totalSnapshots+1) {
            foreach (blk; localFluidBlocks) {
                auto fromName = fluidFilename(iSnap, blk.id);
                auto toName = fluidFilename(iSnap-1, blk.id);
                rename(fromName, toName);
                if (nkCfg.writeLimiterValues) {
                    fromName = limiterFilename(iSnap, blk.id);
                    toName = limiterFilename(iSnap-1, blk.id);
                    rename(fromName, toName);
                }
                if (nkCfg.writeResidualValues) {
                    fromName = residualFilename(iSnap, blk.id);
                    toName = residualFilename(iSnap-1, blk.id);
                    rename(fromName, toName);
                }
                if (nkCfg.writeGradientValues) {
                    fromName = gradientFilename(iSnap, blk.id);
                    toName = gradientFilename(iSnap-1, blk.id);
                    rename(fromName, toName);
                }
            }
        }

        writeFluidValues(nkCfg.totalSnapshots);
        if (nkCfg.writeLimiterValues) {
            writeLimiterValues(nkCfg.totalSnapshots);
        }
        if (nkCfg.writeResidualValues) {
            writeResidualValues(nkCfg.totalSnapshots);
        }
        if (nkCfg.writeGradientValues) {
            writeGradientValues(nkCfg.totalSnapshots);
        }

	// Shuffle the restart info
	snapshots[0 .. nkCfg.totalSnapshots-1] = snapshots[1 .. nkCfg.totalSnapshots];
	// Add new info at end
	snapshots[$-1].step = step;
	snapshots[$-1].dt = dt;
	snapshots[$-1].cfl = cfl;
	snapshots[$-1].phase = currentPhase;
	snapshots[$-1].stepsIntoPhase = stepsIntoCurrentPhase;
	snapshots[$-1].globalResidual = globalResidual;
	snapshots[$-1].prevGlobalResidual = prevGlobalResidual;
	snapshots[$-1].residuals = currentResiduals.dup;
    }

    if (GlobalConfig.is_master_task) {
	writeRestartMetadata(snapshots);
        writeln(); // write a blank line for aesthetic purposes
    }
}

void writeRestartMetadata(RestartInfo[] snapshots)
{
    auto f = File(lmrCfg.restartFile, "w");
    auto timeNow =  cast(DateTime)(Clock.currTime());
    f.writefln("# Restart metadata written at: %s", timeNow.toSimpleString());
    f.writefln("# step,  dt,  cfl, phase, steps-taken-into-phase, global-residual, pre-global-residual, n-conserved quantities residuals");
    foreach (snap; snapshots) {
	f.writef("%04d %.18e %.18e %04d %04d %.18e %.18e", snap.step, snap.dt, snap.cfl, snap.phase, snap.stepsIntoPhase, snap.globalResidual, snap.prevGlobalResidual);
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
            snapshots[$-1].phase = to!int(tks[3]);
            snapshots[$-1].stepsIntoPhase = to!int(tks[4]);
	    snapshots[$-1].globalResidual = to!double(tks[5]);
	    snapshots[$-1].prevGlobalResidual = to!double(tks[6]);
	    size_t start_idx = 7;
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
    formattedWrite(writer, "%15s\t\tRELATIVE\t\tABSOLUTE\n", "RESIDUALS");
    formattedWrite(writer, "%15s\t\t%10.6e\t\t%10.6e\n", "global", globalResidual.re/referenceGlobalResidual.re, globalResidual.re);
    foreach (ivar; 0 .. cqi.n) {
        formattedWrite(writer, "%15s\t\t%10.6e\t\t%10.6e\n", cqi.names[ivar], currentResiduals[ivar].re/referenceResiduals[ivar].re, currentResiduals[ivar].re);
    }
    writeln(writer.data);
}

/**
 * Write fluid field values to disk.
 *
 * Authors: RJG and KAD
 * Date: 2024-03-08
 */
void writeFluidValues(int iSnap)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("    |");
        writefln("    |-->  Writing fluid values");
    }

    foreach (blk; localFluidBlocks) {
        auto fileName = fluidFilename(iSnap, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        fluidBlkIO.writeVariablesToFile(fileName, cells);
    }
}

/**
 * Write limiter values to disk.
 *
 * Authors: RJG and KAD
 * Date: 2023-08-13
 */
void writeLimiterValues(int iSnap)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("    |");
        writefln("    |-->  Writing limiter values");
    }

    foreach (blk; localFluidBlocks) {
        auto fileName = limiterFilename(iSnap, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        limBlkIO.writeVariablesToFile(fileName, cells);
    }
}

/**
 * Write residual values to disk.
 *
 * Authors: KAD and RJG
 * Date: 2024-03-07
 */
void writeResidualValues(int iSnap)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("    |");
        writefln("    |-->  Writing residual values");
    }

    foreach (blk; localFluidBlocks) {
        auto fileName = residualFilename(iSnap, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        residBlkIO.writeVariablesToFile(fileName, cells);
    }
}

/**
 * Write gradient values to disk.
 *
 * Authors: RJG and KAD
 * Date: 2024-08-01
 */
void writeGradientValues(int iSnap)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("    |");
        writefln("    |-->  Writing gradient values");
    }

    foreach (blk; localFluidBlocks) {
        auto fileName = gradientFilename(iSnap, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        gradBlkIO.writeVariablesToFile(fileName, cells);
    }
}

/**
 * Write loads to disk.
 *
 * Authors: RJG and KAD
 * Date: 2023-11-19
 */
void writeLoads(int step, ref int nWrittenLoads)
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("+++++++++++++++++++++++++++++++++++++++");
        writefln("+   Writing loads at step = %6d    +", step);
        writefln("+++++++++++++++++++++++++++++++++++++++");
        writeln();
    }
    if (cfg.is_master_task) {
        init_current_loads_indx_dir(nWrittenLoads);
    }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    wait_for_current_indx_dir(nWrittenLoads);
    writeLoadsToFile(-1, nWrittenLoads);
    if (cfg.is_master_task) {
        update_loads_metadata_file(-1, step, nWrittenLoads);
    }

    nWrittenLoads++;

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
               blk.rhs[] = 0.0;
               nm.smla.multiply_block_upper_triangular(blk.flowJacobian.local, blk."~lhs_vec~"[], blk.rhs[], blk.cells.length, nConserved);
               nm.smla.multiply_block_lower_triangular(blk.flowJacobian.local, blk."~lhs_vec~"[], blk.rhs[], blk.cells.length, nConserved);
               blk.rhs[] = blk."~rhs_vec~"[] - blk.rhs[];
               blk."~lhs_vec~"[] = 0.0;
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
             blk.rhs[] = 0.0;
             nm.smla.multiply_block_upper_triangular(blk.flowJacobian.local, blk."~lhs_vec~", blk.rhs, blk.cells.length, nConserved);
             blk.rhs[] = blk."~rhs_vec~"[] - blk.rhs[];
             nm.smla.block_lower_triangular_solve(blk.flowJacobian.local, blk.rhs[], blk."~lhs_vec~"[], blk.cells.length, nConserved);

             // backward sweep
             blk.rhs[] = 0.0;
             nm.smla.multiply_block_lower_triangular(blk.flowJacobian.local, blk."~lhs_vec~"[], blk.rhs[], blk.cells.length, nConserved);
             blk.rhs[] = blk."~rhs_vec~"[] - blk.rhs[];
             nm.smla.block_upper_triangular_solve(blk.flowJacobian.local, blk.rhs, blk."~lhs_vec~"[], blk.cells.length, nConserved);
         }
    }

    }";
    return code;
}

