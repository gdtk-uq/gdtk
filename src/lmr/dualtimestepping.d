/**
 * Core set of functions used to advance the flow field in time using a dual time stepping scheme.
 *
 * Authors: RJG and KAD
 * Date: 2025-009-06
 * History:
 *   2025-09-06 First pass implementation
 */

module lmr.dualtimestepping;

import core.memory : GC;
import core.stdc.stdlib : exit;
import std.algorithm : min;
import std.algorithm.searching : countUntil;
import std.array : appender;
import std.conv : to;
import std.datetime : DateTime, Clock;
import std.file;
import std.format : formattedWrite;
import std.json : JSONValue;
import std.math;
import std.parallelism : parallel, defaultPoolThreads;
import std.range : walkLength;
import std.stdio : File, write, writeln, writefln, stdout;
import std.string;
import std.typecons : Tuple, tuple;

import dyaml;

import geom;
import nm.bbla;
import nm.number : number;
import nm.smla;
import ntypes.complex;
import util.json_helper;
import util.lua;
import util.lua_service;
import util.time_utils : timeStringToSeconds;

import lmr.newtonkrylovsolver;
import lmr.history : initHistoryCells, writeHistoryCellsToFiles;
import lmr.timemarching : addToTimesFile, prepareTimesFileOnRestart, nextXxxxTime, writeSnapshotFiles_timemarching, writeLoadsFiles_timemarching;
import lmr.bc;
import lmr.blockio;
import lmr.conservedquantities : ConservedQuantities, copy_values_from;
import lmr.fileutil : ensure_directory_is_present;
import lmr.fluidblock : FluidBlock;
import lmr.fvcell : FVCell;
import lmr.fvcellio;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.init;
import lmr.lmrconfig;
import lmr.lmrexceptions;
import lmr.lmrerrors : LmrError, lmrErrorExit;
import lmr.lmrwarnings;
import lmr.lua_helper;
import lmr.sfluidblock : SFluidBlock;
import lmr.simcore : compute_mass_balance, call_UDF_at_write_to_file;
import lmr.simcore_exchange;
import lmr.simcore_gasdynamic_step : detect_shocks;
import lmr.ufluidblock : UFluidBlock;
import lmr.user_defined_source_terms : getUDFSourceTermsForCell;
import lmr.loads : writeLoadsToFile,
                   init_current_loads_indx_dir,
                   wait_for_current_indx_dir,
                   writeLoadsToFile,
                   count_written_loads,
                   update_loads_metadata_file;
import lmr.grid_motion;
import lmr.grid_motion_udf;
import lmr.grid_motion_shock_fitting;
import lmr.lmrwarnings;

version(mpi_parallel) {
    import mpi;
}

/**
 * Time step controller for the physical time step based on the procedure from Howard (2010).
 *
 * REFERENCES:   M. A. Howard
 *               PhD thesis,
 *               Univesity of Colorado (2010)
 *
 *               P. M. Gresho, R. L. Lee, R. L. Sani, and T. W. Stullich
 *               On the Time-Dependent FEM Solution of the Incompressible Navier-Stokes Equations in Two- and Three-Dimensions
 *               International Conference on Numerical Methods in Laminar and Turbulent Flow (1978)
 *
 *
 * Authors: RJG and KAD
 * Date: 2025-09-07
 */

struct DualTimeStepController {

    double nextTimeStep()
    {
        size_t nConserved = GlobalConfig.cqi.n;

        // model parameters used in eqn. 2.112
        double a,b;
        if (bdfOrder == 2) {
            a = 1.0/3.0;
            b = 3.0*(1.0+dtPrev/dt);
        } else { // we assume BDF1
            a = 0.5;
            b = 2.0;
        }

        // save the previous physical time-step
        dtPrev = dt;

        // Reference conserved quantity magnitude from Howard (2010) eqn 2.113
        double U_inf_mag = maxConservedQuantityVectorMagnitude();

        // L2 norm of conserved quantity update from Howard (2010) eqn 2.113
        double dU_L2 = conservedQuantityUpdateNorm();

        // this is an estimate of the (relative) norm of the local truncation error for the step just completed from Howard (2010) eqn. 2.113
        double d = (dU_L2/(sqrt(dof)*U_inf_mag));

        // This input value is the maximum permissible value of the relative error in a single step,
        // we set this value to 1.0e-03 as recommended in Gresho et al. (1978)
        double eps = 1.0e-03;

        // time-step update from Howard (2010) eqn. 2.112
        dt = dtPrev*pow((b*eps/d),a);

        // Howard mentions that the time-step is typically limited to grow by no more than 10% to 20%
        dt = fmin(1.1*dtPrev, dt);
        // we will also limit the time-step to some maximum allowable value
        dt = fmin(dt, dtMax);

        return dt;
    }

    @nogc double evalDtInvEff()
    {
        double dtInv;
        if (bdfOrder == 2) {
            dtInv = 1.0/dt + 1.0/dtPrev - dt/(dtPrev*(dt+dtPrev));
        } else { // we assume BDF1
            dtInv = 1.0/dt;
        }
        return dtInv;
    }

private:
    double dt;
    double dtPrev;
    double dtMax;
    int bdfOrder;
    int targetBdfOrder;
    int bdfStartupSteps = 10;
    double dof = 0; // degrees of freedom from Howard (2010) eqn 2.113
}

DualTimeStepController dtsController;

void initDualTimeNewtonKrylovSimulation(int snapshotStart, int maxCPUs, int threadsPerMPITask, string maxWallClock)
{
    alias cfg = GlobalConfig;
    nkCfg.dualTime = true;

    if (cfg.verbosity_level > 0 && cfg.is_master_task) {
        writeln("lmr run: Begin initDualTimeNewtonKrylovSimulation()...");
    }

    // Start timer
    SimState.maxWallClockSeconds = timeStringToSeconds(maxWallClock);
    SimState.wall_clock_start = Clock.currTime();

    // Initialise baseline configuration
    JSONValue cfgData= initConfiguration();

    // After loading configuration, we can check if there are any settings
    // that are not particularly appropriate for the steady-state solver.
    if (cfg.extrema_clipping) {
        writeln(warningMessage(LmrWarning.useOfExtremaClipping));
    }
    if (cfg.nFluidBlocks == 0 && cfg.is_master_task) {
        throw new NewtonKrylovException("No FluidBlocks; no point in continuing with simulation initialisation.");
    }

    // if we have grid motion, we need two grid time levels
    if (cfg.grid_motion != GridMotion.none) {
        cfg.n_grid_time_levels = 2;
    }

    initLocalBlocks();

    /* RJG, 2024-03-05
     * The shared memory model is not playing nicely with how the Newton-Krylov module
     * is set up. For the present, we limit the threads to 1.
     */
    initThreadPool(1, 1);

    initFluidBlocksBasic(cfgData);
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

    initFullFaceDataExchange(cfgData);
    initMappedCellDataExchange();
    initGhostCellGeometry();
    initLeastSquaresStencils();
    initStructuredStencilData();

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

    if (cfg.grid_motion == cfg.grid_motion.shock_fitting) {
        initShockFitting(cfgData);
    }

    orderBlocksBySize();
    initMasterLuaState();
    initCornerCoordinates();
    if (cfg.turb_model.needs_dwall) initWallDistances();

    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

    // [TODO] Add in electric field solver initialisation.

    // Read the control because it sets up part of initial configuration
    readControl();

    // Create progress file and add 0th step entry.
    initTimeMarchingProgressFile();

    SimState.is_restart = (snapshotStart > 0);
    SimState.current_tindx = snapshotStart;
    if (SimState.current_tindx == 0) {
        // On a fresh start, set loads tindx to 0
        SimState.current_loads_tindx = 0;
    } else {
        // We'll need to work a little harder and find out how many loads we've already written.
        SimState.current_loads_tindx = count_written_loads();
    }

    initHistoryCells();
    if (cfg.write_loads && (SimState.current_loads_tindx == 0)) {
        if (cfg.is_master_task) {
            initLoadsFiles();
        }
    }

    // Flags to indicate that the saved output is fresh.
    // On startup or restart, it is assumed to be so.
    SimState.output_just_written = true;
    SimState.history_just_written = true;
    SimState.loads_just_written = true;

    // When starting a new calculation,
    // set the global time step to the initial value.
    // If we have elected to run with a variable time step,
    // this value may get revised on the very first step.
    SimState.dt_global = cfg.dt_init;
    dtsController.dt = cfg.dt_init;
    dtsController.dtPrev = cfg.dt_init;
    dtsController.dtMax = cfg.dt_max;
    dtsController.dof = computeDoFs();
    SimState.cfl_max = 0.0;
    // BDF schemes are not self-starting, we will take some low order start up steps
    if (cfg.dualtimestepping_update_scheme == DualTimeSteppingUpdate.bdf2) {
        dtsController.bdfOrder = 1;
        dtsController.targetBdfOrder = 2;
    } else {
        dtsController.bdfOrder = 1;
        dtsController.targetBdfOrder = 1;
    }

    // dt_global will be updated if restart branch is selected below.
    if (!SimState.is_restart) {
        SimState.time = 0.0;
        SimState.step = 0;
        if (cfg.is_master_task) {
            // Clean out any existing times file.
            if (lmrCfg.timesFile.exists) lmrCfg.timesFile.remove;
            // We can put an initial entry in the times file now.
            addToTimesFile();
        }
    } else {
        // SimState.time, SimState.step and SimState.dt_global set in prepareTimesFileOnRestart
        // It's convenient to do that while we have the times file loaded.
        if (cfg.is_master_task) {
            // RJG, 2024-11-06
            // For large-scale work (many processes on parallel filesystems),
            // we need to be careful about attempting to read one small file from many processes.
            // The solution is to read only from master, then broadcast to other ranks.
            // This is the same solution NNG implemented in the steady-state codepath for
            // broadcasting restart information.
            //
            // This issue was noted by Sebastiaan van Oeveren when working on the UQ Bunya cluster.
            prepareTimesFileOnRestart(SimState.current_tindx);
            writeln("*** RESTARTING SIMULATION ***");
            writefln("RESTART-SNAPSHOT: %d", SimState.current_tindx);
            writefln("RESTART-STEP: %d", SimState.step);
            writefln("RESTART-TIME: %8.3e", SimState.time);
        }
        version(mpi_parallel) {
            // After setting SimState on master, broadcast to other ranks
            // 0. Make non-shared temporaries
            double simTime = SimState.time;
            int step = SimState.step;
            double dtGlobal = SimState.dt_global;
            // 1. Broadcast
            MPI_Bcast(&simTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&dtGlobal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // 2. Write to shared variables
            SimState.time = simTime;
            SimState.step = step;
            SimState.dt_global = dtGlobal;
        }
    }

    // Do some memory clean-up and reporting.
    GC.collect();
    GC.minimize();
    if (cfg.verbosity_level > 0 && cfg.is_master_task) {
        // For reporting wall-clock time, convert to seconds with precision of milliseconds.
        double wall_clock_elapsed = to!double((Clock.currTime() - SimState.wall_clock_start).total!"msecs"())/1000.0;
        writefln("lmr run: Done initDualTimeNewtonKrylovSimulation() at wall-clock(WC)= %.1f sec", wall_clock_elapsed);
        stdout.flush();
    }
}

void performDualTimeNewtonKrylovUpdates(int snapshotStart, double startCFL, int maxCPUs, int threadsPerMPITask)
{
    alias cfg = GlobalConfig;

    //----------------------------------------------
    // Some Newton-Krylov specific initialisations
    //----------------------------------------------
    readNewtonKrylovConfig();

    // Safety check: ensure grid motion is disabled before proceeding
    if (cfg.grid_motion != GridMotion.none && cfg.is_master_task) {
        throw new NewtonKrylovException("Grid motion is not compatible with dual time stepping.");
    }

    // safety check: we anticipate each nonlinear solve only requiring very few steps to converge
    if (nkCfg.numberOfStepsForSettingReferenceResiduals != 0) {
        throw new NewtonKrylovException("Reference residuals should be evaluated at step 0 with dual time stepping.");
    }

    int cumulativeNewtonSteps = 0;
    int startStep = 1;
    bool finalStep = false;
    double cfl;
    double dt;
    int currentPhase = 0;
    int stepsIntoCurrentPhase = 0;
    bool updatePreconditionerThisStep = false;
    CFLSelector cflSelector;
    int numberBadSteps = 0;
    bool startOfNewPhase = false;
    double omega = 1.0;
    bool terminalPhase = false;

    if (nkCfg.usePreconditioner) initPreconditioner();
    size_t nConserved = cfg.cqi.n;
    referenceResiduals = new ConservedQuantities(nConserved);
    currentResiduals = new ConservedQuantities(nConserved);
    rowScale = new ScaleFactors(nConserved, cfg.dimensions);
    invColScale = new ScaleFactors(nConserved, cfg.dimensions);
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
    if (nkCfg.useFGMRES) { allocateGlobalFGMRESWorkspace(); }

    foreach (blk; localFluidBlocks) {
        if (nkCfg.useFGMRES) {
            blk.allocate_FGMRES_workspace(nkCfg.maxLinearSolverIterations,
                                          nkCfg.maxFgmresPreconditioningIterations,
                                          nkCfg.useRealValuedFrechetDerivative);
        } else {
            blk.allocate_GMRES_workspace(nkCfg.maxLinearSolverIterations,
                                         nkCfg.useRealValuedFrechetDerivative);
        }
    }

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

    // We need to create the diagnostics file
    // note: upon a restart we write a new diagnostics file --> old file will be lost
    // TODO: make a copy of old diagnotistics files
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

    //----------------------------------------------
    // Some time marching specific initialisations
    //----------------------------------------------
    // The next time for output...
    SimState.t_plot = nextXxxxTime(SimState.time, cfg.dt_plot_schedule);
    SimState.t_history = nextXxxxTime(SimState.time, cfg.dt_history);
    SimState.t_loads = nextXxxxTime(SimState.time, cfg.dt_loads);

    // fill conserved quantity arrays with valid data
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (cell; blk.cells) {
            cell.encode_conserved(0, 0, 0.0);
            foreach (n; 0 .. cfg.n_flow_time_levels) {
                cell.U[n].copy_values_from(cell.U[0]);
            }
        }
    }

    //----------------------------------------------------------------
    // Top of main time-stepping loop
    //----------------------------------------------------------------
    double wall_clock_elapsed;
    SimState.wall_clock_start = Clock.currTime();
    SimState.target_time = cfg.max_time;

    // Normally, we can terminate upon either reaching a maximum time or upon reaching a maximum iteration count.
    shared bool finished_time_stepping = (SimState.time >= SimState.target_time) || (SimState.step >= cfg.max_step);
    while ( !finished_time_stepping ) {

        // we need to evaluate time-dependent boundary conditions and source terms at the future state
        SimState.step = SimState.step + 1;
        SimState.time = SimState.time + dtsController.dt;
        if (SimState.step >= dtsController.bdfStartupSteps) { dtsController.bdfOrder = dtsController.targetBdfOrder; }

        // reset Newton-Krylov solver parameters
        finalStep = false;
        currentPhase = 0;
        setPhaseSettings(0);
        cfl = activePhase.startCFL;
        stepsIntoCurrentPhase = 0;
        updatePreconditionerThisStep = true;
        terminalPhase = false;
        if (currentPhase == nkCfg.numberOfPhases-1) terminalPhase = true; // we are beginning in the terminal phase

        // evaluate a reference residual for the nonlinear solve
        evalResidual(0);
        setResiduals();
        addUnsteadyTermToResiduals();
        computeGlobalResidual();
        referenceGlobalResidual = globalResidual;
        computeResiduals(referenceResiduals);
        // Add value of 1.0 to each residaul.
        // If values are very large, 1.0 makes no difference.
        // If values are zero, the 1.0 should mean the reference residual
        // asymptotes to an absolute residual.
        foreach (ref residual; referenceResiduals) residual += to!number(1.0);
        if (nkCfg.numberOfStepsForSettingReferenceResiduals == 0) {
            if (cfg.is_master_task) {
                writefln("***  Reference global residual: %.12e", referenceGlobalResidual);
            }
        }

        // solve nonlinear system using Newton-Krylov solver with pseudo time stepping
        foreach (step; 1 .. nkCfg.maxNewtonSteps+1) {
            cumulativeNewtonSteps += 1;
            //----
            // 0. Check for any special actions based on step to perform at START of step
            //----
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
            evalDtInvEff();

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

            //----
            // 1. Perforn Newton update
            //----
            prevGlobalResidual = globalResidual;
            if (nkCfg.useFGMRES) {
                solveNewtonStepFGMRES();
            } else {
                solveNewtonStep(updatePreconditionerThisStep);
            }

            // 1a. perform a physicality check if required
            omega = nkCfg.usePhysicalityCheck ? determineRelaxationFactor() : 1.0;

            // 1b. do a line search if required
            if ( (omega > nkCfg.minRelaxationFactorForUpdate) && nkCfg.useLineSearch ) {
                omega = applyLineSearch(omega);
            }

            // 1c. check if we achived the allowable linear solver tolerance
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
                    if (cfg.is_master_task) {
                        string exitMsg = "Update failure in Newton step.\n";
                        exitMsg ~= format("step= %d, CFL= %e, dt= %e, global-residual= %e \n", step, cfl, dt, globalResidual);
                        exitMsg ~= "Error message from failed update:\n";
                        exitMsg ~= format("%s", e.msg);
                        exitMsg ~= "You might be able to try a smaller CFL.\n";
                        exitMsg ~= "Bailing out!\n";
                        lmrErrorExit(LmrError.unrecoverableUpdate, exitMsg);
                    }
                }
                numberBadSteps = 0;
            }
            else if ((omega < nkCfg.minRelaxationFactorForUpdate || failedToAchieveAllowableLinearSolverTolerance) && activePhase.useAutoCFL) {
                numberBadSteps++;
                if (numberBadSteps == nkCfg.maxConsecutiveBadSteps) {
                    string exitMsg = "Too many consecutive bad steps while trying to update flow state.\n";
                    exitMsg ~= format("Number of bad steps = %d\n", numberBadSteps);
                    exitMsg ~= format("Last attempted CFL = %e\n", cfl);
                    exitMsg ~= "Bailing out!\n";
                    lmrErrorExit(LmrError.unrecoverableUpdate, exitMsg);
                }
                // Return flow states to their original state for next attempt.
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (cell; blk.cells) {
                        cell.decode_conserved(0, 0, 0.0);
                    }
                }
            }
            else {
                if (cfg.is_master_task) {
                    string exitMsg = "ERROR: relaxation factor for Newton update is too small.\n";
                    exitMsg ~= format("step= %d, relaxation factor= %f\n", step, omega);
                    exitMsg ~= "Bailing out!\n";
                    lmrErrorExit(LmrError.hitNumericalRailGuard, exitMsg);
                }
            }

            //----
            // 2. Post-update actions
            //----
            // Here we need to do some house-keeping and see if we continue with iterations.

            // 2a. Reporting (to files and screen)
            if (((step % nkCfg.stepsBetweenDiagnostics) == 0) || (finalStep && nkCfg.writeDiagnosticsOnLastStep)) {
                writeDiagnostics(cumulativeNewtonSteps, dt, cfl, wall_clock_elapsed, omega, currentPhase, residualsUpToDate);
            }
            version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
            // Reporting to screen on progress.
            if (((step % nkCfg.stepsBetweenStatus) == 0) || finalStep) {
                printNewtonStepStatusToScreen(step, cfl, dt, residualsUpToDate);
            }
            // 2b. Stopping checks.
            string reasonForStop;
            wall_clock_elapsed = 1.0e-3*(Clock.currTime() - SimState.wall_clock_start).total!"msecs"();
            if (step == nkCfg.maxNewtonSteps) {
                finalStep = true;
                if (cfg.is_master_task) {
                    writeln("*** INNER SOLVER COMPLETE: Reached maximum number of steps.");
                    reasonForStop = "STOP-REASON: maximum-steps";
                }
            }
            if (!activePhase.ignoreStoppingCriteria) {
                if (globalResidual <= nkCfg.stopOnAbsoluteResidual) {
                    finalStep = true;
                    if (cfg.is_master_task) {
                        writefln("*** INNER SOLVE COMPLETE: The absolute global residual is below target value of %.6e", nkCfg.stopOnAbsoluteResidual);
                        reasonForStop = "STOP-REASON: absolute-global-residual-target";
                    }
                }
                if ((globalResidual/referenceGlobalResidual) <= nkCfg.stopOnRelativeResidual) {
                    finalStep = true;
                    if (cfg.is_master_task) {
                        writefln("*** INNER SOLVE COMPLETE: The relative global residual is below target value of %.6e", nkCfg.stopOnRelativeResidual);
                        reasonForStop = "STOP-REASON: relative-global-residual-target";
                    }
                }
            }
            if (finalStep) {
                break;
            }

        }

        // Update the time record and (occasionally) print status.
        if (cfg.is_master_task) {
            try {
                std.file.write(lmrCfg.progFile, format("%d\n", SimState.step));
            } catch (Exception e) {
                lmrErrorExit(LmrError.inputOutput, "Couldn't write to file: " ~ lmrCfg.progFile);
            }
        }
        //
        SimState.output_just_written = false;
        SimState.history_just_written = false;
        SimState.loads_just_written = false;
        if ((SimState.step % cfg.print_count) == 0) {
            // Print the current time-stepping status.
            auto writer = appender!string();
            formattedWrite(writer, "Step=%7d t=%10.3e dt=%10.3e cfl=%.2f ",
                           SimState.step, SimState.time, SimState.dt_global, SimState.cfl_max);

            // For reporting wall-clock time, convert to seconds with precision of milliseconds.
            wall_clock_elapsed = to!double((Clock.currTime()-SimState.wall_clock_start).total!"msecs"())/1000.0;
            double wall_clock_per_step = wall_clock_elapsed / SimState.step;
            double WCtFT = (cfg.max_time - SimState.time) / SimState.dt_global * wall_clock_per_step;
            double WCtMS = (cfg.max_step - SimState.step) * wall_clock_per_step;
            formattedWrite(writer, "WC=%.1f WCtFT=%.1f WCtMS=%.1f \n",
                           wall_clock_elapsed, WCtFT, WCtMS);
            if (cfg.verbosity_level >= 0 && cfg.is_master_task) {
                writeln(writer.data);
                stdout.flush();
            }
            version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

        }

        //----
        // 3. Post-nonlinear solve actions
        //----

        // 3a. (Occasionally) Write out an intermediate solution
        if ( SimState.step == cfg.write_flow_solution_at_step ) {
            writeSnapshotFiles_timemarching();
            GC.collect();
            GC.minimize();
        }
        if ((SimState.time >= SimState.t_plot) && !SimState.output_just_written) {
            writeSnapshotFiles_timemarching();
            if (cfg.udf_supervisor_file.length > 0) { call_UDF_at_write_to_file(); }
            SimState.output_just_written = true;
            SimState.t_plot = nextXxxxTime(SimState.time, cfg.dt_plot_schedule);
            GC.collect();
            GC.minimize();
        }

        // 3b. (Occasionally) Write out the cell history data and loads on boundary groups data
        if ((SimState.time >= SimState.t_history) && !SimState.history_just_written) {
            writeHistoryCellsToFiles(SimState.time);
            SimState.history_just_written = true;
            SimState.t_history = nextXxxxTime(SimState.time, cfg.dt_history);
            GC.collect();
            GC.minimize();
        }
        if (cfg.write_loads &&
            ( ((SimState.time >= SimState.t_loads) && !SimState.loads_just_written) ||
              SimState.step == cfg.write_loads_at_step )) {
            writeLoadsFiles_timemarching();
            SimState.loads_just_written = true;
            SimState.current_loads_tindx = SimState.current_loads_tindx + 1;
            SimState.t_loads = nextXxxxTime(SimState.time, cfg.dt_loads);
            GC.collect();
            GC.minimize();
        }

        // check if we should stop the simulation
        wall_clock_elapsed = (Clock.currTime() - SimState.wall_clock_start).total!"seconds"();
        if (SimState.time >= SimState.target_time) { finished_time_stepping = true; }
        if (SimState.step >= cfg.max_step) { finished_time_stepping = true; }
        if (cfg.halt_now == 1) { finished_time_stepping = true; }
        if (SimState.maxWallClockSeconds > 0 && (wall_clock_elapsed > SimState.maxWallClockSeconds)) {
            finished_time_stepping = true;
        }
        version(mpi_parallel) {
            // If one task is finished time-stepping, all tasks have to finish.
            int localFinishFlag = to!int(finished_time_stepping);
            MPI_Allreduce(MPI_IN_PLACE, &localFinishFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            finished_time_stepping = to!bool(localFinishFlag);
        }

        if(finished_time_stepping && cfg.is_master_task) {
            // Make an announcement about why we are finishing time-stepping.
            write("STOP-REASON: ");
            if (SimState.time >= SimState.target_time) {
                writefln("maximum-time");
                writefln("Reached target simulation time of %g seconds.", SimState.target_time);
            }
            if (SimState.step >= cfg.max_step) {
                writefln("maximum-steps");
                writefln("Reached maximum number of steps with step=%d.", SimState.step);
            }
            if (cfg.halt_now == 1) { writeln("Halt set in control file."); }
            if (SimState.maxWallClockSeconds > 0 && (wall_clock_elapsed > SimState.maxWallClockSeconds)) {
                writefln("maximum-wall-clock");
                writefln("Reached maximum wall-clock time with elapsed time %s.", to!string(wall_clock_elapsed));
            }
            writefln("FINAL-STEP: %d", SimState.step);
            writefln("FINAL-TIME: %g", SimState.time);
            stdout.flush();
        }

        //----
        // 4. Prepare for next physical time step
        //----

        // evaluate the next (hopefully larger) physical time step and evaluate the equivalent CFL
        if (!GlobalConfig.fixed_time_step && SimState.step >= dtsController.bdfStartupSteps) {
            SimState.dt_global = dtsController.nextTimeStep();
            dt = setDtInCells(1.0, false);
            SimState.cfl_max = SimState.dt_global/dt;
        }

        // shuffle conserved quantities
        // --> note that U[1] is used by the Newton-Krylov solver,
        //     so previous time-step solutions are stored in indexes >1
        if (dtsController.targetBdfOrder == 2) {
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (cell; blk.cells) {
                    // U_n => U_n-1
                    cell.U[3].copy_values_from(cell.U[2]);
                    // U_n+1 => U_n
                    cell.U[2].copy_values_from(cell.U[0]);
                }
            }
        } else { // we assume BDF1
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (cell; blk.cells) {
                    // U_n+1 => U_n
                    cell.U[2].copy_values_from(cell.U[0]);
                }
            }
        }

    }
}

/**
 * Add the unsteady term arising in the dual time stepping scheme to the residual vector.
 *
 * Note that U[1] is used by the Newton-Krylov solver,
 * so previous time-step solutions are stored in indexes >1
 *
 * Authors: RJG and KAD
 * Date: 2025-09-07
 */
void addUnsteadyTermToResiduals()
{
    size_t nConserved = GlobalConfig.cqi.n;
    int bdfOrder = dtsController.bdfOrder;
    double dt = dtsController.dt;
    double dtPrev = dtsController.dtPrev;
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t startIdx = 0;
        foreach (i, cell; blk.cells) {
            foreach (ivar; 0 .. nConserved) {
                if (bdfOrder == 2) {
                    // note that for a fixed time-step, the coefficients should reduce down to the standard BDF2 coefficients, i.e. c1 = 3/2, c2 = 2, and c3 = 1/2
                    double c1 = 1.0/dt + 1.0/dtPrev - dt/(dtPrev*(dt+dtPrev));
                    double c2 = 1.0/dt + 1.0/dtPrev;
                    double c3 = dt/(dtPrev*(dt+dtPrev));
                    blk.R[startIdx+ivar] += -c1*cell.U[0][ivar].re + c2*cell.U[2][ivar].re - c3*cell.U[3][ivar].re;
                } else { // we assume BDF1
                    blk.R[startIdx+ivar] += -(1.0/dt)*cell.U[0][ivar].re + (1.0/dt)*cell.U[2][ivar].re;
                }
            }
            startIdx += nConserved;
        }
    }
}

/**
 * Adds the effective physical dt inverse to the local psuedo dt inverse.
 *
 *
 * Authors: RJG and KAD
 * Date: 2025-09-07
 */
void evalDtInvEff()
{
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (cell; blk.cells) {
            cell.dt_inv += dtsController.evalDtInvEff();
        }
    }
}

/**
 * Computes the magnitude of all the conserved quantity vectors and returns the maximum.
 *
 *
 * Authors: RJG and KAD
 * Date: 2025-09-11
 */
double maxConservedQuantityVectorMagnitude() {
    size_t nConserved = GlobalConfig.cqi.n;
    double U_inf_mag = 0.0;
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            double U_inf_mag_tmp = 0.0;
            foreach (j; 0 .. nConserved) {
                U_inf_mag_tmp += cell.U[0][j].re*cell.U[0][j].re;
            }
            U_inf_mag_tmp = sqrt(U_inf_mag_tmp);
            U_inf_mag = fmax(U_inf_mag, U_inf_mag_tmp);
        }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(U_inf_mag), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    return U_inf_mag;
}

/**
 * Computes the norm of the conserved quantity update vector.
 *
 * Note: we don't have the relaxation factor at the point that this
 *       routine is called, otherwise we could compute it explicitly.
 *
 * Authors: RJG and KAD
 * Date: 2025-09-11
 */
double conservedQuantityUpdateNorm() {
    size_t nConserved = GlobalConfig.cqi.n;
    double dU_L2 = 0.0;
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            foreach (j; 0 .. nConserved) {
                dU_L2 += pow((cell.U[0][j].re - cell.U[2][j].re), 2.0);
            }
        }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(dU_L2), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    dU_L2 = sqrt(dU_L2);
    return dU_L2;
}

/**
 * Computes the degrees of freedom for the simulation.
 *
 *
 * Authors: RJG and KAD
 * Date: 2025-09-11
 */
double computeDoFs() {
    size_t nConserved = GlobalConfig.cqi.n;
    double dof = 0.0;
    foreach (blk; localFluidBlocks) {
        dof = blk.cells.length * nConserved;
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(dof), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    return dof;
}

void printNewtonStepStatusToScreen(int step, double cfl, double dt, ref bool residualsUpToDate)
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
    formattedWrite(writer, "\tinner-iteration= %6d\tpseudo-cfl=%10.3e\tpseudo-dt=%10.3e\tglobal-relative-residual=%.12e\tglobal-absolute-residual=%.12e",
                   step, cfl, dt, globalResidual.re/referenceGlobalResidual.re, globalResidual.re);
    writeln(writer.data);
}
