/** steadystate_core.d
 * Core set of functions used in the Newton-Krylov updates for steady-state convergence.
 *
 * Author: Rowan G.
 * Date: 2016-10-09
 */

module steadystate_core;

import core.stdc.stdlib : exit;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.array;
import std.math;
import std.datetime;

import nm.smla;
import nm.bbla;
import nm.complex;
import nm.number;

import gas;
import fvcell;
import fvinterface;
import geom;
import special_block_init;
import fluidblock;
import fluidblockio;
import fluidblockio_old;
import fluidblockio_new;
import sfluidblock;
import globaldata;
import globalconfig;
import simcore;
import simcore_gasdynamic_step;
import simcore_solid_step;
import simcore_exchange;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess : readTimesFile;
import loads;
//import shape_sensitivity_core : sss_preconditioner_initialisation, sss_preconditioner;
import jacobian;
import solid_loose_coupling_update;

version(mpi_parallel) {
    import mpi;
}


static int fnCount = 0;

// Module-local, global memory arrays and matrices
number[] g0;
number[] g1;
number[] h;
number[] hR;
Matrix!number H0;
Matrix!number H1;
Matrix!number Gamma;
Matrix!number Q0;
Matrix!number Q1;

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

void extractRestartInfoFromTimesFile(string jobName, ref RestartInfo[] times)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;
    size_t SPECIES = GlobalConfig.cqi.species;
    size_t MODES = GlobalConfig.cqi.modes;
    auto cqi = GlobalConfig.cqi;

    auto gmodel = GlobalConfig.gmodel_master;
    RestartInfo restartInfo = RestartInfo(GlobalConfig.cqi.n);
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
            restartInfo.residuals.vec[cqi.mass] = to!double(tokens[6+MASS]);
            restartInfo.residuals.vec[cqi.xMom] = to!double(tokens[6+X_MOM]);
            restartInfo.residuals.vec[cqi.yMom] = to!double(tokens[6+Y_MOM]);
            if ( GlobalConfig.dimensions == 3 )
                restartInfo.residuals.vec[cqi.zMom] = to!double(tokens[6+Z_MOM]);
            restartInfo.residuals.vec[cqi.totEnergy] = to!double(tokens[6+TOT_ENERGY]);
            foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
                restartInfo.residuals.vec[cqi.rhoturb+it] = to!double(tokens[6+TKE+it]);
            }
            version(multi_species_gas){
            foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                restartInfo.residuals.vec[cqi.species+sp] = to!double(tokens[6+SPECIES+sp]);
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                restartInfo.residuals.vec[cqi.modes+imode] = to!double(tokens[6+MODES+imode]);
            }
            }
            times ~= restartInfo;
        }
        line = timesFile.readln().strip();
    }
    timesFile.close();
    return;
}


void iterate_to_steady_state(int snapshotStart, int maxCPUs, int threadsPerMPITask)
{
    auto wallClockStart = Clock.currTime();
    string jobName = GlobalConfig.base_file_name;
    int nsteps = GlobalConfig.sssOptions.nTotalSteps;
    int nIters = 0;
    int nRestarts;
    double linSolResid = 0;
    double relGlobalResidReduction = GlobalConfig.sssOptions.stopOnRelGlobalResid;
    double absGlobalResidReduction = GlobalConfig.sssOptions.stopOnAbsGlobalResid;
    double cfl_max = GlobalConfig.sssOptions.cfl_max;
    int cfl_schedule_current_index = 0;
    auto cfl_schedule_value_list = GlobalConfig.sssOptions.cfl_schedule_value_list;
    auto cfl_schedule_iter_list = GlobalConfig.sssOptions.cfl_schedule_iter_list;
    bool residual_based_cfl_scheduling = GlobalConfig.sssOptions.residual_based_cfl_scheduling;
    int kmax = GlobalConfig.sssOptions.maxSubIterations;
    // Settings for start-up phase
    int LHSeval0 = GlobalConfig.sssOptions.LHSeval0;
    int RHSeval0 = GlobalConfig.sssOptions.RHSeval0;
    double cfl0 = GlobalConfig.sssOptions.cfl0;
    double tau0 = GlobalConfig.sssOptions.tau0;
    double eta0 = GlobalConfig.sssOptions.eta0;
    double sigma0 = GlobalConfig.sssOptions.sigma0;
    // Settings for inexact Newton phase
    int LHSeval1;
    if (GlobalConfig.interpolation_order < 2) { LHSeval1 = GlobalConfig.interpolation_order; }
    else { LHSeval1 = GlobalConfig.sssOptions.LHSeval1; }
    int RHSeval1 = GlobalConfig.interpolation_order; // GlobalConfig.sssOptions.RHSeval1
    double cfl1 = GlobalConfig.sssOptions.cfl1;
    double tau1 = GlobalConfig.sssOptions.tau1;
    double sigma1 = GlobalConfig.sssOptions.sigma1;
    EtaStrategy etaStrategy = GlobalConfig.sssOptions.etaStrategy;
    double eta1 = GlobalConfig.sssOptions.eta1;
    double eta1_max = GlobalConfig.sssOptions.eta1_max;
    double eta1_min = GlobalConfig.sssOptions.eta1_min;
    double etaRatioPerStep = GlobalConfig.sssOptions.etaRatioPerStep;
    double gamma = GlobalConfig.sssOptions.gamma;
    double alpha = GlobalConfig.sssOptions.alpha;
    double limiterFreezingResidReduction = GlobalConfig.sssOptions.limiterFreezingResidReduction;
    int limiterFreezingCount = GlobalConfig.sssOptions.limiterFreezingCount;
    int countsBeforeFreezing = 0;
    bool limiterFreezingCondition = false;
    int interpOrderSave = GlobalConfig.interpolation_order;

    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;
    size_t SPECIES = GlobalConfig.cqi.species;
    size_t MODES = GlobalConfig.cqi.modes;

    ConservedQuantities maxResiduals = new ConservedQuantities(nConserved);
    ConservedQuantities currResiduals = new ConservedQuantities(nConserved);
    number mass_balance = 0.0;
    number omega = 1.0;
    number theta = GlobalConfig.sssOptions.physicalityCheckTheta;
    bool pc_matrix_evaluated = false;

    double cfl, cflTrial;
    double dt;
    double etaTrial;
    double pseudoSimTime = 0.0;
    double normOld, normNew;
    int snapshotsCount = GlobalConfig.sssOptions.snapshotsCount;
    int nTotalSnapshots = GlobalConfig.sssOptions.nTotalSnapshots;
    int nWrittenSnapshots = 0;
    int writeDiagnosticsCount = GlobalConfig.sssOptions.writeDiagnosticsCount;
    int writeLoadsCount = GlobalConfig.sssOptions.writeLoadsCount;

    int startStep;
    int nStartUpSteps = GlobalConfig.sssOptions.nStartUpSteps;
    bool inexactNewtonPhase = false;
    bool dangerousExceptionsAreFatal = GlobalConfig.dangerous_exceptions_are_fatal;

    // No need to have more task threads than blocks
    int extraThreadsInPool;
    auto nBlocksInThreadParallel = localFluidBlocks.length;
    version(mpi_parallel) {
        extraThreadsInPool = min(threadsPerMPITask-1, nBlocksInThreadParallel-1);
    } else {
        extraThreadsInPool = min(maxCPUs-1, nBlocksInThreadParallel-1);
    }
    defaultPoolThreads(extraThreadsInPool);
    version(mpi_parallel) {
        writefln("MPI-task %d : running with %d threads.", GlobalConfig.mpi_rank_for_local_task, extraThreadsInPool+1);
    }
    else {
        writefln("Single process running with %d threads.", extraThreadsInPool+1); // +1 for main thread.
    }
    double normRef = 0.0;
    bool residualsUpToDate = false;
    bool finalStep = false;
    bool usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
    if (usePreconditioner) {
        evalRHS(0.0, 0);
        // initialize the flow Jacobians used as local precondition matrices for GMRES
        final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
            case PreconditionMatrixType.jacobi:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(-1, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.ilu:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.sgs:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.sgs_relax:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.lu_sgs:
                // do nothing
                break;
        } // end switch
        //foreach (blk; localFluidBlocks) { blk.verify_jacobian(GlobalConfig.sssOptions.preconditionerSigma); }
    }

    // We need to calculate the initial residual if we are starting from scratch.
    if ( snapshotStart == 0 ) {
        evalRHS(0.0, 0);
        max_residuals(maxResiduals);
        foreach (blk; parallel(localFluidBlocks, 1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            size_t nsp = blk.myConfig.gmodel.n_species;
            size_t nmodes = blk.myConfig.gmodel.n_modes;
                auto cqi = blk.myConfig.cqi;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.FU[cellCount+MASS] = -cell.dUdt[0].vec[cqi.mass];
                    blk.FU[cellCount+X_MOM] = -cell.dUdt[0].vec[cqi.xMom];
                    blk.FU[cellCount+Y_MOM] = -cell.dUdt[0].vec[cqi.yMom];
                    if ( GlobalConfig.dimensions == 3 )
                        blk.FU[cellCount+Z_MOM] = -cell.dUdt[0].vec[cqi.zMom];
                    blk.FU[cellCount+TOT_ENERGY] = -cell.dUdt[0].vec[cqi.totEnergy];
                    foreach(it; 0 .. nturb) {
                        blk.FU[cellCount+TKE+it] = -cell.dUdt[0].vec[cqi.rhoturb+it];
                    }
                    version(multi_species_gas){
                        if ( nsp > 1 ) {
                            foreach(sp; 0 .. nsp) { blk.FU[cellCount+SPECIES+sp] = -cell.dUdt[0].vec[cqi.species+sp]; }
                        }
                    }
                    version(multi_T_gas){
                        foreach(imode; 0 .. nmodes) { blk.FU[cellCount+MODES+imode] = -cell.dUdt[0].vec[cqi.modes+imode]; }
                    }
                    cellCount += nConserved;
                }
        }
        mixin(dot_over_blocks("normRef", "FU", "FU"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(normRef), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        normRef = sqrt(normRef);
    }

    // if we are restarting a simulation we need read the initial residuals from a file
    RestartInfo[] times;
    if (snapshotStart > 0) { // && GlobalConfig.is_master_task) {
        extractRestartInfoFromTimesFile(jobName, times);
        normOld = times[snapshotStart].globalResidual;
        // We need to read in the reference residual values from a file.
        string refResidFname = jobName ~ "-ref-residuals.saved";
        auto refResid = File(refResidFname, "r");
        auto line = refResid.readln().strip();
        auto tokens = line.split();
        normRef = to!double(tokens[0]);
        auto cqi = GlobalConfig.cqi;
        maxResiduals.vec[cqi.mass] = to!double(tokens[1+MASS]);
        maxResiduals.vec[cqi.xMom] = to!double(tokens[1+X_MOM]);
        maxResiduals.vec[cqi.yMom] = to!double(tokens[1+Y_MOM]);
        if ( GlobalConfig.dimensions == 3 )
            maxResiduals.vec[cqi.zMom] = to!double(tokens[1+Z_MOM]);
        maxResiduals.vec[cqi.totEnergy] = to!double(tokens[1+TOT_ENERGY]);
        foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
            maxResiduals.vec[cqi.rhoturb+it] = to!double(tokens[1+TKE+it]);
        }
        version(multi_species_gas){
            if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                    maxResiduals.vec[cqi.species+sp] = to!double(tokens[1+SPECIES+sp]);
                }
            }
        }
        version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                maxResiduals.vec[cqi.modes+imode] = to!double(tokens[1+MODES+imode]);
            }
        }

        // We also need to determine how many snapshots have already been written
        auto timesFile = File("./config/" ~ jobName ~ ".times");
        line = timesFile.readln().strip();
        while (line.length > 0) {
            if (line[0] != '#') {
                nWrittenSnapshots++;
            }
            line = timesFile.readln().strip();
        }
        timesFile.close();
        nWrittenSnapshots--; // We don't count the initial solution as a written snapshot

        // Where should we be in the CFL schedule list?
        if (!residual_based_cfl_scheduling) {
            while ((cfl_schedule_current_index < cfl_schedule_value_list.length) &&
                   (times[snapshotStart].step > cfl_schedule_iter_list[cfl_schedule_current_index])) {
                cfl_schedule_current_index++;
            }
        }
    }

    if (GlobalConfig.is_master_task) {
        // report the initial residuals
        auto cqi = GlobalConfig.cqi;
        writeln("Reference residuals are established as:");
        writefln("GLOBAL:         %.12e", normRef);
        writefln("MASS:           %.12e", maxResiduals.vec[cqi.mass].re);
        writefln("X-MOMENTUM:     %.12e", maxResiduals.vec[cqi.xMom].re);
        writefln("Y-MOMENTUM:     %.12e", maxResiduals.vec[cqi.yMom].re);
        if ( GlobalConfig.dimensions == 3 )
            writefln("Z-MOMENTUM:     %.12e", maxResiduals.vec[cqi.zMom].re);
        writefln("ENERGY:         %.12e", maxResiduals.vec[cqi.totEnergy].re);
        foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
            string tvname = capitalize(GlobalConfig.turb_model.primitive_variable_name(it));
            writefln("%s:            %.12e",tvname, maxResiduals.vec[cqi.rhoturb+it].re);
        }
        version(multi_species_gas){
            if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                    string spname = capitalize(GlobalConfig.gmodel_master.species_name(sp));
                    writefln("%s:            %.12e",spname, maxResiduals.vec[cqi.species+sp].re);
                }
            }
        }
        version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                string modename = "T_MODES["~to!string(imode)~"]"; //capitalize(GlobalConfig.gmodel_master.energy_mode_name(imode));
                writefln("%s:            %.12e",modename, maxResiduals.vec[cqi.modes+imode].re);
            }
        }
        // store the initial residuals
        string refResidFname = jobName ~ "-ref-residuals.saved";
        auto refResid = File(refResidFname, "w");
        if ( GlobalConfig.dimensions == 2 ) {
            refResid.writef("%.18e %.18e %.18e %.18e %.18e",
                            normRef, maxResiduals.vec[cqi.mass].re, maxResiduals.vec[cqi.xMom].re,
                            maxResiduals.vec[cqi.yMom].re, maxResiduals.vec[cqi.totEnergy].re);
        }
        else {
            refResid.writef("%.18e %.18e %.18e %.18e %.18e %.18e",
                            normRef, maxResiduals.vec[cqi.mass].re, maxResiduals.vec[cqi.xMom].re,
                            maxResiduals.vec[cqi.yMom].re, maxResiduals.vec[cqi.zMom].re,
                            maxResiduals.vec[cqi.totEnergy].re);
        }
        foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
            refResid.writef(" %.18e", maxResiduals.vec[cqi.rhoturb+it].re);
        }
        version(multi_species_gas){
            if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                    refResid.writef(" %.18e", maxResiduals.vec[cqi.species+sp].re);
                }
            }
        }
        version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                refResid.writef(" %.18e", maxResiduals.vec[cqi.modes+imode].re);
            }
        }
        refResid.write("\n");
        refResid.close();
    }

    double wallClockElapsed;
    RestartInfo restartInfo;
    // We need to do some configuration based on whether we are starting from scratch,
    // or attempting to restart from an earlier snapshot.
    if ( snapshotStart == 0 ) {
        startStep = 1;
        restartInfo.pseudoSimTime = 0.0;
        restartInfo.dt = dt;
        restartInfo.cfl = cfl;
        restartInfo.step = 0;
        restartInfo.globalResidual = normRef;
        restartInfo.residuals = maxResiduals;
        times ~= restartInfo;
    } else {
        restartInfo = times[snapshotStart];
        dt = restartInfo.dt;
        cfl = restartInfo.cfl;
        startStep = restartInfo.step + 1;
        pseudoSimTime = restartInfo.pseudoSimTime;
        if ( GlobalConfig.is_master_task ) {
            writefln("Restarting steps from step= %d", startStep);
            writefln("   pseudo-sim-time= %.6e dt= %.6e", pseudoSimTime, dt);
        }
    }

    // setup the diagnostics file
    auto residFname = "e4-nk.diagnostics.dat";
    File fResid;
    if (GlobalConfig.is_master_task) {
        if ( snapshotStart == 0 ) {
            // Open file for writing diagnostics
            fResid = File(residFname, "w");
            fResid.writeln("#  1: step");
            fResid.writeln("#  2: pseudo-time");
            fResid.writeln("#  3: dt");
            fResid.writeln("#  4: CFL");
            fResid.writeln("#  5: eta");
            fResid.writeln("#  6: nRestarts");
            fResid.writeln("#  7: nIters");
            fResid.writeln("#  8: nFnCalls");
            fResid.writeln("#  9: wall-clock, s");
            fResid.writeln("#  10: global-residual-abs");
            fResid.writeln("# 11: global-residual-rel");
            fResid.writefln("#  %02d: mass-abs", 12+2*MASS);
            fResid.writefln("# %02d: mass-rel", 12+2*MASS+1);
            fResid.writefln("# %02d: x-mom-abs", 12+2*X_MOM);
            fResid.writefln("# %02d: x-mom-rel", 12+2*X_MOM+1);
            fResid.writefln("# %02d: y-mom-abs", 12+2*Y_MOM);
            fResid.writefln("# %02d: y-mom-rel", 12+2*Y_MOM+1);
            if ( GlobalConfig.dimensions == 3 ) {
                fResid.writefln("# %02d: z-mom-abs", 12+2*Z_MOM);
                fResid.writefln("# %02d: z-mom-rel", 12+2*Z_MOM+1);
            }
            fResid.writefln("# %02d: energy-abs", 12+2*TOT_ENERGY);
            fResid.writefln("# %02d: energy-rel", 12+2*TOT_ENERGY+1);
            auto nt = GlobalConfig.turb_model.nturb;
            foreach(it; 0 .. nt) {
                string tvname = GlobalConfig.turb_model.primitive_variable_name(it);
                fResid.writefln("# %02d: %s-abs", 12+2*(TKE+it), tvname);
                fResid.writefln("# %02d: %s-rel", 12+2*(TKE+it)+1, tvname);
            }
            auto nsp = GlobalConfig.gmodel_master.n_species;
            if ( nsp > 1) {
                foreach(sp; 0 .. nsp) {
                    string spname = GlobalConfig.gmodel_master.species_name(sp);
                    fResid.writefln("# %02d: %s-abs", 12+2*(SPECIES+sp), spname);
                    fResid.writefln("# %02d: %s-rel", 12+2*(SPECIES+sp)+1, spname);
                }
            }
            auto nmodes = GlobalConfig.gmodel_master.n_modes;
            foreach(imode; 0 .. nmodes) {
                string modename = "T_MODES["~to!string(imode)~"]"; //GlobalConfig.gmodel_master.energy_mode_name(imode);
                fResid.writefln("# %02d: %s-abs", 12+2*(MODES+imode), modename);
                fResid.writefln("# %02d: %s-rel", 12+2*(MODES+imode)+1, modename);
            }
            fResid.writefln("# %02d: mass-balance", 12+2*(max(TOT_ENERGY+1, TKE+nt,SPECIES+nsp,MODES+nmodes)));
            fResid.writefln("# %02d: linear-solve-residual", 1+12+2*(max(TOT_ENERGY+1, TKE+nt,SPECIES+nsp,MODES+nmodes)));
            if (GlobalConfig.sssOptions.useLineSearch || GlobalConfig.sssOptions.usePhysicalityCheck) {
                fResid.writefln("# %02d: omega", 2+12+2*(max(TOT_ENERGY+1, TKE+nt,SPECIES+nsp,MODES+nmodes)));
            }
            fResid.writefln("# %02d: PC", 3+12+2*(max(TOT_ENERGY+1, TKE+nt,SPECIES+nsp,MODES+nmodes)));
            fResid.close();
        }
    }

    // calculate an initial timestep
    cfl = cfl0;
    dt = determine_dt(cfl);
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }

    // Begin Newton steps
    int LHSeval;
    int RHSeval;
    double eta;
    double tau;
    double sigma;
    foreach (step; startStep .. nsteps+1) {
        SimState.step = step;
        if ( (step/GlobalConfig.control_count)*GlobalConfig.control_count == step ) {
            read_control_file(); // Reparse the time-step control parameters occasionally.
        }
        residualsUpToDate = false;
        if ( step <= nStartUpSteps ) {
            LHSeval = LHSeval0;
            RHSeval = RHSeval0;
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(RHSeval);
            eta = eta0;
            tau = tau0;
            sigma = sigma0;
            usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
        }
        else {
            LHSeval = LHSeval1;
            RHSeval = RHSeval1;
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(RHSeval);
            eta = eta1;
            tau = tau1;
            sigma = sigma1;
            usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
            if (step > GlobalConfig.freeze_limiter_on_step) {
                // compute the limiter value once more before freezing it
                if (GlobalConfig.frozen_limiter == false) {
                    evalRHS(pseudoSimTime, 0);
                    GlobalConfig.frozen_limiter = true;
                }
            }
        }

        // solve linear system
        version(lu_sgs) { lusgs_solve(step, pseudoSimTime, dt, normNew, startStep); }
        else { rpcGMRES_solve(step, pseudoSimTime, dt, eta, sigma, usePreconditioner, normNew, nRestarts, nIters, linSolResid, startStep, LHSeval, RHSeval, pc_matrix_evaluated); }

        // calculate a relaxation factor for the nonlinear update via a physicality check
        omega = 1.0;
        if (GlobalConfig.sssOptions.usePhysicalityCheck) {
            foreach (blk; parallel(localFluidBlocks,1)) {
                int cellCount = 0;
                auto cqi = blk.myConfig.cqi;
                foreach (cell; blk.cells) {
                    // limit the change in the conserved mass by some user defined amount (theta)
                    number rel_diff_limit = fabs(blk.dU[cellCount+MASS]/(theta*cell.U[0].vec[cqi.mass]));
                    omega = 1.0/(max(rel_diff_limit,1.0/omega));
                    // ensure thermodynamic state variables are positive
                    bool failed_decode = true;
                    while (failed_decode) {
                        failed_decode = false;
                        cell.U[1].copy_values_from(cell.U[0]);
                        foreach (j; 0 .. cqi.n) {
                            cell.U[1].vec[j] = cell.U[0].vec[j] + omega*blk.dU[cellCount+j];
                        }
                        try {
                            cell.decode_conserved(0, 1, 0.0);
                        }
                        catch (FlowSolverException e) {
                            failed_decode = true;
                        }
                        // reduce relaxation factor if appropriate (0.7 seems like a resonable reduction factor)
                        if (failed_decode) omega *= 0.7;
                        // return cell to original state
                        cell.decode_conserved(0, 0, 0.0);
                    }
                    cellCount += nConserved;
                }
            }
            // communicate minimum omega to all threads
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &(omega.re), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
        } // end physicality check

        if (omega < 0.01)  {
            if ( GlobalConfig.is_master_task ) {
                writefln("WARNING: nonlinear update relaxation factor too small for step= %d", step);
            }

            // update isn't good, reduce CFL and try again (a reduction factor of 0.5 is aggressie, could lower it)
            cfl = 0.5*cfl;
            dt = determine_dt(cfl);
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }

            // return cell flow-states to their original state
            foreach (blk; parallel(localFluidBlocks,1)) {
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    cell.decode_conserved(0, 0, 0.0);
                }
            }

            if (dangerousExceptionsAreFatal) {
                writefln("         dangerousExceptionsAreFatal set to true, bailing out!");
                exit(1);
            } else {
                if ( GlobalConfig.is_master_task ) {
                    writefln("         Taking nonlinear step again with the CFL reduced by a factor of 0.5");
                }
            }

            // we don't procede with the nonlinear update for this step
            continue;
        }

        // If we get here, things are good
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            size_t nsp = blk.myConfig.gmodel.n_species;
            size_t nmodes = blk.myConfig.gmodel.n_modes;
            int cellCount = 0;
            auto cqi = blk.myConfig.cqi;
            foreach (cell; blk.cells) {
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].vec[cqi.mass] = cell.U[0].vec[cqi.mass] + omega*blk.dU[cellCount+MASS];
                cell.U[1].vec[cqi.xMom] = cell.U[0].vec[cqi.xMom] + omega*blk.dU[cellCount+X_MOM];
                cell.U[1].vec[cqi.yMom] = cell.U[0].vec[cqi.yMom] + omega*blk.dU[cellCount+Y_MOM];
                if ( blk.myConfig.dimensions == 3 )
                    cell.U[1].vec[cqi.zMom] = cell.U[0].vec[cqi.zMom] + omega*blk.dU[cellCount+Z_MOM];
                cell.U[1].vec[cqi.totEnergy] = cell.U[0].vec[cqi.totEnergy] + omega*blk.dU[cellCount+TOT_ENERGY];
                foreach(it; 0 .. nturb){
                    cell.U[1].vec[cqi.rhoturb+it] = cell.U[0].vec[cqi.rhoturb+it] + omega*blk.dU[cellCount+TKE+it];
                }
                version(multi_species_gas){
                    if (blk.myConfig.n_species > 1) {
                        foreach(sp; 0 .. nsp) { cell.U[1].vec[cqi.species+sp] = cell.U[0].vec[cqi.species+sp] + omega*blk.dU[cellCount+SPECIES+sp]; }
                    } else {
                        // enforce mass fraction of 1 for single species gas
                        cell.U[1].vec[cqi.species+0] = cell.U[1].vec[cqi.mass];
                    }
                }
                version(multi_T_gas){
                    foreach(imode; 0 .. nmodes) { cell.U[1].vec[cqi.modes+imode] = cell.U[0].vec[cqi.modes+imode] + omega*blk.dU[cellCount+MODES+imode]; }
                }
                cell.decode_conserved(0, 1, 0.0);
                cellCount += nConserved;
            }
        }

        // Put flow state into U[0] ready for next iteration.
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (cell; blk.cells) {
                swap(cell.U[0], cell.U[1]);
            }
        }

        // after a successful fluid domain update, proceed to perform a solid domain update
        if (localSolidBlocks.length > 0) { solid_update(step, pseudoSimTime, cfl, eta, sigma); }

        pseudoSimTime += dt;
        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();

	if (!limiterFreezingCondition && (normNew/normRef <= limiterFreezingResidReduction)) {
	    countsBeforeFreezing++;
	    if (countsBeforeFreezing >= limiterFreezingCount) {
                if (GlobalConfig.frozen_limiter == false) {
                    evalRHS(pseudoSimTime, 0);
                    GlobalConfig.frozen_limiter = true;
                }
                limiterFreezingCondition = true;
                writefln("=== limiter freezing condition met at step: %d ===", step);
            }
	}
        // Check on some stopping criteria
        if ( step == nsteps ) {
            if (GlobalConfig.is_master_task) {
                writeln("STOPPING: Reached maximum number of steps.");
            }
            finalStep = true;
        }
        if ( normNew <= absGlobalResidReduction && step > nStartUpSteps ) {
            if (GlobalConfig.is_master_task) {
                writeln("STOPPING: The absolute global residual is below target value.");
                writefln("          current value= %.12e   target value= %.12e", normNew, absGlobalResidReduction);
            }
            finalStep = true;
        }
        if ( normNew/normRef <= relGlobalResidReduction && step > nStartUpSteps ) {
            if (GlobalConfig.is_master_task) {
                writeln("STOPPING: The relative global residual is below target value.");
                writefln("          current value= %.12e   target value= %.12e", normNew/normRef, relGlobalResidReduction);
            }
            finalStep = true;
        }
	if (GlobalConfig.halt_now == 1) {
            if (GlobalConfig.is_master_task) {
                writeln("STOPPING: Halt set in control file.");
            }
            finalStep = true;
        }

        // Now do some output and diagnostics work
        if ( (step % writeDiagnosticsCount) == 0 || finalStep ) {
            mass_balance = 0.0;
            compute_mass_balance(mass_balance);
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &(mass_balance.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            }
            // Write out residuals
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            if (GlobalConfig.is_master_task) {
                auto cqi = GlobalConfig.cqi;
                fResid = File(residFname, "a");
                fResid.writef("%8d  %20.16e  %20.16e %20.16e %20.16e %3d %3d %5d %.8f %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e ",
                              step, pseudoSimTime, dt, cfl, eta, nRestarts, nIters, fnCount, wallClockElapsed,
                              normNew, normNew/normRef,
                              currResiduals.vec[cqi.mass].re, currResiduals.vec[cqi.mass].re/maxResiduals.vec[cqi.mass].re,
                              currResiduals.vec[cqi.xMom].re, currResiduals.vec[cqi.xMom].re/maxResiduals.vec[cqi.xMom].re,
                              currResiduals.vec[cqi.yMom].re, currResiduals.vec[cqi.yMom].re/maxResiduals.vec[cqi.yMom].re);
                if ( GlobalConfig.dimensions == 3 )
                    fResid.writef("%20.16e  %20.16e  ", currResiduals.vec[cqi.zMom].re, currResiduals.vec[cqi.zMom].re/maxResiduals.vec[cqi.zMom].re);
                fResid.writef("%20.16e  %20.16e  ",
                              currResiduals.vec[cqi.totEnergy].re, currResiduals.vec[cqi.totEnergy].re/maxResiduals.vec[cqi.totEnergy].re);
                foreach(it; 0 .. GlobalConfig.turb_model.nturb){
                    fResid.writef("%20.16e  %20.16e  ",
                                  currResiduals.vec[cqi.rhoturb+it].re, currResiduals.vec[cqi.rhoturb+it].re/maxResiduals.vec[cqi.rhoturb+it].re);
                }
                version(multi_species_gas){
                if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                    foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                        fResid.writef("%20.16e  %20.16e  ",
                                      currResiduals.vec[cqi.species+sp].re, currResiduals.vec[cqi.species+sp].re/maxResiduals.vec[cqi.species+sp].re);
                    }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
                    fResid.writef("%20.16e  %20.16e  ",
                                  currResiduals.vec[cqi.modes+imode].re, currResiduals.vec[cqi.modes+imode].re/maxResiduals.vec[cqi.modes+imode].re);
                }
                }
                fResid.writef("%20.16e ", fabs(mass_balance.re));
                fResid.writef("%20.16e ", linSolResid);
                if (GlobalConfig.sssOptions.useLineSearch || GlobalConfig.sssOptions.usePhysicalityCheck) {
                    fResid.writef("%20.16e ", omega.re);
                }
                fResid.writef("%d ", pc_matrix_evaluated);
                fResid.write("\n");
                fResid.close();
            }
        }

        // write out the loads
        if ( (step % writeLoadsCount) == 0 || finalStep || step == GlobalConfig.write_loads_at_step) {
            if (GlobalConfig.is_master_task) {
                init_current_loads_tindx_dir(step);
            }
            version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
            wait_for_current_tindx_dir(step);
            write_boundary_loads_to_file(pseudoSimTime, step);
            if (GlobalConfig.is_master_task) {
                update_loads_times_file(pseudoSimTime, step);
            }
        }

        if ( (step % GlobalConfig.print_count) == 0 || finalStep ) {
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            if (GlobalConfig.is_master_task) {
                auto cqi = GlobalConfig.cqi;
                auto writer = appender!string();
                formattedWrite(writer, "STEP= %7d  pseudo-time=%10.3e dt=%10.3e cfl=%10.3e  WC=%.1f \n", step, pseudoSimTime, dt, cfl, wallClockElapsed);
                formattedWrite(writer, "RESIDUALS        absolute        relative\n");
                formattedWrite(writer, "  global         %10.6e    %10.6e\n", normNew, normNew/normRef);
                formattedWrite(writer, "  mass           %10.6e    %10.6e\n", currResiduals.vec[cqi.mass].re, currResiduals.vec[cqi.mass].re/maxResiduals.vec[cqi.mass].re);
                formattedWrite(writer, "  x-mom          %10.6e    %10.6e\n", currResiduals.vec[cqi.xMom].re, currResiduals.vec[cqi.xMom].re/maxResiduals.vec[cqi.xMom].re);
                formattedWrite(writer, "  y-mom          %10.6e    %10.6e\n", currResiduals.vec[cqi.yMom].re, currResiduals.vec[cqi.yMom].re/maxResiduals.vec[cqi.yMom].re);
                if ( GlobalConfig.dimensions == 3 )
                    formattedWrite(writer, "  z-mom          %10.6e    %10.6e\n", currResiduals.vec[cqi.zMom].re, currResiduals.vec[cqi.zMom].re/maxResiduals.vec[cqi.zMom].re);
                formattedWrite(writer, "  total-energy   %10.6e    %10.6e\n", currResiduals.vec[cqi.totEnergy].re, currResiduals.vec[cqi.totEnergy].re/maxResiduals.vec[cqi.totEnergy].re);
                foreach(it; 0 .. GlobalConfig.turb_model.nturb){
                    auto tvname = GlobalConfig.turb_model.primitive_variable_name(it);
                    formattedWrite(writer, "  %s            %10.6e    %10.6e\n", tvname, currResiduals.vec[cqi.rhoturb+it].re, currResiduals.vec[cqi.rhoturb+it].re/maxResiduals.vec[cqi.rhoturb+it].re);
                }
                version(multi_species_gas){
                if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                    foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                        auto spname = GlobalConfig.gmodel_master.species_name(sp);
                        formattedWrite(writer, "  %s            %10.6e    %10.6e\n", spname, currResiduals.vec[cqi.species+sp].re, currResiduals.vec[cqi.species+sp].re/maxResiduals.vec[cqi.species+sp].re);
                    }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
                    auto modename = "T_MODES["~to!string(imode)~"]"; //GlobalConfig.gmodel_master.energy_mode_name(imode);
                    formattedWrite(writer, "  %s            %10.6e    %10.6e\n", modename, currResiduals.vec[cqi.modes+imode].re, currResiduals.vec[cqi.modes+imode].re/maxResiduals.vec[cqi.modes+imode].re);
                }
                }
                writeln(writer.data);
            }
        }

        // Write out the flow field, if required
        if ( (step % snapshotsCount) == 0 || finalStep || step == GlobalConfig.write_flow_solution_at_step) {
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            if (GlobalConfig.is_master_task) {
                writefln("-----------------------------------------------------------------------");
                writefln("Writing flow solution at step= %4d; pseudo-time= %6.3e", step, pseudoSimTime);
                writefln("-----------------------------------------------------------------------\n");
            }
            FluidBlockIO[] io_list = localFluidBlocks[0].block_io;
            bool legacy = is_legacy_format(GlobalConfig.flow_format);
            nWrittenSnapshots++;
            if ( nWrittenSnapshots <= nTotalSnapshots ) {
                if (GlobalConfig.is_master_task){
                    if (legacy) {
                        ensure_directory_is_present(make_path_name!"flow"(nWrittenSnapshots));
                        ensure_directory_is_present(make_path_name!"solid"(nWrittenSnapshots));
                    } else {
                        foreach(io; io_list) {
                            string path = "CellData/"~io.tag;
                            if (io.do_save()) ensure_directory_is_present(make_path_name(path, nWrittenSnapshots));
                        }
                    }
                }
                version(mpi_parallel) {
                    MPI_Barrier(MPI_COMM_WORLD);
                }
                foreach (blk; localFluidBlocks) {
                    if (legacy) {
                        auto fileName = make_file_name!"flow"(jobName, blk.id, nWrittenSnapshots, "gz");
                        blk.write_solution(fileName, pseudoSimTime);
                    } else {
                        foreach(io; blk.block_io) {
                            auto fileName = make_file_name("CellData", io.tag, jobName, blk.id, nWrittenSnapshots, GlobalConfig.flowFileExt);
                            if (io.do_save()) io.save_to_file(fileName, pseudoSimTime);
                        }
                    }
                }
                foreach (sblk; localSolidBlocks) {
                    auto fileName = make_file_name!"solid"(jobName, sblk.id, nWrittenSnapshots, "gz");
                    sblk.writeSolution(fileName, pseudoSimTime);
                }
                restartInfo.pseudoSimTime = pseudoSimTime;
                restartInfo.dt = dt;
                restartInfo.cfl = cfl;
                restartInfo.step = step;
                restartInfo.globalResidual = normNew;
                restartInfo.residuals = currResiduals;
                times ~= restartInfo;
                if (GlobalConfig.is_master_task) { rewrite_times_file(times); }
            }
            else {
                // We need to shuffle all of the fluid snapshots...
                foreach ( iSnap; 2 .. nTotalSnapshots+1) {
                    foreach (blk; localFluidBlocks) {
                        if (legacy) {
                            auto fromName = make_file_name!"flow"(jobName, blk.id, iSnap, "gz");
                            auto toName = make_file_name!"flow"(jobName, blk.id, iSnap-1, "gz");
                            rename(fromName, toName);
                        } else {
                            foreach (io; io_list) {
                                auto fromName = make_file_name("CellData", io.tag, jobName, blk.id, iSnap, GlobalConfig.flowFileExt);
                                auto toName = make_file_name("CellData", io.tag, jobName, blk.id, iSnap-1, GlobalConfig.flowFileExt);
                                rename(fromName, toName);
                            }
                        }
                    }
                }
                // ... and add the new fluid snapshot.
                foreach (blk; localFluidBlocks) {
                    if (legacy) {
                        auto fileName = make_file_name!"flow"(jobName, blk.id, nTotalSnapshots, "gz");
                        blk.write_solution(fileName, pseudoSimTime);
                    } else {
                        foreach(io; blk.block_io) {
                            auto fileName = make_file_name("CellData", io.tag, jobName, blk.id, nTotalSnapshots, GlobalConfig.flowFileExt);
                            if (io.do_save()) io.save_to_file(fileName, pseudoSimTime);

                        }
                    }
                }
                // We need to shuffle all of the solid snapshots...
                foreach ( iSnap; 2 .. nTotalSnapshots+1) {
                    foreach (blk; localSolidBlocks) {
                        auto fromName = make_file_name!"solid"(jobName, blk.id, iSnap, "gz");
                        auto toName = make_file_name!"solid"(jobName, blk.id, iSnap-1, "gz");
                        rename(fromName, toName);
                    }
                }
                // ... and add the new solid snapshot.
                foreach (sblk; localSolidBlocks) {
                    auto fileName = make_file_name!"solid"(jobName, sblk.id, nTotalSnapshots, "gz");
                    sblk.writeSolution(fileName, pseudoSimTime);
                }
                remove(times, 1);
                restartInfo.pseudoSimTime = pseudoSimTime;
                restartInfo.dt = dt;
                restartInfo.cfl = cfl;                
                restartInfo.step = step;
                restartInfo.globalResidual = normNew;
                restartInfo.residuals = currResiduals;
                times[$-1] = restartInfo;
                if (GlobalConfig.is_master_task) { rewrite_times_file(times); }
            }
        }

        if (finalStep) break;
        
        if (!inexactNewtonPhase && normNew/normRef < tau ) {
            // Switch to inexactNewtonPhase
            inexactNewtonPhase = true;
        }

        // Choose a new timestep and eta value.
        auto normRatio = normOld/normNew;
        if (residual_based_cfl_scheduling) { 
            if (inexactNewtonPhase) {
                if (omega >= 0.1) {
                    if (step < nStartUpSteps) {
                        // Let's assume we're still letting the shock settle
                        // when doing low order steps, so we use a power of 0.75 as a default
                        double p0 =  GlobalConfig.sssOptions.p0;
                        cflTrial = cfl*pow(normOld/normNew, p0);
                    }
                    else {
                        // We use a power of 1.0 as a default
                        double p1 =  GlobalConfig.sssOptions.p1;
                        cflTrial = cfl*pow(normOld/normNew, p1);
                    }
                    // Apply safeguards to dt
                    cflTrial = fmin(cflTrial, 2.0*cfl);
                    cflTrial = fmax(cflTrial, 0.1*cfl);
                    cfl = cflTrial;
                    cfl = fmin(cflTrial, cfl_max);
                }
            }
        } else { // user defined CFL growth
            if (cfl_schedule_iter_list.canFind(step)) {
                cfl = cfl_schedule_value_list[cfl_schedule_current_index];
                cfl_schedule_current_index += 1;
            }
        }

        // Update dt based on new CFL
        dt = determine_dt(cfl);
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        }
                
        if (step == nStartUpSteps) {
            // At the swap-over point from start-up phase to main phase
            // we need to do a few special things.
            // 1. Reset dt to user's choice for this new phase based on cfl1.
            if (GlobalConfig.is_master_task) { writefln("step= %d dt= %e  cfl1= %f", step, dt, cfl1); }
            cfl = cfl1;
            dt = determine_dt(cfl);
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
            if (GlobalConfig.is_master_task) { writefln("after choosing new timestep: %e", dt); }
            // 2. Reset the inexact Newton phase.
            //    We'll take some constant timesteps at the new dt
            //    until the residuals have dropped.
            inexactNewtonPhase = false;
        }

        if (step > nStartUpSteps) {
            // Adjust eta1 according to eta update strategy
            final switch (etaStrategy) {
            case EtaStrategy.constant:
                break;
            case EtaStrategy.geometric:
                eta1 *= etaRatioPerStep;
                // but don't go past eta1_min
                eta1 = fmax(eta1, eta1_min);
                break;
            case EtaStrategy.adaptive, EtaStrategy.adaptive_capped:
                // Use Eisenstat & Walker adaptive strategy no. 2
                auto etaOld = eta1;
                eta1 = gamma*pow(1./normRatio, alpha);
                // Apply Eisenstat & Walker safeguards
                auto testVal = gamma*pow(etaOld, alpha);
                if ( testVal > 0.1 ) {
                    eta1 = fmax(eta1, testVal);
                }
                // and don't let eta1 get larger than a max value
                eta1 = fmin(eta1, eta1_max);
                if ( etaStrategy == EtaStrategy.adaptive_capped ) {
                    // Additionally, we cap the maximum so that we never
                    // retreat on the eta1 value.
                    eta1_max = fmin(eta1, eta1_max);
                    // but don't let our eta1_max cap get too tiny
                    eta1_max = fmax(eta1, eta1_min);
                }
                break;
            }
        }

        normOld = normNew;
    }
    
    // for some simulations we freeze the limiter to assist in reducing the relative residuals
    // to machine precision - in these instances it is typically necessary to store the limiter 
    // values for further analysis.
    if (GlobalConfig.frozen_limiter) {
        string limValDir = "limiter-values";
        if (GlobalConfig.is_master_task) { ensure_directory_is_present(limValDir); }
        version(mpi_parallel) {
            MPI_Barrier(MPI_COMM_WORLD);
        }
	foreach (blk; localFluidBlocks) {
            string fileName = format("%s/limiter-values-b%04d.dat", limValDir, blk.id);
            auto outFile = File(fileName, "w");
	    foreach (cell; blk.cells) {
		outFile.writef("%.16e \n", cell.gradients.rhoPhi.re);
		outFile.writef("%.16e \n", cell.gradients.velxPhi.re);
		outFile.writef("%.16e \n", cell.gradients.velyPhi.re);
		if (blk.myConfig.dimensions == 3) {
		    outFile.writef("%.16e \n", cell.gradients.velzPhi.re);
		}
		outFile.writef("%.16e \n", cell.gradients.pPhi.re);
                foreach(it; 0 .. blk.myConfig.turb_model.nturb){
		    outFile.writef("%.16e \n", cell.gradients.turbPhi[it].re);
		}
                version(multi_species_gas){
                if ( blk.myConfig.gmodel.n_species > 1 ) {
                    foreach(sp; 0 .. blk.myConfig.gmodel.n_species) {
                        outFile.writef("%.16e \n", cell.gradients.massfPhi[sp].re);
                    }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. blk.myConfig.gmodel.n_modes) {
                    outFile.writef("%.16e \n", cell.gradients.T_modesPhi[imode].re);
                }
                foreach(imode; 0 .. blk.myConfig.gmodel.n_modes) {
                    outFile.writef("%.16e \n", cell.gradients.u_modesPhi[imode].re);
                }
                }

                
            }
            outFile.close();
        }
    } // end if (GlobalConfig.frozen_limiter)    
}

void allocate_global_fluid_workspace()
{
    size_t mOuter = to!size_t(GlobalConfig.sssOptions.maxOuterIterations);
    g0.length = mOuter+1;
    g1.length = mOuter+1;
    h.length = mOuter+1;
    hR.length = mOuter+1;
    H0 = new Matrix!number(mOuter+1, mOuter);
    H1 = new Matrix!number(mOuter+1, mOuter);
    Gamma = new Matrix!number(mOuter+1, mOuter+1);
    Q0 = new Matrix!number(mOuter+1, mOuter+1);
    Q1 = new Matrix!number(mOuter+1, mOuter+1);
}

double determine_dt(double cflInit)
{
    double signal, dt;
    bool first = true;
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            signal = cell.signal_frequency();
            if (blk.myConfig.sssOptions.inviscidCFL) {
                // overwrite signal with inviscid signal
                signal = cell.signal_hyp.re;
            }
            cell.dt_local = cflInit / signal;
            if (first) {
                dt = cell.dt_local;
                first = false;
            }
            else {
                dt = fmin(dt, cell.dt_local);
            }
        }
    }
    return dt;
}

double determine_min_cfl(double dt)
{
    double signal, cfl_local, cfl;
    bool first = true;

    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            signal = cell.signal_frequency();
            cfl_local = dt * signal;
            if (first) {
                cfl = cfl_local;
                first = false;
            }
            else {
                cfl = fmin(cfl, cfl_local);
            }
        }
    }
    return cfl;
}

void evalRHS(double pseudoSimTime, int ftl)
{
    fnCount++;
    int gtl = 0;

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
    }
    exchange_ghost_cell_boundary_data(pseudoSimTime, 0, ftl);
    foreach (blk; localFluidBlocks) {
        blk.applyPreReconAction(pseudoSimTime, 0, ftl);
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
                boundary.preSpatialDerivActionAtBndryFaces[0].apply(pseudoSimTime, gtl, ftl);
            }
        }
    }

    bool allow_high_order_interpolation = true;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase0(allow_high_order_interpolation, gtl);
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(pseudoSimTime, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase1(allow_high_order_interpolation, gtl);
    }
    foreach (blk; localFluidBlocks) {
        blk.applyPostConvFluxAction(pseudoSimTime, 0, ftl);
    }

    if (GlobalConfig.viscous) {
        foreach (blk; localFluidBlocks) {
            blk.applyPreSpatialDerivActionAtBndryFaces(pseudoSimTime, 0, ftl);
            blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.flow_property_spatial_derivatives(0);
        }
        // for unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
        // we need to transfer the viscous gradients before the flux calc
        exchange_ghost_cell_boundary_viscous_gradient_data(pseudoSimTime, to!int(gtl), to!int(ftl));
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
            blk.applyPostDiffFluxAction(pseudoSimTime, 0, ftl);
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
                getUDFSourceTermsForCell(blk.myL, cell, 0, pseudoSimTime, blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                cell.add_udf_source_vector();
            }
            cell.time_derivatives(0, ftl);
        }
    }
}

void evalJacobianVecProd(double pseudoSimTime, double sigma, int LHSeval, int RHSeval)
{
    if (GlobalConfig.sssOptions.useComplexMatVecEval)
        evalComplexMatVecProd(pseudoSimTime, sigma, LHSeval, RHSeval);
    else
        evalRealMatVecProd(pseudoSimTime, sigma, LHSeval, RHSeval);
}

void evalRealMatVecProd(double pseudoSimTime, double sigma, int LHSeval, int RHSeval)
{
    foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(LHSeval); }
    bool frozen_limiter_save = GlobalConfig.frozen_limiter;
    if (GlobalConfig.sssOptions.frozenLimiterOnLHS) {
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = true; }
    }
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;
    size_t SPECIES = GlobalConfig.cqi.species;
    size_t MODES = GlobalConfig.cqi.modes;

    // We perform a Frechet derivative to evaluate J*D^(-1)v
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        size_t nsp = blk.myConfig.gmodel.n_species;
        size_t nmodes = blk.myConfig.gmodel.n_modes;
        auto cqi = blk.myConfig.cqi;
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
        int cellCount = 0;
        foreach (cell; blk.cells) {
            cell.U[1].copy_values_from(cell.U[0]);
            cell.U[1].vec[cqi.mass] += sigma*blk.zed[cellCount+MASS];
            cell.U[1].vec[cqi.xMom] += sigma*blk.zed[cellCount+X_MOM];
            cell.U[1].vec[cqi.yMom] += sigma*blk.zed[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                cell.U[1].vec[cqi.zMom] += sigma*blk.zed[cellCount+Z_MOM];
            cell.U[1].vec[cqi.totEnergy] += sigma*blk.zed[cellCount+TOT_ENERGY];
            foreach(it; 0 .. nturb){
                cell.U[1].vec[cqi.rhoturb+it] += sigma*blk.zed[cellCount+TKE+it];
            }
            version(multi_species_gas){
            if ( nsp > 1 ) {
                foreach(sp; 0 .. nsp) { cell.U[1].vec[cqi.species+sp] += sigma*blk.zed[cellCount+SPECIES+sp]; }
            } else {
                // enforce mass fraction of 1 for single species gas
                if (blk.myConfig.n_species == 1) {
                    cell.U[1].vec[cqi.species+0] = cell.U[1].vec[cqi.mass];
                }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes) { cell.U[1].vec[cqi.modes+imode] += sigma*blk.zed[cellCount+MODES+imode]; }
            }
            cell.decode_conserved(0, 1, 0.0);
            cellCount += nConserved;
        }
    }
    evalRHS(pseudoSimTime, 1);
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        size_t nsp = blk.myConfig.gmodel.n_species;
        size_t nmodes = blk.myConfig.gmodel.n_modes;
        auto cqi = blk.myConfig.cqi;
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.zed[cellCount+MASS] = (cell.dUdt[1].vec[cqi.mass] - blk.FU[cellCount+MASS])/(sigma);
            blk.zed[cellCount+X_MOM] = (cell.dUdt[1].vec[cqi.xMom] - blk.FU[cellCount+X_MOM])/(sigma);
            blk.zed[cellCount+Y_MOM] = (cell.dUdt[1].vec[cqi.yMom] - blk.FU[cellCount+Y_MOM])/(sigma);
            if ( blk.myConfig.dimensions == 3 )
                blk.zed[cellCount+Z_MOM] = (cell.dUdt[1].vec[cqi.zMom] - blk.FU[cellCount+Z_MOM])/(sigma);
            blk.zed[cellCount+TOT_ENERGY] = (cell.dUdt[1].vec[cqi.totEnergy] - blk.FU[cellCount+TOT_ENERGY])/(sigma);
            foreach(it; 0 .. nturb){
                blk.zed[cellCount+TKE+it] = (cell.dUdt[1].vec[cqi.rhoturb+it] - blk.FU[cellCount+TKE+it])/(sigma);
            }
            version(multi_species_gas){
            if ( nsp > 1 ) {
                foreach(sp; 0 .. nsp){ blk.zed[cellCount+SPECIES+sp] = (cell.dUdt[1].vec[cqi.species+sp] - blk.FU[cellCount+SPECIES+sp])/(sigma); }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes){ blk.zed[cellCount+MODES+imode] = (cell.dUdt[1].vec[cqi.modes+imode] - blk.FU[cellCount+MODES+imode])/(sigma); }
            }
            cell.decode_conserved(0, 0, 0.0);
            cellCount += nConserved;
        }
    }
    foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(RHSeval); }
    if (GlobalConfig.sssOptions.frozenLimiterOnLHS) {
        foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = frozen_limiter_save; }
    }
}

void evalComplexMatVecProd(double pseudoSimTime, double sigma, int LHSeval, int RHSeval)
{
    version(complex_numbers) {
        foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(LHSeval); }
        bool frozen_limiter_save = GlobalConfig.frozen_limiter;
        if (GlobalConfig.sssOptions.frozenLimiterOnLHS) {
            foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = true; }
        }
        // Make a stack-local copy of conserved quantities info
        size_t nConserved = GlobalConfig.cqi.n;
        size_t MASS = GlobalConfig.cqi.mass;
        size_t X_MOM = GlobalConfig.cqi.xMom;
        size_t Y_MOM = GlobalConfig.cqi.yMom;
        size_t Z_MOM = GlobalConfig.cqi.zMom;
        size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
        size_t TKE = GlobalConfig.cqi.rhoturb;
        size_t SPECIES = GlobalConfig.cqi.species;
        size_t MODES = GlobalConfig.cqi.modes;

        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.clear_fluxes_of_conserved_quantities();
            foreach (cell; blk.cells) cell.clear_source_vector();

            size_t nturb = blk.myConfig.turb_model.nturb;
            size_t nsp = blk.myConfig.gmodel.n_species;
            size_t nmodes = blk.myConfig.gmodel.n_modes;
            auto cqi = blk.myConfig.cqi;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].vec[cqi.mass] += complex(0.0, sigma*blk.zed[cellCount+MASS].re);
                cell.U[1].vec[cqi.xMom] += complex(0.0, sigma*blk.zed[cellCount+X_MOM].re);
                cell.U[1].vec[cqi.yMom] += complex(0.0, sigma*blk.zed[cellCount+Y_MOM].re);
                if ( blk.myConfig.dimensions == 3 )
                    cell.U[1].vec[cqi.zMom] += complex(0.0, sigma*blk.zed[cellCount+Z_MOM].re);
                cell.U[1].vec[cqi.totEnergy] += complex(0.0, sigma*blk.zed[cellCount+TOT_ENERGY].re);
                foreach(it; 0 .. nturb){
                    cell.U[1].vec[cqi.rhoturb+it] += complex(0.0, sigma*blk.zed[cellCount+TKE+it].re);
                }
                version(multi_species_gas){
                if ( nsp > 1 ) {
                    foreach(sp; 0 .. nsp){ cell.U[1].vec[cqi.species+sp] += complex(0.0, sigma*blk.zed[cellCount+SPECIES+sp].re); }
                } else {
                    // enforce mass fraction of 1 for single species gas
                    cell.U[1].vec[cqi.species+0] = cell.U[1].vec[cqi.mass];
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. nmodes){ cell.U[1].vec[cqi.modes+imode] += complex(0.0, sigma*blk.zed[cellCount+MODES+imode].re); }
                }
                cell.decode_conserved(0, 1, 0.0);
                cellCount += nConserved;
            }
        }
        evalRHS(pseudoSimTime, 1);
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            size_t nsp = blk.myConfig.gmodel.n_species;
            size_t nmodes = blk.myConfig.gmodel.n_modes;
            auto cqi = blk.myConfig.cqi;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                blk.zed[cellCount+MASS] = cell.dUdt[1].vec[cqi.mass].im/(sigma);
                blk.zed[cellCount+X_MOM] = cell.dUdt[1].vec[cqi.xMom].im/(sigma);
                blk.zed[cellCount+Y_MOM] = cell.dUdt[1].vec[cqi.yMom].im/(sigma);
                if ( blk.myConfig.dimensions == 3 )
                    blk.zed[cellCount+Z_MOM] = cell.dUdt[1].vec[cqi.zMom].im/(sigma);
                blk.zed[cellCount+TOT_ENERGY] = cell.dUdt[1].vec[cqi.totEnergy].im/(sigma);
                foreach(it; 0 .. nturb){
                    blk.zed[cellCount+TKE+it] = cell.dUdt[1].vec[cqi.rhoturb+it].im/(sigma);
                }
                version(multi_species_gas){
                if ( nsp > 1 ) {
                    foreach(sp; 0 .. nsp){ blk.zed[cellCount+SPECIES+sp] = cell.dUdt[1].vec[cqi.species+sp].im/(sigma); }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. nmodes){ blk.zed[cellCount+MODES+imode] = cell.dUdt[1].vec[cqi.modes+imode].im/(sigma); }
                }
                cellCount += nConserved;
            }
            // we must explicitly remove the imaginary components from the cell and interface flowstates
            foreach(cell; blk.cells) { cell.fs.clear_imaginary_components(); }
            foreach(bc; blk.bc) {
                foreach(ghostcell; bc.ghostcells) { ghostcell.fs.clear_imaginary_components(); }
            }
            foreach(face; blk.faces) { face.fs.clear_imaginary_components(); }
        }
        foreach (blk; parallel(localFluidBlocks,1)) { blk.set_interpolation_order(RHSeval); }
        if (GlobalConfig.sssOptions.frozenLimiterOnLHS) {
            foreach (blk; parallel(localFluidBlocks,1)) { GlobalConfig.frozen_limiter = frozen_limiter_save; }
        }
    } else {
        throw new Error("Oops. Steady-State Solver setting: useComplexMatVecEval is not compatible with real-number version of the code.");
    }
}


string dot_over_blocks(string dot, string A, string B)
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

string norm2_over_blocks(string norm2, string blkMember)
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

void lusgs_solve(int step, double pseudoSimTime, double dt, ref double residual, int startStep)
{
    // Data-parallel implementation of a matrix-free LU-SGS method
    //
    // Wright (1997),
    // A family of data-parallel relaxation methods for the Navier-Stokes equations,
    // University of Minesota
    //
    // Sharov, Luo, and Baum (2000),
    // Implementation of Unstructured Grid GMRES + LU-SGS Method on Shared-Memory, Cache-Based Parallel Computer,
    // AIAA conference paper
    //

    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;

    // Evaluate RHS residual (R)
    evalRHS(pseudoSimTime, 0);

    // calculate global residual
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        auto cqi = blk.myConfig.cqi;
        int cellCount = 0;
        foreach (i, cell; blk.cells) {
            blk.FU[cellCount+MASS] = cell.dUdt[0].vec[cqi.mass];
            blk.FU[cellCount+X_MOM] = cell.dUdt[0].vec[cqi.xMom];
            blk.FU[cellCount+Y_MOM] = cell.dUdt[0].vec[cqi.yMom];
            if ( GlobalConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = cell.dUdt[0].vec[cqi.zMom];
            blk.FU[cellCount+TOT_ENERGY] = cell.dUdt[0].vec[cqi.totEnergy];
            foreach(it; 0 .. nturb) {
                blk.FU[cellCount+TKE+it] = cell.dUdt[0].vec[cqi.rhoturb+it];
            }
            cellCount += nConserved;
        }
    }
    mixin(dot_over_blocks("residual", "FU", "FU"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(residual), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    residual = sqrt(residual);

    mixin(lusgs_solve("dU", "FU"));
} // end lusgs_solve()

string lusgs_solve(string lhs_vec, string rhs_vec)
{
    string code = "{

    int kmax = GlobalConfig.sssOptions.maxSubIterations;
    bool matrix_based = false;
    double omega = 2.0;

    // 2. initial subiteration
    foreach (blk; parallel(localFluidBlocks,1)) {
        int cellCount = 0;
        foreach (cell; blk.cells) {
            number dtInv;
            if (blk.myConfig.with_local_time_stepping) { dtInv = 1.0/cell.dt_local; }
            else { dtInv = 1.0/dt; }
            auto z_local = blk."~lhs_vec~"[cellCount..cellCount+nConserved]; // this is actually a reference not a copy
            cell.lusgs_startup_iteration(dtInv, omega, z_local, blk."~rhs_vec~"[cellCount..cellCount+nConserved]);
            cellCount += nConserved;
        }
    }

    // 3. kmax subiterations
    foreach (k; 0 .. kmax) {
         // shuffle dU values
         foreach (blk; parallel(localFluidBlocks,1)) {
             int cellCount = 0;
             foreach (cell; blk.cells) {
                 cell.dUk[0..nConserved] = blk."~lhs_vec~"[cellCount..cellCount+nConserved];
                 cellCount += nConserved;
             }
         }

         // exchange boundary dU values
         exchange_ghost_cell_boundary_data(pseudoSimTime, 0, 0);

         // perform subiteraion
         foreach (blk; parallel(localFluidBlocks,1)) {
             int cellCount = 0;
             foreach (cell; blk.cells) {
                  auto z_local = blk."~lhs_vec~"[cellCount..cellCount+nConserved]; // this is actually a reference not a copy
                  cell.lusgs_relaxation_iteration(omega, matrix_based, z_local, blk."~rhs_vec~"[cellCount..cellCount+nConserved]);
                  cellCount += nConserved;
             }
         }
    }

    }";
    return code;
}

void rpcGMRES_solve(int step, double pseudoSimTime, double dt, double eta, double sigma, bool usePreconditioner,
                    ref double residual, ref int nRestarts, ref int nIters, ref double linSolResid, int startStep, int LHSeval, int RHSeval, ref bool pc_matrix_evaluated)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;
    size_t SPECIES = GlobalConfig.cqi.species;
    size_t MODES = GlobalConfig.cqi.modes;

    number resid;

    int interpOrderSave = GlobalConfig.interpolation_order;
    // Presently, just do one block
    int maxIters = GlobalConfig.sssOptions.maxOuterIterations;
    // We add 1 because the user thinks of "re"starts, so they
    // might legitimately ask for no restarts. We still have
    // to execute at least once.
    int maxRestarts = GlobalConfig.sssOptions.maxRestarts + 1;
    size_t m = to!size_t(maxIters);
    size_t r;
    size_t iterCount;

    // Variables for max rates of change
    // Use these for equation scaling.
    double minNonDimVal = 1.0; // minimum value used for non-dimensionalisation
                               // when our time rates of change are very small
                               // then we'll avoid non-dimensionalising by
                               // values close to zero.

    // 1. Evaluate r0, beta, v1
    evalRHS(pseudoSimTime, 0);

    // Store dUdt[0] as F(U)
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        size_t nsp = blk.myConfig.gmodel.n_species;
        size_t nmodes = blk.myConfig.gmodel.n_modes;
        auto cqi = blk.myConfig.cqi;
        int cellCount = 0;
        blk.maxRate.vec[cqi.mass] = 0.0;
        blk.maxRate.vec[cqi.xMom] = 0.0;
        blk.maxRate.vec[cqi.yMom] = 0.0;
        if ( blk.myConfig.dimensions == 3 )
            blk.maxRate.vec[cqi.zMom] = 0.0;
        blk.maxRate.vec[cqi.totEnergy] = 0.0;
        foreach(it; 0 .. nturb){
            blk.maxRate.vec[cqi.rhoturb+it] = 0.0;
        }
        version(multi_species_gas){
        if ( nsp > 1 ) {
            foreach(sp; 0 .. nsp){ blk.maxRate.vec[cqi.species+sp] = 0.0; }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. nmodes){ blk.maxRate.vec[cqi.modes+imode] = 0.0; }
        }
        foreach (i, cell; blk.cells) {
            blk.FU[cellCount+MASS] = cell.dUdt[0].vec[cqi.mass];
            blk.FU[cellCount+X_MOM] = cell.dUdt[0].vec[cqi.xMom];
            blk.FU[cellCount+Y_MOM] = cell.dUdt[0].vec[cqi.yMom];
            if ( blk.myConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = cell.dUdt[0].vec[cqi.zMom];
            blk.FU[cellCount+TOT_ENERGY] = cell.dUdt[0].vec[cqi.totEnergy];
            foreach(it; 0 .. nturb){
                blk.FU[cellCount+TKE+it] = cell.dUdt[0].vec[cqi.rhoturb+it];
            }
            version(multi_species_gas){
            if ( nsp > 1 ) {
                foreach(sp; 0 .. nsp){ blk.FU[cellCount+SPECIES+sp] = cell.dUdt[0].vec[cqi.species+sp]; }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes){ blk.FU[cellCount+MODES+imode] = cell.dUdt[0].vec[cqi.modes+imode]; }
            }
            cellCount += nConserved;
            /*
            if (blk.id == 0) {
                writefln("i= %d, dUdt.mass= %e", i, cell.dUdt[0].vec[cqi.mass].re);
            }
            */
            blk.maxRate.vec[cqi.mass] = fmax(blk.maxRate.vec[cqi.mass], fabs(cell.dUdt[0].vec[cqi.mass]));
            /*
            if (blk.id == 0) {
                writefln("i= %d, maxRate.vec[cqi.mass]= %e", i, blk.maxRate.vec[cqi.mass].re);
            }
            */
            blk.maxRate.vec[cqi.xMom] = fmax(blk.maxRate.vec[cqi.xMom], fabs(cell.dUdt[0].vec[cqi.xMom]));
            blk.maxRate.vec[cqi.yMom] = fmax(blk.maxRate.vec[cqi.yMom], fabs(cell.dUdt[0].vec[cqi.yMom]));
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.vec[cqi.zMom] = fmax(blk.maxRate.vec[cqi.zMom], fabs(cell.dUdt[0].vec[cqi.zMom]));
            blk.maxRate.vec[cqi.totEnergy] = fmax(blk.maxRate.vec[cqi.totEnergy], fabs(cell.dUdt[0].vec[cqi.totEnergy]));
            foreach(it; 0 .. nturb){
                blk.maxRate.vec[cqi.rhoturb+it] = fmax(blk.maxRate.vec[cqi.rhoturb+it], fabs(cell.dUdt[0].vec[cqi.rhoturb+it]));
            }
            version(multi_species_gas){
            if ( nsp > 1 ) {
                foreach(sp; 0 .. nsp){ blk.maxRate.vec[cqi.species+sp] = fmax(blk.maxRate.vec[cqi.species+sp], fabs(cell.dUdt[0].vec[cqi.species+sp])); }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes){ blk.maxRate.vec[cqi.modes+imode] = fmax(blk.maxRate.vec[cqi.modes+imode], fabs(cell.dUdt[0].vec[cqi.modes+imode])); }
            }

        }
    }

    number maxMass = 0.0;
    number maxMomX = 0.0;
    number maxMomY = 0.0;
    number maxMomZ = 0.0;
    number maxEnergy = 0.0;
    number[2] maxTurb; maxTurb[0] = 0.0; maxTurb[1] = 0.0;
    // we currently expect no more than 32 species
    number[32] maxSpecies;
    foreach (sp; 0 .. maxSpecies.length) { maxSpecies[sp] = 0.0; }
    number[2] maxModes;
    foreach (imode; 0 .. maxModes.length) { maxModes[imode] = 0.0; }

    foreach (blk; localFluidBlocks) {
        auto cqi = blk.myConfig.cqi;
        maxMass = fmax(maxMass, blk.maxRate.vec[cqi.mass]);
        maxMomX = fmax(maxMomX, blk.maxRate.vec[cqi.xMom]);
        maxMomY = fmax(maxMomY, blk.maxRate.vec[cqi.yMom]);
        if ( blk.myConfig.dimensions == 3 )
            maxMomZ = fmax(maxMomZ, blk.maxRate.vec[cqi.zMom]);
        maxEnergy = fmax(maxEnergy, blk.maxRate.vec[cqi.totEnergy]);
        foreach(it; 0 .. blk.myConfig.turb_model.nturb){
            maxTurb[it] = fmax(maxTurb[it], blk.maxRate.vec[cqi.rhoturb+it]);
        }
        version(multi_species_gas){
        if ( blk.myConfig.gmodel.n_species > 1 ) {
            foreach(sp; 0 .. blk.myConfig.gmodel.n_species){
                maxSpecies[sp] = fmax(maxSpecies[sp], blk.maxRate.vec[cqi.species+sp]);
            }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. blk.myConfig.gmodel.n_modes){
            maxModes[imode] = fmax(maxModes[imode], blk.maxRate.vec[cqi.modes+imode]);
        }
        }
    }
    // In distributed memory, reduce the max values and ensure everyone has a copy
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(maxMass.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &(maxMomX.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &(maxMomY.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (GlobalConfig.dimensions == 3) { MPI_Allreduce(MPI_IN_PLACE, &(maxMomZ.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); }
        MPI_Allreduce(MPI_IN_PLACE, &(maxEnergy.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        foreach (it; 0 .. GlobalConfig.turb_model.nturb) {
            MPI_Allreduce(MPI_IN_PLACE, &(maxTurb[it].re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        if ( GlobalConfig.gmodel_master.n_species > 1 ) {
            foreach (sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                MPI_Allreduce(MPI_IN_PLACE, &(maxSpecies[sp].re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            }
        }
        if ( GlobalConfig.gmodel_master.n_modes > 1 ) {
            foreach (imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                MPI_Allreduce(MPI_IN_PLACE, &(maxModes[imode].re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            }
        }
    }
    // Place some guards when time-rate-of-changes are very small.
    maxMass = fmax(maxMass, minNonDimVal);
    maxMomX = fmax(maxMomX, minNonDimVal);
    maxMomY = fmax(maxMomY, minNonDimVal);
    if ( GlobalConfig.dimensions == 3 )
        maxMomZ = fmax(maxMomZ, minNonDimVal);
    maxEnergy = fmax(maxEnergy, minNonDimVal);
    number maxMom = max(maxMomX, maxMomY, maxMomZ);
    foreach(it; 0 .. GlobalConfig.turb_model.nturb){
        maxTurb[it] = fmax(maxTurb[it], minNonDimVal);
    }
    foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
        maxSpecies[sp] = fmax(maxSpecies[sp], minNonDimVal);
    }
    foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
        maxModes[imode] = fmax(maxModes[imode], minNonDimVal);
    }

    // Get a copy of the maxes out to each block
    foreach (blk; parallel(localFluidBlocks,1)) {
        auto cqi = blk.myConfig.cqi;
        if (blk.myConfig.sssOptions.useScaling) {
            blk.maxRate.vec[cqi.mass] = maxMass;
            blk.maxRate.vec[cqi.xMom] = maxMomX;
            blk.maxRate.vec[cqi.yMom] = maxMomY;
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.vec[cqi.zMom] = maxMomZ;
            blk.maxRate.vec[cqi.totEnergy] = maxEnergy;
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                blk.maxRate.vec[cqi.rhoturb+it] = maxTurb[it];
            }
            version(multi_species_gas){
            if ( blk.myConfig.gmodel.n_species > 1 ) {
                foreach(sp; 0 .. blk.myConfig.gmodel.n_species){
                    blk.maxRate.vec[cqi.species+sp] = maxSpecies[sp];
                }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. blk.myConfig.gmodel.n_modes){
                blk.maxRate.vec[cqi.modes+imode] = maxModes[imode];
            }
            }
        }
        else { // just scale by 1
            blk.maxRate.vec[cqi.mass] = 1.0;
            blk.maxRate.vec[cqi.xMom] = 1.0;
            blk.maxRate.vec[cqi.yMom] = 1.0;
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.vec[cqi.zMom] = 1.0;
            blk.maxRate.vec[cqi.totEnergy] = 1.0;
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                blk.maxRate.vec[cqi.rhoturb+it] = 1.0;
            }
            version(multi_species_gas){
            if ( blk.myConfig.gmodel.n_species ) {
                foreach(sp; 0 .. blk.myConfig.gmodel.n_species){
                    blk.maxRate.vec[cqi.species+sp] = 1.0;
                }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. blk.myConfig.gmodel.n_modes){
                blk.maxRate.vec[cqi.modes+imode] = 1.0;
            }
            }

        }
    }

    // scale turbulence variables before calculating norm
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            foreach(it; 0 .. nturb){
                number fac = blk.myConfig.turb_model.gmres_scaling_factor(it);
                blk.FU[cellCount+TKE+it] = 1.0/fac*blk.FU[cellCount+TKE+it];
            }
            cellCount += nConserved;
        }
    }

    double unscaledNorm2;
    mixin(dot_over_blocks("unscaledNorm2", "FU", "FU"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &unscaledNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    unscaledNorm2 = sqrt(unscaledNorm2);

    // unscale turbulence variables after calculating norm
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            foreach(it; 0 .. nturb){
                number fac = blk.myConfig.turb_model.gmres_scaling_factor(it);
                blk.FU[cellCount+TKE+it] = fac*blk.FU[cellCount+TKE+it];
            }
            cellCount += nConserved;
        }
    }

    // Initialise some arrays and matrices that have already been allocated
    g0[] = to!number(0.0);
    g1[] = to!number(0.0);
    H0.zeros();
    H1.zeros();

    // We'll scale r0 against these max rates of change.
    // r0 = b - A*x0
    // Taking x0 = [0] (as is common) gives r0 = b = FU
    // apply scaling
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        size_t nsp = blk.myConfig.gmodel.n_species;
        size_t nmodes = blk.myConfig.gmodel.n_modes;
        auto cqi = blk.myConfig.cqi;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.r0[cellCount+MASS] = (1./blk.maxRate.vec[cqi.mass])*blk.FU[cellCount+MASS];
            blk.r0[cellCount+X_MOM] = (1./blk.maxRate.vec[cqi.xMom])*blk.FU[cellCount+X_MOM];
            blk.r0[cellCount+Y_MOM] = (1./blk.maxRate.vec[cqi.yMom])*blk.FU[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                blk.r0[cellCount+Z_MOM] = (1./blk.maxRate.vec[cqi.zMom])*blk.FU[cellCount+Z_MOM];
            blk.r0[cellCount+TOT_ENERGY] = (1./blk.maxRate.vec[cqi.totEnergy])*blk.FU[cellCount+TOT_ENERGY];
            foreach(it; 0 .. nturb){
                blk.r0[cellCount+TKE+it] = (1./blk.maxRate.vec[cqi.rhoturb+it])*blk.FU[cellCount+TKE+it];
            }
            version(multi_species_gas){
            if ( nsp > 1 ) {
                foreach(sp; 0 .. nsp){
                    blk.r0[cellCount+SPECIES+sp] = (1./blk.maxRate.vec[cqi.species+sp])*blk.FU[cellCount+SPECIES+sp];
                }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes){
                blk.r0[cellCount+MODES+imode] = (1./blk.maxRate.vec[cqi.modes+imode])*blk.FU[cellCount+MODES+imode];
            }
            }

            cellCount += nConserved;
        }
    }

    // Then compute v = r0/||r0||
    number beta;
    mixin(dot_over_blocks("beta", "r0", "r0"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
    }
    beta = sqrt(beta);
    number beta0 = beta;
    g0[0] = beta;
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) {
            blk.v[k] = blk.r0[k]/beta;
            blk.V[k,0] = blk.v[k];
        }
    }

    // Compute tolerance
    auto outerTol = eta*beta;

    // Compute precondition matrix
    // if ( usePreconditioner && ( step == startStep || nIters > GlobalConfig.sssOptions.frozenPreconditionerCount) ) {
    pc_matrix_evaluated = false;
    if (usePreconditioner && ( (m == nIters && GlobalConfig.sssOptions.useAdaptivePreconditioner) ||
                               (step == startStep) ||
                               (step%GlobalConfig.sssOptions.frozenPreconditionerCount == 0) )) {
        pc_matrix_evaluated = true;
        final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
        case PreconditionMatrixType.jacobi:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.evaluate_jacobian();
                blk.flowJacobian.prepare_jacobi_preconditioner(blk.cells, dt, blk.cells.length, nConserved);
            }
            break;
        case PreconditionMatrixType.ilu:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.evaluate_jacobian();
                blk.flowJacobian.prepare_ilu_preconditioner(blk.cells, dt, blk.cells.length, nConserved);
            }
            break;
        case PreconditionMatrixType.sgs:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.evaluate_jacobian();
                blk.flowJacobian.prepare_sgs_preconditioner(blk.cells, dt, blk.cells.length, nConserved);
            }
            break;
        case PreconditionMatrixType.sgs_relax:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.evaluate_jacobian();
                blk.flowJacobian.prepare_sgsr_preconditioner(blk.cells, dt, blk.cells.length, nConserved);
            }
            break;
        case PreconditionMatrixType.lu_sgs:
            // do nothing
            break;
        }
    }

    // 2. Start outer-loop of restarted GMRES
    for ( r = 0; r < maxRestarts; r++ ) {
        // 2a. Begin iterations
        foreach (j; 0 .. m) {
            iterCount = j+1;

            // apply scaling
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.turb_model.nturb;
                size_t nsp = blk.myConfig.gmodel.n_species;
                size_t nmodes = blk.myConfig.gmodel.n_modes;
                auto cqi = blk.myConfig.cqi;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.v[cellCount+MASS] *= (blk.maxRate.vec[cqi.mass]);
                    blk.v[cellCount+X_MOM] *= (blk.maxRate.vec[cqi.xMom]);
                    blk.v[cellCount+Y_MOM] *= (blk.maxRate.vec[cqi.yMom]);
                    if ( blk.myConfig.dimensions == 3 )
                        blk.v[cellCount+Z_MOM] *= (blk.maxRate.vec[cqi.zMom]);
                    blk.v[cellCount+TOT_ENERGY] *= (blk.maxRate.vec[cqi.totEnergy]);
                    foreach(it; 0 .. nturb){
                        blk.v[cellCount+TKE+it] *= (blk.maxRate.vec[cqi.rhoturb+it]);
                    }
                    version(multi_species_gas){
                    if ( nsp > 1 ) {
                        foreach(sp; 0 .. nsp){ blk.v[cellCount+SPECIES+sp] *= (blk.maxRate.vec[cqi.species+sp]); }
                    }
                    }
                    version(multi_T_gas){
                    foreach(imode; 0 .. nmodes){ blk.v[cellCount+MODES+imode] *= (blk.maxRate.vec[cqi.modes+imode]); }
                    }
                    cellCount += nConserved;
                }
            }

            // apply preconditioning
            if (usePreconditioner && step >= GlobalConfig.sssOptions.startPreconditioning) {
                final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
                case PreconditionMatrixType.jacobi:
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        blk.flowJacobian.x[] = blk.v[];
                        nm.smla.multiply(blk.flowJacobian.local, blk.flowJacobian.x, blk.zed);
                    }
                    break;
                case PreconditionMatrixType.ilu:
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        blk.zed[] = blk.v[];
                        nm.smla.solve(blk.flowJacobian.local, blk.zed);
                    }
                    break;
                case PreconditionMatrixType.sgs:
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        blk.zed[] = blk.v[];
                        nm.smla.sgs(blk.flowJacobian.local, blk.flowJacobian.diagonal, blk.zed, to!int(nConserved), blk.flowJacobian.D, blk.flowJacobian.Dinv);
                    }
                    break;
                case PreconditionMatrixType.sgs_relax:
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        int local_kmax = GlobalConfig.sssOptions.maxSubIterations;
                        blk.zed[] = blk.v[];
                                nm.smla.sgsr(blk.flowJacobian.local, blk.zed, blk.flowJacobian.x, to!int(nConserved), local_kmax, blk.flowJacobian.Dinv);
                    }
                    break;
                case PreconditionMatrixType.lu_sgs:
                    mixin(lusgs_solve("zed", "v"));
                    break;
                } // end switch
            } else {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    blk.zed[] = blk.v[];
                }
            }

            if (!GlobalConfig.sssOptions.useComplexMatVecEval) {
                // calculate sigma without scaling
                number sumv = 0.0;
                mixin(dot_over_blocks("sumv", "zed", "zed"));
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &(sumv.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(sumv.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
                }

                auto eps0 = GlobalConfig.sssOptions.sigma1;
                number N = 0.0;
                number sume = 0.0;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    int cellCount = 0;
                    foreach (cell; blk.cells) {
                        foreach (val; cell.U[0].vec) {
                            sume += eps0*abs(val) + eps0;
                            N += 1;
                        }
                        cellCount += nConserved;
                    }
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &(sume.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(sume.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &(N.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(N.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
                }

                sigma = (sume/(N*sqrt(sumv))).re;
                //writeln("sigma : ", sigma, ", ", sume, ", ", sumv, ", ", N);
            }

            // Prepare 'w' with (I/dt)(P^-1)v term;
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (i, cell; blk.cells) {
                    foreach (k; 0..nConserved) {
                        ulong idx = i*nConserved + k;
                        number dtInv;
                        if (GlobalConfig.with_local_time_stepping) { dtInv = 1.0/cell.dt_local; }
                        else { dtInv = 1.0/dt; }
                        blk.w[idx] = dtInv*blk.zed[idx];
                    }
                }
            }


            // Evaluate Jz and place in z
            evalJacobianVecProd(pseudoSimTime, sigma, LHSeval, RHSeval);

            // Now we can complete calculation of w
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars)  blk.w[k] = blk.w[k] - blk.zed[k];
            }

            // apply scaling
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.turb_model.nturb;
                size_t nsp = blk.myConfig.gmodel.n_species;
                size_t nmodes = blk.myConfig.gmodel.n_modes;
                auto cqi = blk.myConfig.cqi;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.w[cellCount+MASS] *= (1./blk.maxRate.vec[cqi.mass]);
                    blk.w[cellCount+X_MOM] *= (1./blk.maxRate.vec[cqi.xMom]);
                    blk.w[cellCount+Y_MOM] *= (1./blk.maxRate.vec[cqi.yMom]);
                    if ( blk.myConfig.dimensions == 3 )
                        blk.w[cellCount+Z_MOM] *= (1./blk.maxRate.vec[cqi.zMom]);
                    blk.w[cellCount+TOT_ENERGY] *= (1./blk.maxRate.vec[cqi.totEnergy]);
                    foreach(it; 0 .. nturb){
                        blk.w[cellCount+TKE+it] *= (1./blk.maxRate.vec[cqi.rhoturb+it]);
                    }
                    version(multi_species_gas){
                    if ( nsp > 1 ) {
                        foreach(sp; 0 .. nsp){ blk.w[cellCount+SPECIES+sp] *= (1./blk.maxRate.vec[cqi.species+sp]); }
                    }
                    }
                    version(multi_T_gas){
                    foreach(imode; 0 .. nmodes){ blk.w[cellCount+MODES+imode] *= (1./blk.maxRate.vec[cqi.modes+imode]); }
                    }
                    cellCount += nConserved;
                }
            }


            // The remainder of the algorithm looks a lot like any standard
            // GMRES implementation (for example, see smla.d)
            foreach (i; 0 .. j+1) {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    // Extract column 'i'
                    foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.V[k,i];
                }
                number H0_ij;
                mixin(dot_over_blocks("H0_ij", "w", "v"));
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
            mixin(dot_over_blocks("H0_jp1j", "w", "w"));
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
            if ( j != 0 ) {
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
            if ( j == 0 ) {
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
            resid = fabs(g1[j+1]);
            // DEBUG:
            //      writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
            nIters = to!int(iterCount);
            linSolResid = (resid/beta0).re;
            if ( resid <= outerTol ) {
                m = j+1;
                // DEBUG:
                //      writefln("OUTER: TOL ACHIEVED restart-count= %d iteration-count= %d, resid= %e", r, m, resid);
                //      writefln("RANK %d: tolerance achieved on iteration: %d", GlobalConfig.mpi_rank_for_local_task, m);
                break;
            }
        }

        if (iterCount == maxIters)
            m = maxIters;

        // At end H := R up to row m
        //        g := gm up to row m
        upperSolve!number(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot!number(blk.V, blk.nvars, m, blk.g1, blk.zed);
        }

        // apply scaling
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            size_t nsp = blk.myConfig.gmodel.n_species;
            size_t nmodes = blk.myConfig.gmodel.n_modes;
            auto cqi = blk.myConfig.cqi;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                blk.zed[cellCount+MASS] *= (blk.maxRate.vec[cqi.mass]);
                blk.zed[cellCount+X_MOM] *= (blk.maxRate.vec[cqi.xMom]);
                blk.zed[cellCount+Y_MOM] *= (blk.maxRate.vec[cqi.yMom]);
                if ( blk.myConfig.dimensions == 3 )
                    blk.zed[cellCount+Z_MOM] *= (blk.maxRate.vec[cqi.zMom]);
                blk.zed[cellCount+TOT_ENERGY] *= (blk.maxRate.vec[cqi.totEnergy]);
                foreach(it; 0 .. nturb){
                    blk.zed[cellCount+TKE+it] *= (blk.maxRate.vec[cqi.rhoturb+it]);
                }
                version(multi_species_gas){
                if ( nsp > 1 ) {
                    foreach(sp; 0 .. nsp){ blk.zed[cellCount+SPECIES+sp] *= (blk.maxRate.vec[cqi.species+sp]); }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. nmodes){ blk.zed[cellCount+MODES+imode] *= (blk.maxRate.vec[cqi.modes+imode]); }
                }

                cellCount += nConserved;
            }
        }

        // apply preconditioning
        if (usePreconditioner && step >= GlobalConfig.sssOptions.startPreconditioning) {
            final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
                case PreconditionMatrixType.jacobi:
                    foreach(blk; parallel(localFluidBlocks,1)) {
                        nm.smla.multiply(blk.flowJacobian.local, blk.zed, blk.dU);
                    }
                    break;
                case PreconditionMatrixType.ilu:
                    foreach(blk; parallel(localFluidBlocks,1)) {
                        blk.dU[] = blk.zed[];
                        nm.smla.solve(blk.flowJacobian.local, blk.dU);
                    }
                    break;
                case PreconditionMatrixType.sgs:
                    foreach(blk; parallel(localFluidBlocks,1)) {
                        blk.dU[] = blk.zed[];
                        nm.smla.sgs(blk.flowJacobian.local, blk.flowJacobian.diagonal, blk.dU, to!int(nConserved), blk.flowJacobian.Dinv, blk.flowJacobian.Dinv);
                    }
                    break;
                case PreconditionMatrixType.sgs_relax:
                    int local_kmax = GlobalConfig.sssOptions.maxSubIterations;
                    foreach(blk; parallel(localFluidBlocks,1)) {
                        blk.dU[] = blk.zed[];
                        nm.smla.sgsr(blk.flowJacobian.local, blk.dU, blk.flowJacobian.x, to!int(nConserved), local_kmax, blk.flowJacobian.Dinv);
                    }
                    break;
                case PreconditionMatrixType.lu_sgs:
                    mixin(lusgs_solve("dU", "zed"));
                    break;
            } // end switch
        }
        else {
            foreach(blk; parallel(localFluidBlocks,1)) {
                blk.dU[] = blk.zed[];
            }
        }


        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) blk.dU[k] += blk.x0[k];
        }

        if ( resid <= outerTol || r+1 == maxRestarts ) {
            // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
            // DEBUG:  writefln("Breaking restart loop.");
            // DEBUG:  writefln("RANK %d: breaking restart loop, resid= %e, r+1= %d", GlobalConfig.mpi_rank_for_local_task, resid, r+1);
            break;
        }

        // Else, we prepare for restart by setting x0 and computing r0
        // Computation of r0 as per Fraysee etal (2005)
        foreach (blk; parallel(localFluidBlocks,1)) {
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

        /*
        // apply scaling
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            size_t nsp = blk.myConfig.gmodel.n_species;
            size_t nmodes = blk.myConfig.gmodel.n_modes;
            auto cqi = blk.myConfig.cqi;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                blk.r0[cellCount+MASS] *= (1.0/blk.maxRate.vec[cqi.mass]);
                blk.r0[cellCount+X_MOM] *= (1.0/blk.maxRate.vec[cqi.xMom]);
                blk.r0[cellCount+Y_MOM] *= (1.0/blk.maxRate.vec[cqi.yMom]);
                if ( blk.myConfig.dimensions == 3 )
                    blk.r0[cellCount+Z_MOM] *= (1.0/blk.maxRate.vec[cqi.zMom]);
                blk.r0[cellCount+TOT_ENERGY] *= (1.0/blk.maxRate.vec[cqi.totEnergy]);
                foreach(it; 0 .. nturb){
                    blk.r0[cellCount+TKE+it] *= (1.0/blk.maxRate.vec[cqi.rhoturb+it]);
                }
                version(multi_species_gas){
                if ( nsp > 1 ) {
                    foreach(sp; 0 .. nsp){ blk.r0[cellCount+SPECIES+sp] *= (1.0/blk.maxRate.vec[cqi.species+sp]); }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. nmodes){ blk.r0[cellCount+MODES+imode] *= (1.0/blk.maxRate.vec[cqi.modes+imode]); }
                }
                cellCount += nConserved;
            }
        }
        */

        mixin(dot_over_blocks("beta", "r0", "r0"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
        }
        beta = sqrt(beta);

        // DEBUG: writefln("OUTER: ON RESTART beta= %e", beta);
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.r0[k]/beta;
                blk.V[k,0] = blk.v[k];
            }
        }
        // Re-initialise some vectors and matrices for restart
        g0[] = to!number(0.0);
        g1[] = to!number(0.0);
        H0.zeros();
        H1.zeros();
        // And set first residual entry
        g0[0] = beta;

    }

    residual = unscaledNorm2;
    nRestarts = to!int(r);
}

void max_residuals(ConservedQuantities residuals)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;

    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        size_t nsp = blk.myConfig.gmodel.n_species;
        size_t nmodes = blk.myConfig.gmodel.n_modes;
        auto cqi = blk.myConfig.cqi;
        blk.residuals.copy_values_from(blk.cells[0].dUdt[0]);
        blk.residuals.vec[cqi.mass] = fabs(blk.residuals.vec[cqi.mass]);
        blk.residuals.vec[cqi.xMom] = fabs(blk.residuals.vec[cqi.xMom]);
        blk.residuals.vec[cqi.yMom] = fabs(blk.residuals.vec[cqi.yMom]);
        if ( blk.myConfig.dimensions == 3 )
            blk.residuals.vec[cqi.zMom] = fabs(blk.residuals.vec[cqi.zMom]);
        blk.residuals.vec[cqi.totEnergy] = fabs(blk.residuals.vec[cqi.totEnergy]);
        foreach(it; 0 .. nturb){
            blk.residuals.vec[cqi.rhoturb+it] = fabs(blk.residuals.vec[cqi.rhoturb+it]);
        }
        version(multi_species_gas){
        if ( nsp > 1 ) {
            foreach(sp; 0 .. nsp){ blk.residuals.vec[cqi.species+sp] = fabs(blk.residuals.vec[cqi.species+sp]); }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. nmodes){ blk.residuals.vec[cqi.modes+imode] = fabs(blk.residuals.vec[cqi.modes+imode]); }
        }

        number massLocal, xMomLocal, yMomLocal, zMomLocal, energyLocal;
        number[2] turbLocal;
        // we currently expect no more than 32 species
        number[32] speciesLocal;
        number[2] modesLocal;

        foreach (cell; blk.cells) {
            massLocal = cell.dUdt[0].vec[cqi.mass];
            xMomLocal = cell.dUdt[0].vec[cqi.xMom];
            yMomLocal = cell.dUdt[0].vec[cqi.yMom];
            zMomLocal = cell.dUdt[0].vec[cqi.zMom];
            energyLocal = cell.dUdt[0].vec[cqi.totEnergy];
            foreach(it; 0 .. nturb){
                turbLocal[it] = cell.dUdt[0].vec[cqi.rhoturb+it];
            }
            version(multi_species_gas){
            foreach(sp; 0 .. nsp){ speciesLocal[sp] = cell.dUdt[0].vec[cqi.species+sp]; }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes){ modesLocal[imode] = cell.dUdt[0].vec[cqi.modes+imode]; }
            }
            blk.residuals.vec[cqi.mass] = fmax(blk.residuals.vec[cqi.mass], massLocal);
            blk.residuals.vec[cqi.xMom] = fmax(blk.residuals.vec[cqi.xMom], xMomLocal);
            blk.residuals.vec[cqi.yMom] = fmax(blk.residuals.vec[cqi.yMom], yMomLocal);
            if ( blk.myConfig.dimensions == 3 )
                blk.residuals.vec[cqi.zMom] = fmax(blk.residuals.vec[cqi.zMom], zMomLocal);
            blk.residuals.vec[cqi.totEnergy] = fmax(blk.residuals.vec[cqi.totEnergy], energyLocal);
            foreach(it; 0 .. nturb){
                blk.residuals.vec[cqi.rhoturb+it] = fmax(blk.residuals.vec[cqi.rhoturb+it], turbLocal[it]);
            }
            version(multi_species_gas){
            if ( nsp > 1 ) {
                foreach(sp; 0 .. nsp){ blk.residuals.vec[cqi.species+sp] = fmax(blk.residuals.vec[cqi.species+sp], speciesLocal[sp]); }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. nmodes){ blk.residuals.vec[cqi.modes+imode] = fmax(blk.residuals.vec[cqi.modes+imode], modesLocal[imode]); }
            }
        }
    }
    residuals.copy_values_from(localFluidBlocks[0].residuals);
    foreach (blk; localFluidBlocks) {
        auto cqi = blk.myConfig.cqi;
        residuals.vec[cqi.mass] = fmax(residuals.vec[cqi.mass], blk.residuals.vec[cqi.mass]);
        residuals.vec[cqi.xMom] = fmax(residuals.vec[cqi.xMom], blk.residuals.vec[cqi.xMom]);
        residuals.vec[cqi.yMom] = fmax(residuals.vec[cqi.yMom], blk.residuals.vec[cqi.yMom]);
        if ( blk.myConfig.dimensions == 3 )
            residuals.vec[cqi.zMom] = fmax(residuals.vec[cqi.zMom], blk.residuals.vec[cqi.zMom]);
        residuals.vec[cqi.totEnergy] = fmax(residuals.vec[cqi.totEnergy], blk.residuals.vec[cqi.totEnergy]);
        foreach(it; 0 .. blk.myConfig.turb_model.nturb){
            residuals.vec[cqi.rhoturb+it] = fmax(residuals.vec[cqi.rhoturb+it], blk.residuals.vec[cqi.rhoturb+it]);
        }
        version(multi_species_gas){
        if ( blk.myConfig.gmodel.n_species > 1 ) {
            foreach(sp; 0 .. blk.myConfig.gmodel.n_species ){
                residuals.vec[cqi.species+sp] = fmax(residuals.vec[cqi.species+sp], blk.residuals.vec[cqi.species+sp]);
            }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. blk.myConfig.gmodel.n_modes ){
            residuals.vec[cqi.modes+imode] = fmax(residuals.vec[cqi.modes+imode], blk.residuals.vec[cqi.modes+imode]);
        }
        }

    }
    version(mpi_parallel) {
        // In MPI context, only the master task (rank 0) has collated the residuals
        double maxResid;
        auto cqi = GlobalConfig.cqi;
        MPI_Reduce(&(residuals.vec[cqi.mass].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.vec[cqi.mass].re = maxResid; }
        MPI_Reduce(&(residuals.vec[cqi.xMom].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.vec[cqi.xMom] = to!number(maxResid); }
        MPI_Reduce(&(residuals.vec[cqi.yMom].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.vec[cqi.yMom] = to!number(maxResid); }
        if (GlobalConfig.dimensions == 3) {
            MPI_Reduce(&(residuals.vec[cqi.zMom].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals.vec[cqi.zMom] = to!number(maxResid); }
        }
        MPI_Reduce(&(residuals.vec[cqi.totEnergy].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.vec[cqi.totEnergy].re = maxResid; }
        foreach (it; 0 .. GlobalConfig.turb_model.nturb) {
            MPI_Reduce(&(residuals.vec[cqi.rhoturb+it].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals.vec[cqi.rhoturb+it].re = maxResid; }
        }
        version(multi_species_gas){
        if ( GlobalConfig.gmodel_master.n_species > 1 ) {
            foreach (sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                MPI_Reduce(&(residuals.vec[cqi.species+sp].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                if (GlobalConfig.is_master_task) { residuals.vec[cqi.species+sp].re = maxResid; }
            }
        }
        }
        version(multi_T_gas){
        foreach (imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
            MPI_Reduce(&(residuals.vec[cqi.modes+imode].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals.vec[cqi.modes+imode].re = maxResid; }
        }
        }
    }
}

void rewrite_times_file(RestartInfo[] times)
{
    auto fname = format("./config/%s.times", GlobalConfig.base_file_name);
    auto f = File(fname, "w");
    f.writeln("# tindx sim_time dt_global");
    auto cqi = GlobalConfig.cqi;
    foreach (i, rInfo; times) {
        if ( GlobalConfig.dimensions == 2 ) {
            f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e",
                     i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                     rInfo.globalResidual.re, rInfo.residuals.vec[cqi.mass].re,
                     rInfo.residuals.vec[cqi.xMom].re, rInfo.residuals.vec[cqi.yMom].re,
                     rInfo.residuals.vec[cqi.totEnergy].re);
        }
        else {
            f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e %.18e",
                     i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                     rInfo.globalResidual.re, rInfo.residuals.vec[cqi.mass].re,
                     rInfo.residuals.vec[cqi.xMom].re, rInfo.residuals.vec[cqi.yMom].re, rInfo.residuals.vec[cqi.zMom].re,
                     rInfo.residuals.vec[cqi.totEnergy].re);
        }
        foreach(it; 0 .. GlobalConfig.turb_model.nturb){
            f.writef(" %.18e", rInfo.residuals.vec[cqi.rhoturb+it].re);
        }
        version(multi_species_gas){
        if ( GlobalConfig.gmodel_master.n_species > 1 ) {
            foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                f.writef(" %.18e", rInfo.residuals.vec[cqi.species+sp].re);
            }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
            f.writef(" %.18e", rInfo.residuals.vec[cqi.modes+imode].re);
        }
        }
        f.write("\n");
    }
    f.close();
}
