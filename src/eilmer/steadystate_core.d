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

import nm.bbla;
import fvcell;
import fvinterface;
import geom;
import special_block_init;
import fluidblock;
import sfluidblock;
import globaldata;
import globalconfig;
import simcore;
import fvcore;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess : readTimesFile;
import loads;
import shape_sensitivity_core : sss_preconditioner_initialisation, sss_preconditioner;
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

    this(int n_species, int n_modes)
    {
        residuals = new ConservedQuantities(n_species, n_modes);
    }
}

void extractRestartInfoFromTimesFile(string jobName, ref RestartInfo[] times)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;

    auto gmodel = GlobalConfig.gmodel_master;
    RestartInfo restartInfo = RestartInfo(gmodel.n_species, gmodel.n_modes);
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
            restartInfo.residuals.mass = to!double(tokens[6+MASS]);
            restartInfo.residuals.momentum.refx = to!double(tokens[6+X_MOM]);
            restartInfo.residuals.momentum.refy = to!double(tokens[6+Y_MOM]);
            if ( GlobalConfig.dimensions == 3 ) 
                restartInfo.residuals.momentum.refz = to!double(tokens[6+Z_MOM]);
            restartInfo.residuals.total_energy = to!double(tokens[6+TOT_ENERGY]);
            foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
                restartInfo.residuals.rhoturb[it] = to!double(tokens[6+TKE+it]);
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
    bool withPTC = true;
    string jobName = GlobalConfig.base_file_name;
    int nsteps = GlobalConfig.sssOptions.nTotalSteps;
    int nRestarts;
    int maxNumberAttempts = GlobalConfig.sssOptions.maxNumberAttempts;
    double relGlobalResidReduction = GlobalConfig.sssOptions.stopOnRelGlobalResid;
    double absGlobalResidReduction = GlobalConfig.sssOptions.stopOnAbsGlobalResid;
    double cfl_max = GlobalConfig.sssOptions.cfl_max;
    int cfl_schedule_current_index = 0;
    auto cfl_schedule_value_list = GlobalConfig.sssOptions.cfl_schedule_value_list;
    auto cfl_schedule_iter_list = GlobalConfig.sssOptions.cfl_schedule_iter_list;
    bool residual_based_cfl_scheduling = GlobalConfig.sssOptions.residual_based_cfl_scheduling;
    int kmax = GlobalConfig.sssOptions.maxSubIterations;
    // Settings for start-up phase
    double cfl0 = GlobalConfig.sssOptions.cfl0;
    double tau0 = GlobalConfig.sssOptions.tau0;
    double eta0 = GlobalConfig.sssOptions.eta0;
    double sigma0 = GlobalConfig.sssOptions.sigma0;
    // Settings for inexact Newton phase
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
    
    int n_species = GlobalConfig.gmodel_master.n_species();
    int n_modes = GlobalConfig.gmodel_master.n_modes();
    ConservedQuantities maxResiduals = new ConservedQuantities(n_species, n_modes);
    ConservedQuantities currResiduals = new ConservedQuantities(n_species, n_modes);
    number mass_balance = 0.0;
    
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;

    double cfl, cflTrial;
    double dt;
    double etaTrial;
    int failedAttempt = 0;
    double pseudoSimTime = 0.0;
    double normOld, normNew;
    int snapshotsCount = GlobalConfig.sssOptions.snapshotsCount;
    int nTotalSnapshots = GlobalConfig.sssOptions.nTotalSnapshots;
    int nWrittenSnapshots = 0;
    int writeDiagnosticsCount = GlobalConfig.sssOptions.writeDiagnosticsCount;
    int writeLoadsCount = GlobalConfig.sssOptions.writeLoadsCount;

    int startStep;
    int nPreSteps = GlobalConfig.sssOptions.nPreSteps;
    int nStartUpSteps = GlobalConfig.sssOptions.nStartUpSteps;
    bool inexactNewtonPhase = false;

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
        foreach (blk; localFluidBlocks) {
            // Make a block-local copy of conserved quantities info
            blk.nConserved = nConservedQuantities;
            blk.MASS = massIdx;
            blk.X_MOM = xMomIdx;
            blk.Y_MOM = yMomIdx;
            blk.Z_MOM = zMomIdx;
            blk.TOT_ENERGY = totEnergyIdx;
            blk.TKE = tkeIdx;
            sss_preconditioner_initialisation(blk, nConserved); 
        }
    }
    // Set usePreconditioner to false for pre-steps AND first-order steps.
    usePreconditioner = false;

    // some data structures for use with LU-SGS
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach(c; blk.cells) {
            c.dFdU = new Matrix!number(5,5); // number of conserved variables
            c.dFdU_tmp = new Matrix!number(5,5);
        }
        foreach(f; blk.faces) {
            f.T = new Matrix!number(5,5); // number of conserved variables
            f.Tinv = new Matrix!number(5,5);

            f.T[0,0] = to!number(1.0); f.T[0,1] = to!number(0.0); f.T[0,2] = to!number(0.0); f.T[0,3] = to!number(0.0); f.T[0,4] = to!number(0.0);
            f.T[1,0] = to!number(0.0); f.T[1,1] = f.n.x; f.T[1,2] = f.n.y; f.T[1,3] = f.n.z; f.T[1,4] = to!number(0.0);
            f.T[2,0] = to!number(0.0); f.T[2,1] = f.t1.x; f.T[2,2] = f.t1.y; f.T[2,3] = f.t1.z; f.T[2,4] = to!number(0.0);
            f.T[3,0] = to!number(0.0); f.T[3,1] = f.t2.x; f.T[3,2] = f.t2.y; f.T[3,3] = f.t2.z; f.T[3,4] = to!number(0.0);
            f.T[4,0] = to!number(0.0); f.T[4,1] = to!number(0.0); f.T[4,2] = to!number(0.0); f.T[4,3] = to!number(0.0); f.T[4,4] = to!number(1.0);
            
            f.Tinv = inverse(f.T);
        }
    }
    
    // We only do a pre-step phase if we are starting from scratch.
    if ( snapshotStart == 0 ) {
        cfl = cfl0;
        dt = determine_dt(cfl);
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        }
        // The initial residual is usually a poor choice for basing decisions about how the
        // residual is dropping, particularly when a constant initial condition is given.
        // A constant initial condition gives a zero residual everywhere in the interior
        // of the domain. Therefore, the only residual can come from boundary condition influences.
        // What we do here is use some early steps (I've called them "pre"-steps) to
        // allow the boundary conditions to exert some influence into the interior.
        // We'll find the max residual in that start-up process and use that as our initial value
        // for the "max" residual. Our relative residuals will then be based on that.
        // We can use a fixed timestep (typically small) and low-order reconstruction in this
        // pre-step phase.

        if (GlobalConfig.is_master_task) {
            writeln("---------------------------------------------------------------");
            writeln("  Begin pre-steps to establish sensible max residual.");
            writeln("---------------------------------------------------------------");
        }
        foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(1);
        foreach ( preStep; -nPreSteps .. 0 ) {
            foreach (attempt; 0 .. maxNumberAttempts) {
                try {
                    rpcGMRES_solve(preStep, pseudoSimTime, dt, eta0, sigma0, usePreconditioner, normOld, nRestarts, startStep);
                    //lusgs_solve(preStep, pseudoSimTime, dt, 2.0, normOld, kmax, startStep);
                }
                catch (FlowSolverException e) {
                    version(mpi_parallel) {
                        writefln("Failed when attempting GMRES solve in pre-steps on task %d.", GlobalConfig.mpi_rank_for_local_task);
                    }
                    else {
                        writeln("Failed when attempting GMRES solve in pre-steps.");
                    }
                    writefln("attempt %d: dt= %e", attempt, dt);
                    failedAttempt = 1;
                    cfl = 0.1*cfl;
                    dt = determine_dt(cfl);
                }
                // Coordinate MPI tasks after try-catch statement in case one or more of the tasks encountered an exception.
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &failedAttempt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                }
                if (failedAttempt > 0) { break; }
                // end: coordination of MPI tasks
                foreach (blk; parallel(localFluidBlocks,1)) {
                    size_t nturb = blk.myConfig.turb_model.nturb;
                    int cellCount = 0;
                    foreach (cell; blk.cells) {
                        cell.U[1].copy_values_from(cell.U[0]);
                        cell.U[1].mass = cell.U[0].mass + blk.dU[cellCount+MASS];
                        cell.U[1].momentum.refx = cell.U[0].momentum.x + blk.dU[cellCount+X_MOM];
                        cell.U[1].momentum.refy = cell.U[0].momentum.y + blk.dU[cellCount+Y_MOM];
                        if ( blk.myConfig.dimensions == 3 ) 
                            cell.U[1].momentum.refz = cell.U[0].momentum.z + blk.dU[cellCount+Z_MOM];
                        cell.U[1].total_energy = cell.U[0].total_energy + blk.dU[cellCount+TOT_ENERGY];
                        foreach(it; 0 .. nturb) {
                            cell.U[1].rhoturb[it] = cell.U[0].rhoturb[it] + blk.dU[cellCount+TKE+it];
                        }
                        // enforce mass fraction of 1 for single species gas
                        if (blk.myConfig.n_species == 1) {
                            cell.U[1].massf[0] = cell.U[1].mass;
                        } 
                        try {
                            cell.decode_conserved(0, 1, 0.0);
                        }
                        catch (FlowSolverException e) {
                            version(mpi_parallel) {
                                writefln("Failed to provide sensible update on task %d.", GlobalConfig.mpi_rank_for_local_task);
                            }
                            else {
                                writeln("Failed to provide sensible update.");
                            }
                            writefln("attempt %d: dt= %e", attempt, dt);
                            failedAttempt = 1;
                        }
                        cellCount += nConserved;
                    }
                }
                // Coordinate MPI tasks if one of them had a failed attempt.
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &failedAttempt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }

                if (failedAttempt > 0) {
                    // reduce CFL
                    cfl = 0.1*cfl;
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
                    continue;
		}

                // If we get here, things are good. Put flow state into U[0]
                // ready for next iteration.
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (cell; blk.cells) {
                        swap(cell.U[0], cell.U[1]);
                    }
                }
                // If we got here, we can break the attempts loop.
                break;
            }

            if (failedAttempt > 0) {
                if (GlobalConfig.is_master_task) {
                    writefln("Pre-step failed: %d", SimState.step);
                    writeln("Bailing out!");
                }
                exit(1);
            }
            if (GlobalConfig.is_master_task) {
                writefln("PRE-STEP  %d  ::  dt= %.3e  global-norm= %.12e", preStep, dt, normOld);
            }
            if ( normOld > normRef ) {
                normRef = normOld;
                max_residuals(maxResiduals);
            }
        }

        if (GlobalConfig.is_master_task) {
            writeln("Pre-step phase complete.");
        }

        if (nPreSteps <= 0) {
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(interpOrderSave);
            // Take initial residual as max residual
            evalRHS(0.0, 0);
            max_residuals(maxResiduals);
            foreach (blk; parallel(localFluidBlocks, 1)) {
                size_t nturb = blk.myConfig.turb_model.nturb;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.FU[cellCount+MASS] = -cell.dUdt[0].mass;
                    blk.FU[cellCount+X_MOM] = -cell.dUdt[0].momentum.x;
                    blk.FU[cellCount+Y_MOM] = -cell.dUdt[0].momentum.y;
                    if ( GlobalConfig.dimensions == 3 )
                        blk.FU[cellCount+Z_MOM] = -cell.dUdt[0].momentum.z;
                    blk.FU[cellCount+TOT_ENERGY] = -cell.dUdt[0].total_energy;
                    foreach(it; 0 .. nturb) {
                        blk.FU[cellCount+TKE+it] = -cell.dUdt[0].rhoturb[it];
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
        if (GlobalConfig.is_master_task) {
            writeln("Reference residuals are established as:");
            writefln("GLOBAL:         %.12e", normRef);
            writefln("MASS:           %.12e", maxResiduals.mass.re);
            writefln("X-MOMENTUM:     %.12e", maxResiduals.momentum.x.re);
            writefln("Y-MOMENTUM:     %.12e", maxResiduals.momentum.y.re);
            if ( GlobalConfig.dimensions == 3 )
                writefln("Z-MOMENTUM:     %.12e", maxResiduals.momentum.z.re);
            writefln("ENERGY:         %.12e", maxResiduals.total_energy.re);
            foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
                string tvname = capitalize(GlobalConfig.turb_model.primitive_variable_name(it));
                writefln("%s:            %.12e",tvname, maxResiduals.rhoturb[it].re);
            }
        
            string refResidFname = jobName ~ "-ref-residuals.saved";
            auto refResid = File(refResidFname, "w");
            if ( GlobalConfig.dimensions == 2 ) {
                refResid.writef("%.18e %.18e %.18e %.18e %.18e",
                                normRef, maxResiduals.mass.re, maxResiduals.momentum.x.re,
                                maxResiduals.momentum.y.re, maxResiduals.total_energy.re);
            }
            else {
                refResid.writef("%.18e %.18e %.18e %.18e %.18e %.18e",
                                normRef, maxResiduals.mass.re, maxResiduals.momentum.x.re,
                                maxResiduals.momentum.y.re, maxResiduals.momentum.z.re,
                                maxResiduals.total_energy.re);
            }
            foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
                refResid.writef(" %.18e", maxResiduals.rhoturb[it].re);
            }
            refResid.write("\n");
            refResid.close();
        }
    }

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
        maxResiduals.mass = to!double(tokens[1+MASS]);
        maxResiduals.momentum.refx = to!double(tokens[1+X_MOM]);
        maxResiduals.momentum.refy = to!double(tokens[1+Y_MOM]);
        if ( GlobalConfig.dimensions == 3 ) 
            maxResiduals.momentum.refz = to!double(tokens[1+Z_MOM]);
        maxResiduals.total_energy = to!double(tokens[1+TOT_ENERGY]);
        foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
            maxResiduals.rhoturb[it] = to!double(tokens[1+TKE+it]);
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
            fResid.writeln("#  7: nFnCalls");
            fResid.writeln("#  8: wall-clock, s");
            fResid.writeln("#  9: global-residual-abs");
            fResid.writeln("# 10: global-residual-rel");
            fResid.writefln("#  %02d: mass-abs", 11+2*MASS);
            fResid.writefln("# %02d: mass-rel", 11+2*MASS+1);
            fResid.writefln("# %02d: x-mom-abs", 11+2*X_MOM);
            fResid.writefln("# %02d: x-mom-rel", 11+2*X_MOM+1);
            fResid.writefln("# %02d: y-mom-abs", 11+2*Y_MOM);
            fResid.writefln("# %02d: y-mom-rel", 11+2*Y_MOM+1);
            if ( GlobalConfig.dimensions == 3 ) {
                fResid.writefln("# %02d: z-mom-abs", 11+2*Z_MOM);
                fResid.writefln("# %02d: z-mom-rel", 11+2*Z_MOM+1);
            }
            fResid.writefln("# %02d: energy-abs", 11+2*TOT_ENERGY);
            fResid.writefln("# %02d: energy-rel", 11+2*TOT_ENERGY+1);
            auto nt = GlobalConfig.turb_model.nturb;
            foreach(it; 0 .. nt) {
                string tvname = GlobalConfig.turb_model.primitive_variable_name(it);
                fResid.writefln("# %02d: %s-abs", 11+2*(TKE+it), tvname);
                fResid.writefln("# %02d: %s-rel", 11+2*(TKE+it)+1, tvname);
            }
            fResid.writeln("# %02d: mass-balance", 11+2*(TKE+nt)+2);
            fResid.close();
        }
    }

    // Begin Newton steps
    double eta;
    double tau;
    double sigma;
    foreach (step; startStep .. nsteps+1) {
        if ( (step/GlobalConfig.control_count)*GlobalConfig.control_count == step ) {
            read_control_file(); // Reparse the time-step control parameters occasionally.
        }
        residualsUpToDate = false;
        if ( step <= nStartUpSteps ) {
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(1);
            eta = eta0;
            tau = tau0;
            sigma = sigma0;
            usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
        }
        else {
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(interpOrderSave);
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
        foreach (attempt; 0 .. maxNumberAttempts) {
            failedAttempt = 0;
            try {
                rpcGMRES_solve(step, pseudoSimTime, dt, eta, sigma, usePreconditioner, normNew, nRestarts, startStep);
                //lusgs_solve(step, pseudoSimTime, dt, 2.0, normNew, kmax, startStep);
            }
            catch (FlowSolverException e) {
                writefln("Failed when attempting GMRES solve in main steps.");
                writefln("attempt %d: dt= %e", attempt, dt);
                failedAttempt = 1;
                cfl = 0.1*cfl;
                dt = determine_dt(cfl);
            }
            // Coordinate MPI tasks after try-catch statement in case one or more of the tasks encountered an exception.
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &failedAttempt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
            if (failedAttempt > 0) { continue; }
            
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.turb_model.nturb;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    cell.U[1].copy_values_from(cell.U[0]);
                    cell.U[1].mass = cell.U[0].mass + blk.dU[cellCount+MASS];
                    cell.U[1].momentum.refx = cell.U[0].momentum.x + blk.dU[cellCount+X_MOM];
                    cell.U[1].momentum.refy = cell.U[0].momentum.y + blk.dU[cellCount+Y_MOM];
                    if ( blk.myConfig.dimensions == 3 ) 
                        cell.U[1].momentum.refz = cell.U[0].momentum.z + blk.dU[cellCount+Z_MOM];
                    cell.U[1].total_energy = cell.U[0].total_energy + blk.dU[cellCount+TOT_ENERGY];
                    foreach(it; 0 .. nturb){
                        cell.U[1].rhoturb[it] = cell.U[0].rhoturb[it] + blk.dU[cellCount+TKE+it];
                    }
                    // enforce mass fraction of 1 for single species gas
                    if (blk.myConfig.n_species == 1) {
                        cell.U[1].massf[0] = cell.U[1].mass;
                    } 
                    try {
                        cell.decode_conserved(0, 1, 0.0);
                    }
                    catch (FlowSolverException e) {
                        writefln("Failed attempt %d: dt= %e", attempt, dt);
                        failedAttempt = 1;
                    }
                    cellCount += nConserved;
                }
            }
            
            // Coordinate MPI tasks after try-catch statement in case one or more of the tasks encountered an exception.
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &failedAttempt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            // end: coordination of MPI tasks
            if (failedAttempt > 0) {
                // Reduce CFL
                cfl = 0.1*cfl;
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
                continue;
	    }

            // If we get here, things are good. Put flow state into U[0]
            // ready for next iteration.
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (cell; blk.cells) {
                    swap(cell.U[0], cell.U[1]);
                }
            }
            // If we got here, we can break the attempts loop.
            break;
        }
        if (failedAttempt > 0) {
            if (GlobalConfig.is_master_task) {
                writefln("Step failed: %d", step);
                writeln("Bailing out!");
            }
            exit(1);
        }

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
                fResid = File(residFname, "a");
                fResid.writef("%8d  %20.16e  %20.16e %20.16e %20.16e %3d %5d %.8f %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e ",
                              step, pseudoSimTime, dt, cfl, eta, nRestarts, fnCount, wallClockElapsed, 
                              normNew, normNew/normRef,
                              currResiduals.mass.re, currResiduals.mass.re/maxResiduals.mass.re,
                              currResiduals.momentum.x.re, currResiduals.momentum.x.re/maxResiduals.momentum.x.re,
                              currResiduals.momentum.y.re, currResiduals.momentum.y.re/maxResiduals.momentum.y.re);
                if ( GlobalConfig.dimensions == 3 )
                    fResid.writef("%20.16e  %20.16e  ", currResiduals.momentum.z.re, currResiduals.momentum.z.re/maxResiduals.momentum.z.re);
                fResid.writef("%20.16e  %20.16e  ",
                              currResiduals.total_energy.re, currResiduals.total_energy.re/maxResiduals.total_energy.re);
                foreach(it; 0 .. GlobalConfig.turb_model.nturb){
                    fResid.writef("%20.16e  %20.16e  ",
                                  currResiduals.rhoturb[it].re, currResiduals.rhoturb[it].re/maxResiduals.rhoturb[it].re);
                }
                fResid.writef("%20.16e ", fabs(mass_balance.re));
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
                auto writer = appender!string();
                formattedWrite(writer, "STEP= %7d  pseudo-time=%10.3e dt=%10.3e cfl=%10.3e  WC=%.1f \n", step, pseudoSimTime, dt, cfl, wallClockElapsed);
                formattedWrite(writer, "RESIDUALS        absolute        relative\n");
                formattedWrite(writer, "  global         %10.6e    %10.6e\n", normNew, normNew/normRef);
                formattedWrite(writer, "  mass           %10.6e    %10.6e\n", currResiduals.mass.re, currResiduals.mass.re/maxResiduals.mass.re);
                formattedWrite(writer, "  x-mom          %10.6e    %10.6e\n", currResiduals.momentum.x.re, currResiduals.momentum.x.re/maxResiduals.momentum.x.re);
                formattedWrite(writer, "  y-mom          %10.6e    %10.6e\n", currResiduals.momentum.y.re, currResiduals.momentum.y.re/maxResiduals.momentum.y.re);
                if ( GlobalConfig.dimensions == 3 )
                    formattedWrite(writer, "  z-mom          %10.6e    %10.6e\n", currResiduals.momentum.z.re, currResiduals.momentum.z.re/maxResiduals.momentum.z.re);
                formattedWrite(writer, "  total-energy   %10.6e    %10.6e\n", currResiduals.total_energy.re, currResiduals.total_energy.re/maxResiduals.total_energy.re);
                foreach(it; 0 .. GlobalConfig.turb_model.nturb){
                    auto tvname = GlobalConfig.turb_model.primitive_variable_name(it);
                    formattedWrite(writer, "  %s            %10.6e    %10.6e\n", tvname, currResiduals.rhoturb[it].re, currResiduals.rhoturb[it].re/maxResiduals.rhoturb[it].re);
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
            nWrittenSnapshots++;
            if ( nWrittenSnapshots <= nTotalSnapshots ) {
                if (GlobalConfig.is_master_task) { ensure_directory_is_present(make_path_name!"flow"(nWrittenSnapshots)); }
                version(mpi_parallel) {
                    MPI_Barrier(MPI_COMM_WORLD);
                }
                foreach (blk; localFluidBlocks) {
                    auto fileName = make_file_name!"flow"(jobName, blk.id, nWrittenSnapshots, "gz");
                    blk.write_solution(fileName, pseudoSimTime);
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
                // We need to shuffle all of the snapshots...
                foreach ( iSnap; 2 .. nTotalSnapshots+1) {
                    foreach (blk; localFluidBlocks) {
                        auto fromName = make_file_name!"flow"(jobName, blk.id, iSnap, "gz");
                        auto toName = make_file_name!"flow"(jobName, blk.id, iSnap-1, "gz");
                        rename(fromName, toName);
                    }
                }
                // ... and add the new snapshot.
                foreach (blk; localFluidBlocks) {
                    auto fileName = make_file_name!"flow"(jobName, blk.id, nTotalSnapshots, "gz");        
                    blk.write_solution(fileName, pseudoSimTime);
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
	    }
            outFile.close();
	}
    } // end if (GlobalConfig.frozen_limiter)    
}

void allocate_global_workspace()
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
    if (ftl == 0 && GlobalConfig.do_shock_detect) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.detect_shock_points();
        }
    }

    int gtl = 0;
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
        exchange_ghost_cell_boundary_data(pseudoSimTime, gtl, ftl);
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.viscous_flux();
        }
        foreach (blk; localFluidBlocks) {
            blk.applyPostDiffFluxAction(pseudoSimTime, 0, ftl);
        }
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (i, cell; blk.cells) {
            cell.add_inviscid_source_vector(0, 0.0);
            if (blk.myConfig.viscous) {
                cell.add_viscous_source_vector();
            }
            if (blk.myConfig.udf_source_terms) {
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                if (blk.grid_type == Grid_t.structured_grid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                    i_cell = ijk_indices[0];
                    j_cell = ijk_indices[1];
                    k_cell = ijk_indices[2];
                }
                addUDFSourceTermsToCell(blk.myL, cell, 0, 
                                        pseudoSimTime, blk.myConfig,
                                        blk.id, i_cell, j_cell, k_cell);
            }
            cell.time_derivatives(0, ftl);
        }
    }
}

void evalJacobianVecProd(double pseudoSimTime, double sigma)
{
    if (GlobalConfig.sssOptions.useComplexMatVecEval)
        evalComplexMatVecProd(pseudoSimTime, sigma);
    else
        evalRealMatVecProd(pseudoSimTime, sigma);
}

void evalRealMatVecProd(double pseudoSimTime, double sigma)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;

    // We perform a Frechet derivative to evaluate J*D^(-1)v
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
        int cellCount = 0;
        foreach (cell; blk.cells) {
            cell.U[1].copy_values_from(cell.U[0]);
            cell.U[1].mass += sigma*blk.zed[cellCount+MASS];
            cell.U[1].momentum.refx += sigma*blk.zed[cellCount+X_MOM];
            cell.U[1].momentum.refy += sigma*blk.zed[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                cell.U[1].momentum.refz += sigma*blk.zed[cellCount+Z_MOM];
            cell.U[1].total_energy += sigma*blk.zed[cellCount+TOT_ENERGY];
            foreach(it; 0 .. nturb){
                cell.U[1].rhoturb[it] += sigma*blk.zed[cellCount+TKE+it];
            }
            cell.decode_conserved(0, 1, 0.0);
            cellCount += nConserved;
        }
    }
    evalRHS(pseudoSimTime, 1);
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.zed[cellCount+MASS] = (-cell.dUdt[1].mass - blk.FU[cellCount+MASS])/(sigma);
            blk.zed[cellCount+X_MOM] = (-cell.dUdt[1].momentum.x - blk.FU[cellCount+X_MOM])/(sigma);
            blk.zed[cellCount+Y_MOM] = (-cell.dUdt[1].momentum.y - blk.FU[cellCount+Y_MOM])/(sigma);
            if ( blk.myConfig.dimensions == 3 )
                blk.zed[cellCount+Z_MOM] = (-cell.dUdt[1].momentum.z - blk.FU[cellCount+Z_MOM])/(sigma);
            blk.zed[cellCount+TOT_ENERGY] = (-cell.dUdt[1].total_energy - blk.FU[cellCount+TOT_ENERGY])/(sigma);
            foreach(it; 0 .. nturb){
                blk.zed[cellCount+TKE+it] = (-cell.dUdt[1].rhoturb[it] - blk.FU[cellCount+TKE+it])/(sigma);
            }
            cellCount += nConserved;
        }
    }
}

void evalComplexMatVecProd(double pseudoSimTime, double sigma)
{
    version(complex_numbers) {
        // Make a stack-local copy of conserved quantities info
        size_t nConserved = nConservedQuantities;
        size_t MASS = massIdx;
        size_t X_MOM = xMomIdx;
        size_t Y_MOM = yMomIdx;
        size_t Z_MOM = zMomIdx;
        size_t TOT_ENERGY = totEnergyIdx;
        size_t TKE = tkeIdx;
        
        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.clear_fluxes_of_conserved_quantities();
            foreach (cell; blk.cells) cell.clear_source_vector();

            size_t nturb = blk.myConfig.turb_model.nturb;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                cell.U[1].copy_values_from(cell.U[0]);
                cell.U[1].mass += complex(0.0, sigma*blk.zed[cellCount+MASS].re);
                cell.U[1].momentum.refx += complex(0.0, sigma*blk.zed[cellCount+X_MOM].re);
                cell.U[1].momentum.refy += complex(0.0, sigma*blk.zed[cellCount+Y_MOM].re);
                if ( blk.myConfig.dimensions == 3 )
                    cell.U[1].momentum.refz += complex(0.0, sigma*blk.zed[cellCount+Z_MOM].re);
                cell.U[1].total_energy += complex(0.0, sigma*blk.zed[cellCount+TOT_ENERGY].re);
                foreach(it; 0 .. nturb){
                    cell.U[1].rhoturb[it] += complex(0.0, sigma*blk.zed[cellCount+TKE+it].re);
                }
                cell.decode_conserved(0, 1, 0.0);
                cellCount += nConserved;
            }
        }
        evalRHS(pseudoSimTime, 1);
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                blk.zed[cellCount+MASS] = -cell.dUdt[1].mass.im/(sigma);
                blk.zed[cellCount+X_MOM] = -cell.dUdt[1].momentum.x.im/(sigma);
                blk.zed[cellCount+Y_MOM] = -cell.dUdt[1].momentum.y.im/(sigma);
                if ( blk.myConfig.dimensions == 3 )
                    blk.zed[cellCount+Z_MOM] = -cell.dUdt[1].momentum.z.im/(sigma);
                blk.zed[cellCount+TOT_ENERGY] = -cell.dUdt[1].total_energy.im/(sigma);
                foreach(it; 0 .. nturb){
                    blk.zed[cellCount+TKE+it] = -cell.dUdt[1].rhoturb[it].im/(sigma);
                }
                // we must explicitly remove the complex-perturbation for the turbulent properties - failing to do so
                // will cause imaginary values to leak into the flowstate variables (outside of the Frechet derivative).
                cell.fs.mu_t = cell.fs.mu_t.re;
                cell.fs.k_t = cell.fs.k_t.re;
                cellCount += nConserved;
            }
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


void evalMatrixFreeFluxIncrement(FVCell cell, FVInterface face)
{
    // Matrix-Free Flux vector increment
    //
    // As defined on right column, pg 4 of
    // Rieger and Jameson (1988),
    // Solution of Steady Three-Dimensional Compressible Euler and Navier-Stokes Equations by an Implicit LU Scheme
    // AIAA conference paper
    //
    // Uses Roe's split flux scheme for LHS Jacobian as per
    // Luo, Baum, and Lohner (1998)
    // A Fast, Matrix-free Implicit Method for Compressible Flows on Unstructured Grids,
    // Journal of computational physics
    // 

    // make sure cells (particularly ghost cells) have conserved quantities filled
    cell.encode_conserved(0, 0, 0.0);

    // peturb conserved quantities by current approximation of dU
    cell.U[1].copy_values_from(cell.U[0]);
    cell.U[1].mass += cell.dU[0].mass;
    cell.U[1].momentum.refx += cell.dU[0].momentum.x;
    cell.U[1].momentum.refy += cell.dU[0].momentum.y;
    if (GlobalConfig.dimensions == 3 )
        cell.U[1].momentum.refz += cell.dU[0].momentum.z;
    cell.U[1].total_energy += cell.dU[0].total_energy;

    // update primitive variables
    cell.decode_conserved(0, 1, 0.0);          

    // Peturbed state flux
    auto gmodel = GlobalConfig.gmodel_master;
    number rho = cell.fs.gas.rho;
    number velx = cell.fs.vel.dot(face.n);
    number vely = cell.fs.vel.dot(face.t1);
    number velz = cell.fs.vel.dot(face.t2);
    number p = cell.fs.gas.p;
    number e = gmodel.internal_energy(cell.fs.gas);

    cell.dF.mass = rho*velx;
    cell.dF.momentum.refx = p + rho*velx*velx;
    cell.dF.momentum.refy = rho*velx*vely;
    cell.dF.momentum.refz = rho*velx*velz;
    cell.dF.total_energy = (rho*e + rho*(velx^^2 + vely^^2 + velz^^2)/2.0 + p)*velx;

    // reset primitive variables to unperturbed state
    cell.decode_conserved(0, 0, 0.0);          

    // original state flux
    rho = cell.fs.gas.rho;
    velx = cell.fs.vel.dot(face.n);
    vely = cell.fs.vel.dot(face.t1);
    velz = cell.fs.vel.dot(face.t2);
    p = cell.fs.gas.p;
    e = gmodel.internal_energy(cell.fs.gas);

    // flux vector increment
    cell.dF.mass -= rho*velx;
    cell.dF.momentum.refx -= p + rho*velx*velx;
    cell.dF.momentum.refy -= rho*velx*vely;
    cell.dF.momentum.refz -= rho*velx*velz;
    cell.dF.momentum.transform_to_global_frame(face.n, face.t1, face.t2);
    cell.dF.total_energy -= (rho*e + rho*(velx^^2 + vely^^2 + velz^^2)/2.0 + p)*velx;
}


void evalMatrixBasedFluxIncrement(FVCell cell, FVInterface face)
{
    // Matrix-Based Flux vector increment
    //
    // Hand differentiation of Roe's split flux scheme for LHS Jacobian as per
    // Luo, Baum, and Lohner (1998)
    // A Fast, Matrix-free Implicit Method for Compressible Flows on Unstructured Grids,
    // Journal of computational physics
    // 

    // primitive variables
    auto gmodel = GlobalConfig.gmodel_master;
    number gam = gmodel.gamma(cell.fs.gas);
    number rho = cell.fs.gas.rho;
    number u = cell.fs.vel.dot(face.n);
    number v = cell.fs.vel.dot(face.t1);
    number w = cell.fs.vel.dot(face.t2);
    number p = cell.fs.gas.p;
    number e = gmodel.internal_energy(cell.fs.gas);

    // conserved variables
    number U1 = rho; 
    number U2 = rho*u;
    number U3 = rho*v;
    number U4 = rho*w;
    number U5 = rho*e + rho*(u^^2 + v^^2 + w^^2)/2.0;

    // Roe's split flux Jacobian
    cell.dFdU[0,0] = to!number(0.0);
    cell.dFdU[0,1] = to!number(1.0);
    cell.dFdU[0,2] = to!number(0.0);
    cell.dFdU[0,3] = to!number(0.0);
    cell.dFdU[0,4] = to!number(0.0); 

    cell.dFdU[1,0] = -(U2*U2)/(U1*U1) + (gam-1.0)*(U2*U2+U3*U3+U4*U4)/(2.0*U1*U1);
    cell.dFdU[1,1] = (2.0*U2)/U1 + (1.0-gam)*(U2/U1);
    cell.dFdU[1,2] = (1.0-gam)*(U3/U1);
    cell.dFdU[1,3] = (1.0-gam)*(U4/U1);
    cell.dFdU[1,4] = (gam-1.0);
    
    cell.dFdU[2,0] = to!number(0.0);
    cell.dFdU[2,1] = U3/U1;
    cell.dFdU[2,2] = U2/U1;
    cell.dFdU[2,3] = to!number(0.0);
    cell.dFdU[2,4] = to!number(0.0); 

    cell.dFdU[3,0] = to!number(0.0);
    cell.dFdU[3,1] = U4/U1;
    cell.dFdU[3,2] = to!number(0.0);
    cell.dFdU[3,3] = U2/U1;
    cell.dFdU[3,4] = to!number(0.0);
    
    cell.dFdU[4,0] = -gam*(U5*U2)/(U1*U1) + (gam-1.0)*(U2*U2*U2+U2*U3*U3+U2*U4*U4)/(U1*U1);
    cell.dFdU[4,1] = gam*(U5/U1) + (1.0-gam)*(3*U2*U2+U3*U3+U4*U4)/(2*U1*U1);
    cell.dFdU[4,2] = (1.0-gam)*(U3*U2)/(U1*U1);
    cell.dFdU[4,3] = (1.0-gam)*(U4*U2)/(U1*U1);
    cell.dFdU[4,4] = gam*(U2/U1); 

    // rotate matrix back to global frame
    dot(face.Tinv, cell.dFdU, cell.dFdU_tmp);
    dot(cell.dFdU_tmp, face.T, cell.dFdU);

    // flux vector increment (i.e. dFdU * dU)
    cell.dF.mass = cell.dFdU[0,0]*cell.dU[0].mass
        + cell.dFdU[0,1]*cell.dU[0].momentum.x
        + cell.dFdU[0,2]*cell.dU[0].momentum.y
        + cell.dFdU[0,3]*cell.dU[0].momentum.z
        + cell.dFdU[0,4]*cell.dU[0].total_energy;
    
    cell.dF.momentum.refx = cell.dFdU[1,0]*cell.dU[0].mass
        + cell.dFdU[1,1]*cell.dU[0].momentum.x
        + cell.dFdU[1,2]*cell.dU[0].momentum.y
        + cell.dFdU[1,3]*cell.dU[0].momentum.z
        + cell.dFdU[1,4]*cell.dU[0].total_energy;
    
    cell.dF.momentum.refy = cell.dFdU[2,0]*cell.dU[0].mass
        + cell.dFdU[2,1]*cell.dU[0].momentum.x
        + cell.dFdU[2,2]*cell.dU[0].momentum.y
        + cell.dFdU[2,3]*cell.dU[0].momentum.z
        + cell.dFdU[2,4]*cell.dU[0].total_energy;

    cell.dF.momentum.refz = cell.dFdU[3,0]*cell.dU[0].mass
        + cell.dFdU[3,1]*cell.dU[0].momentum.x
        + cell.dFdU[3,2]*cell.dU[0].momentum.y
        + cell.dFdU[3,3]*cell.dU[0].momentum.z
        + cell.dFdU[3,4]*cell.dU[0].total_energy;
    
    cell.dF.total_energy = cell.dFdU[4,0]*cell.dU[0].mass
        + cell.dFdU[4,1]*cell.dU[0].momentum.x
        + cell.dFdU[4,2]*cell.dU[0].momentum.y
        + cell.dFdU[4,3]*cell.dU[0].momentum.z
        + cell.dFdU[4,4]*cell.dU[0].total_energy;
}

number spectral_radius(FVInterface face, bool viscous, double omega) {
    auto gmodel = GlobalConfig.gmodel_master;
    number lam = 0.0;
    auto fvel = fabs(face.fs.vel.dot(face.n));
    lam += omega*(fvel + face.fs.gas.a);
    if (viscous) {
        auto dpos = face.left_cell.pos[0] - face.right_cell.pos[0];
        number dr = sqrt(dpos.dot(dpos));
        number Prandtl = face.fs.gas.mu * gmodel.Cp(face.fs.gas) / face.fs.gas.k;
        lam += (1.0/dr)*fmax(4.0/(3.0*face.fs.gas.rho), gmodel.gamma(face.fs.gas)/face.fs.gas.rho) * (face.fs.gas.mu/Prandtl);
    }
    return lam;
}

void lusgs_solve(int step, double pseudoSimTime, double dt, double omega, ref double residual, int kmax, int startStep)
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
    
    // 1. Evaluate RHS residual (R)
    evalRHS(pseudoSimTime, 0);
    
    // 2. Compute scalar diagonal (D) and dQ^0 = D**(-1)*R
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (cell; blk.cells) {
            number dtInv;
            if (GlobalConfig.with_local_time_stepping) { dtInv = 1.0/cell.dt_local; }
            else { dtInv = 1.0/dt; }
            number tmp = 0.0;
            foreach (f; cell.iface) {
                tmp += f.area[0]*spectral_radius(f, GlobalConfig.viscous, omega);
            }
            // diagonal scalar
            cell.D = cell.volume[0]*dtInv + 0.5*tmp;

            cell.dU[0].mass = (1.0/cell.D)*cell.dUdt[0].mass*cell.volume[0];
            cell.dU[0].momentum.refx = (1.0/cell.D)*cell.dUdt[0].momentum.x*cell.volume[0];
            cell.dU[0].momentum.refy = (1.0/cell.D)*cell.dUdt[0].momentum.y*cell.volume[0];
            if ( blk.myConfig.dimensions == 3 )
                cell.dU[0].momentum.refz = (1.0/cell.D)*cell.dUdt[0].momentum.z*cell.volume[0];
            cell.dU[0].total_energy = (1.0/cell.D)*cell.dUdt[0].total_energy*cell.volume[0];
        }
    }
    
    // 3. perform kmax subiterations    
    foreach (k; 0 .. kmax) {
        // exchange dU values
        exchange_ghost_cell_boundary_data(pseudoSimTime, 0, 0);
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (cell; blk.cells) {
                number LU_mass = 0.0;
                number LU_momentumx = 0.0; number LU_momentumy = 0.0; number LU_momentumz = 0.0;
                number LU_totalenergy = 0.0;
                // approximate off-diagonal terms (L+U)
                foreach (i; 1..cell.cell_cloud.length) {
                    FVCell ncell = cell.cell_cloud[i]; 
                    FVInterface f = cell.iface[i-1];
                    number lij = spectral_radius(f, GlobalConfig.viscous, omega);
                    
                    evalMatrixFreeFluxIncrement(ncell, f);
                    
                    LU_mass += (ncell.dF.mass*cell.outsign[i-1] - lij*ncell.dU[0].mass)*f.area[0];
                    LU_momentumx += (ncell.dF.momentum.x*cell.outsign[i-1] - lij*ncell.dU[0].momentum.x)*f.area[0];
                    LU_momentumy += (ncell.dF.momentum.y*cell.outsign[i-1] - lij*ncell.dU[0].momentum.y)*f.area[0];
                    if ( blk.myConfig.dimensions == 3 )
                        LU_momentumz += (ncell.dF.momentum.z*cell.outsign[i-1] - lij*ncell.dU[0].momentum.z)*f.area[0];
                    LU_totalenergy += (ncell.dF.total_energy*cell.outsign[i-1] - lij*ncell.dU[0].total_energy)*f.area[0];
                }
                // update dU
                cell.dU[1].mass = (1.0/cell.D) * (cell.dUdt[0].mass*cell.volume[0] - 0.5*LU_mass);
                cell.dU[1].momentum.refx = (1.0/cell.D) * (cell.dUdt[0].momentum.x*cell.volume[0] - 0.5*LU_momentumx);
                cell.dU[1].momentum.refy = (1.0/cell.D) * (cell.dUdt[0].momentum.y*cell.volume[0] - 0.5*LU_momentumy);
                if ( blk.myConfig.dimensions == 3 )
                    cell.dU[1].momentum.refz = (1.0/cell.D) * (cell.dUdt[0].momentum.z*cell.volume[0] - 0.5*LU_momentumz);
                cell.dU[1].total_energy = (1.0/cell.D) * (cell.dUdt[0].total_energy*cell.volume[0] - 0.5*LU_totalenergy);                
            }
        }
        // update most recent values of dU
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (cell; blk.cells) {
                cell.dU[0].copy_values_from(cell.dU[1]);
            }
        }
    }

    // pack up dU vector
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;

    foreach (blk; parallel(localFluidBlocks,1)) {
        int cellCount = 0;
        foreach (i, cell; blk.cells) {
            // used later to update conserved quantities
            blk.dU[cellCount+MASS] = cell.dU[0].mass;
            blk.dU[cellCount+X_MOM] = cell.dU[0].momentum.x;
            blk.dU[cellCount+Y_MOM] = cell.dU[0].momentum.y;
            if ( blk.myConfig.dimensions == 3 )
                blk.dU[cellCount+Z_MOM] = cell.dU[0].momentum.z;
            blk.dU[cellCount+TOT_ENERGY] = cell.dU[0].total_energy;

            // used later to compute residual
            blk.FU[cellCount+MASS] = -cell.dUdt[0].mass;
            blk.FU[cellCount+X_MOM] = -cell.dUdt[0].momentum.x;
            blk.FU[cellCount+Y_MOM] = -cell.dUdt[0].momentum.y;
            if ( GlobalConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = -cell.dUdt[0].momentum.z;
            blk.FU[cellCount+TOT_ENERGY] = -cell.dUdt[0].total_energy;

            cellCount += nConserved;
        }
    }

    // calculate residual
    mixin(dot_over_blocks("residual", "FU", "FU"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(residual), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    residual = sqrt(residual);
    
} // end lusgs_solve()

void rpcGMRES_solve(int step, double pseudoSimTime, double dt, double eta, double sigma, bool usePreconditioner,
                    ref double residual, ref int nRestarts, int startStep)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;

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
        int cellCount = 0;
        blk.maxRate.mass = 0.0;
        blk.maxRate.momentum.refx = 0.0;
        blk.maxRate.momentum.refy = 0.0;
        if ( blk.myConfig.dimensions == 3 )
            blk.maxRate.momentum.refz = 0.0;
        blk.maxRate.total_energy = 0.0;
        foreach(it; 0 .. nturb){
            blk.maxRate.rhoturb[it] = 0.0;
        }
        foreach (i, cell; blk.cells) {
            blk.FU[cellCount+MASS] = -cell.dUdt[0].mass;
            blk.FU[cellCount+X_MOM] = -cell.dUdt[0].momentum.x;
            blk.FU[cellCount+Y_MOM] = -cell.dUdt[0].momentum.y;
            if ( blk.myConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = -cell.dUdt[0].momentum.z;
            blk.FU[cellCount+TOT_ENERGY] = -cell.dUdt[0].total_energy;
            foreach(it; 0 .. nturb){
                blk.FU[cellCount+TKE+it] = -cell.dUdt[0].rhoturb[it];
            }
            cellCount += nConserved;
            /*
            if (blk.id == 0) {
                writefln("i= %d, dUdt.mass= %e", i, cell.dUdt[0].mass.re);
            }
            */
            blk.maxRate.mass = fmax(blk.maxRate.mass, fabs(cell.dUdt[0].mass));
            /*
            if (blk.id == 0) {
                writefln("i= %d, maxRate.mass= %e", i, blk.maxRate.mass.re);
            }
            */
            blk.maxRate.momentum.refx = fmax(blk.maxRate.momentum.x, fabs(cell.dUdt[0].momentum.x));
            blk.maxRate.momentum.refy = fmax(blk.maxRate.momentum.y, fabs(cell.dUdt[0].momentum.y));
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.momentum.refz = fmax(blk.maxRate.momentum.z, fabs(cell.dUdt[0].momentum.z));
            blk.maxRate.total_energy = fmax(blk.maxRate.total_energy, fabs(cell.dUdt[0].total_energy));
            foreach(it; 0 .. nturb){
                blk.maxRate.rhoturb[it] = fmax(blk.maxRate.rhoturb[it], fabs(cell.dUdt[0].rhoturb[it]));
            }
        }
    }

    number maxMass = 0.0;
    number maxMomX = 0.0;
    number maxMomY = 0.0;
    number maxMomZ = 0.0;
    number maxEnergy = 0.0;
    number[2] maxTurb; maxTurb[0] = 0.0; maxTurb[1] = 0.0;
    foreach (blk; localFluidBlocks) {
        maxMass = fmax(maxMass, blk.maxRate.mass);
        maxMomX = fmax(maxMomX, blk.maxRate.momentum.x);
        maxMomY = fmax(maxMomY, blk.maxRate.momentum.y);
        if ( blk.myConfig.dimensions == 3 )
            maxMomZ = fmax(maxMomZ, blk.maxRate.momentum.z);
        maxEnergy = fmax(maxEnergy, blk.maxRate.total_energy);
        foreach(it; 0 .. blk.myConfig.turb_model.nturb){
            maxTurb[it] = fmax(maxTurb[it], blk.maxRate.rhoturb[it]);
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

    // Get a copy of the maxes out to each block
    foreach (blk; parallel(localFluidBlocks,1)) {
        if (blk.myConfig.sssOptions.useScaling) {
            blk.maxRate.mass = maxMass;
            blk.maxRate.momentum.refx = maxMomX;
            blk.maxRate.momentum.refy = maxMomY;
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.momentum.refz = maxMomZ;
            blk.maxRate.total_energy = maxEnergy;
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                blk.maxRate.rhoturb[it] = maxTurb[it];
            }
        }
        else { // just scale by 1
            blk.maxRate.mass = 1.0;
            blk.maxRate.momentum.refx = 1.0;
            blk.maxRate.momentum.refy = 1.0;
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.momentum.refz = 1.0;
            blk.maxRate.total_energy = 1.0;
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                blk.maxRate.rhoturb[it] = 1.0; 
            }
        }
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.FU[cellCount+MASS] = -blk.FU[cellCount+MASS];
            blk.FU[cellCount+X_MOM] = -blk.FU[cellCount+X_MOM];
            blk.FU[cellCount+Y_MOM] = -blk.FU[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = -blk.FU[cellCount+Z_MOM];
            blk.FU[cellCount+TOT_ENERGY] = -blk.FU[cellCount+TOT_ENERGY];
            foreach(it; 0 .. nturb){
                number fac = blk.myConfig.turb_model.gmres_scaling_factor(it);
                blk.FU[cellCount+TKE+it] = -1.0/fac*blk.FU[cellCount+TKE+it];
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

    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.FU[cellCount+MASS] = -blk.FU[cellCount+MASS];
            blk.FU[cellCount+X_MOM] = -blk.FU[cellCount+X_MOM];
            blk.FU[cellCount+Y_MOM] = -blk.FU[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = -blk.FU[cellCount+Z_MOM];
            blk.FU[cellCount+TOT_ENERGY] = -blk.FU[cellCount+TOT_ENERGY];
            foreach(it; 0 .. nturb){
                number fac = blk.myConfig.turb_model.gmres_scaling_factor(it);
                blk.FU[cellCount+TKE+it] = -fac*blk.FU[cellCount+TKE+it];
            }
            cellCount += nConserved;
        }
    }

    // Initialise some arrays and matrices that have already been allocated
    g0[] = to!number(0.0);
    g1[] = to!number(0.0);
    H0.zeros();
    H1.zeros();

    //double dtInv = 1.0/dt;

    // We'll scale r0 against these max rates of change.
    // r0 = b - A*x0
    // Taking x0 = [0] (as is common) gives r0 = b = FU
    // apply scaling
    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.r0[cellCount+MASS] = -(1./blk.maxRate.mass)*blk.FU[cellCount+MASS];
            blk.r0[cellCount+X_MOM] = -(1./blk.maxRate.momentum.x)*blk.FU[cellCount+X_MOM];
            blk.r0[cellCount+Y_MOM] = -(1./blk.maxRate.momentum.y)*blk.FU[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                blk.r0[cellCount+Z_MOM] = -(1./blk.maxRate.momentum.z)*blk.FU[cellCount+Z_MOM];
            blk.r0[cellCount+TOT_ENERGY] = -(1./blk.maxRate.total_energy)*blk.FU[cellCount+TOT_ENERGY];
            foreach(it; 0 .. nturb){
                blk.r0[cellCount+TKE+it] = -(1./blk.maxRate.rhoturb[it])*blk.FU[cellCount+TKE+it];
            }
            cellCount += nConserved;
        }
    }

    // Then compute v = r0/||r0||
    number beta;
    mixin(dot_over_blocks("beta", "r0", "r0"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    beta = sqrt(beta);
    g0[0] = beta;
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) {
            blk.v[k] = blk.r0[k]/beta;
            blk.V[k,0] = blk.v[k];
        }
    }

    // Compute tolerance
    auto outerTol = eta*beta;

    // 2. Start outer-loop of restarted GMRES
    for ( r = 0; r < maxRestarts; r++ ) {
        // 2a. Begin iterations
        foreach (j; 0 .. m) {
            iterCount = j+1;

            // apply scaling
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.turb_model.nturb;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.v[cellCount+MASS] *= (blk.maxRate.mass);
                    blk.v[cellCount+X_MOM] *= (blk.maxRate.momentum.x);
                    blk.v[cellCount+Y_MOM] *= (blk.maxRate.momentum.y);
                    if ( blk.myConfig.dimensions == 3 )
                        blk.v[cellCount+Z_MOM] *= (blk.maxRate.momentum.z);
                    blk.v[cellCount+TOT_ENERGY] *= (blk.maxRate.total_energy);
                    foreach(it; 0 .. nturb){
                        blk.v[cellCount+TKE+it] *= (blk.maxRate.rhoturb[it]);
                    }
                    cellCount += nConserved;
                }
            }
            
            // apply preconditioning
            if (usePreconditioner && step >= GlobalConfig.sssOptions.startPreconditioning) {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    int n = blk.myConfig.sssOptions.frozenPreconditionerCount; //GlobalConfig.sssOptions.frozenPreconditionerCount;
                    // We compute the precondition matrix on the very first step after the start up steps
                    // We then only update the precondition matrix once per GMRES call on every nth flow solver step. 
                    if (r == 0 && j == 0 && (step == blk.myConfig.sssOptions.startPreconditioning ||
                                             step%n == 0 ||
                                             step == startStep ||
                                             step == blk.myConfig.sssOptions.nStartUpSteps+1))
                        sss_preconditioner(blk, nConserved, dt);
                    final switch (blk.myConfig.sssOptions.preconditionMatrixType) {
                        case PreconditionMatrixType.block_diagonal:
                            int cellCount = 0;
                            number[] tmp;
                            tmp.length = nConserved;
                            foreach (cell; blk.cells) { 
                                nm.bbla.dot(cell.dConservative, blk.v[cellCount..cellCount+nConserved], tmp);
                                blk.zed[cellCount..cellCount+nConserved] = tmp[];
                                cellCount += nConserved;
                            }
                            break;
                        case PreconditionMatrixType.ilu:
                            blk.zed[] = blk.v[];
			    nm.smla.transpose_solve(blk.P, blk.zed);
			    break;
                    } // end switch
                }
            }
            else {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    blk.zed[] = blk.v[];
                }
            }
            
            // Prepare 'w' with (I/dt)(P^-1)v term;
            //foreach (blk; parallel(localFluidBlocks,1)) {
                //double dtInv = 1.0/dt;
                //foreach (idx; 0..blk.w.length) blk.w[idx] = dtInv*blk.zed[idx];
            //}

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
            evalJacobianVecProd(pseudoSimTime, sigma);

            // Now we can complete calculation of w
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars)  blk.w[k] += blk.zed[k];
            }
            
            // apply scaling
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.turb_model.nturb;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.w[cellCount+MASS] *= (1./blk.maxRate.mass);
                    blk.w[cellCount+X_MOM] *= (1./blk.maxRate.momentum.x);
                    blk.w[cellCount+Y_MOM] *= (1./blk.maxRate.momentum.y);
                    if ( blk.myConfig.dimensions == 3 )
                        blk.w[cellCount+Z_MOM] *= (1./blk.maxRate.momentum.z);
                    blk.w[cellCount+TOT_ENERGY] *= (1./blk.maxRate.total_energy);
                    foreach(it; 0 .. nturb){
                        blk.w[cellCount+TKE+it] *= (1./blk.maxRate.rhoturb[it]);
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
                    MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
                MPI_Allreduce(MPI_IN_PLACE, &(H0_jp1j.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
            int cellCount = 0;
            foreach (cell; blk.cells) {
                blk.zed[cellCount+MASS] *= (blk.maxRate.mass);
                blk.zed[cellCount+X_MOM] *= (blk.maxRate.momentum.x);
                blk.zed[cellCount+Y_MOM] *= (blk.maxRate.momentum.y);
                if ( blk.myConfig.dimensions == 3 )
                    blk.zed[cellCount+Z_MOM] *= (blk.maxRate.momentum.z);
                blk.zed[cellCount+TOT_ENERGY] *= (blk.maxRate.total_energy);
                foreach(it; 0 .. nturb){
                    blk.zed[cellCount+TKE+it] *= (blk.maxRate.rhoturb[it]);
                }
                cellCount += nConserved;
            }
        }

        // apply preconditioning
        if (usePreconditioner && step >= GlobalConfig.sssOptions.startPreconditioning) {
            foreach(blk; parallel(localFluidBlocks,1)) {
                final switch (blk.myConfig.sssOptions.preconditionMatrixType) {
                case PreconditionMatrixType.block_diagonal:
                    int cellCount = 0;
                    number[] tmp;
                    tmp.length = nConserved;
                    foreach (cell; blk.cells) {
                        nm.bbla.dot(cell.dConservative, blk.zed[cellCount..cellCount+nConserved], tmp);
                        blk.dU[cellCount..cellCount+nConserved] = tmp[];
                        cellCount += nConserved;
                    }
                    break;
                case PreconditionMatrixType.ilu:
                    blk.dU[] = blk.zed[];
		    nm.smla.transpose_solve(blk.P, blk.dU);
                    break;
                } // end switch
            }
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

        // apply scaling
        foreach (blk; parallel(localFluidBlocks,1)) {
            size_t nturb = blk.myConfig.turb_model.nturb;
            int cellCount = 0;
            foreach (cell; blk.cells) {
                blk.r0[cellCount+MASS] *= (1.0/blk.maxRate.mass);
                blk.r0[cellCount+X_MOM] *= (1.0/blk.maxRate.momentum.x);
                blk.r0[cellCount+Y_MOM] *= (1.0/blk.maxRate.momentum.y);
                if ( blk.myConfig.dimensions == 3 )
                    blk.r0[cellCount+Z_MOM] *= (1.0/blk.maxRate.momentum.z);
                blk.r0[cellCount+TOT_ENERGY] *= (1.0/blk.maxRate.total_energy);
                foreach(it; 0 .. nturb){
                    blk.r0[cellCount+TKE+it] *= (1.0/blk.maxRate.rhoturb[it]);
                }
                cellCount += nConserved;
            }
        }

        mixin(dot_over_blocks("beta", "r0", "r0"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
    size_t nConserved = nConservedQuantities;

    foreach (blk; parallel(localFluidBlocks,1)) {
        size_t nturb = blk.myConfig.turb_model.nturb;
        blk.residuals.copy_values_from(blk.cells[0].dUdt[0]);
        blk.residuals.mass = fabs(blk.residuals.mass);
        blk.residuals.momentum.refx = fabs(blk.residuals.momentum.x);
        blk.residuals.momentum.refy = fabs(blk.residuals.momentum.y);
        if ( blk.myConfig.dimensions == 3 )
            blk.residuals.momentum.refz = fabs(blk.residuals.momentum.z);
        blk.residuals.total_energy = fabs(blk.residuals.total_energy);
        foreach(it; 0 .. nturb){
            blk.residuals.rhoturb[it] = fabs(blk.residuals.rhoturb[it]);
        }
        number massLocal, xMomLocal, yMomLocal, zMomLocal, energyLocal;
        number[2] turbLocal;
        foreach (cell; blk.cells) {
            massLocal = cell.dUdt[0].mass;
            xMomLocal = cell.dUdt[0].momentum.x;
            yMomLocal = cell.dUdt[0].momentum.y;
            zMomLocal = cell.dUdt[0].momentum.z;
            energyLocal = cell.dUdt[0].total_energy;
            foreach(it; 0 .. nturb){
                turbLocal[it] = cell.dUdt[0].rhoturb[it];
            }
            blk.residuals.mass = fmax(blk.residuals.mass, massLocal);
            blk.residuals.momentum.refx = fmax(blk.residuals.momentum.x, xMomLocal);
            blk.residuals.momentum.refy = fmax(blk.residuals.momentum.y, yMomLocal);
            if ( blk.myConfig.dimensions == 3 )
                blk.residuals.momentum.refz = fmax(blk.residuals.momentum.z, zMomLocal);
            blk.residuals.total_energy = fmax(blk.residuals.total_energy, energyLocal);
            foreach(it; 0 .. nturb){
                blk.residuals.rhoturb[it] = fmax(blk.residuals.rhoturb[it], turbLocal[it]);
            }
        }
    }
    residuals.copy_values_from(localFluidBlocks[0].residuals);
    foreach (blk; localFluidBlocks) {
        residuals.mass = fmax(residuals.mass, blk.residuals.mass);
        residuals.momentum.refx = fmax(residuals.momentum.x, blk.residuals.momentum.x);
        residuals.momentum.refy = fmax(residuals.momentum.y, blk.residuals.momentum.y);
        if ( blk.myConfig.dimensions == 3 )
            residuals.momentum.refz = fmax(residuals.momentum.z, blk.residuals.momentum.z);
        residuals.total_energy = fmax(residuals.total_energy, blk.residuals.total_energy);
        foreach(it; 0 .. blk.myConfig.turb_model.nturb){
            residuals.rhoturb[it] = fmax(residuals.rhoturb[it], blk.residuals.rhoturb[it]);
        }
    }
    version(mpi_parallel) {
        // In MPI context, only the master task (rank 0) has collated the residuals
        double maxResid;
        MPI_Reduce(&(residuals.mass.re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.mass.re = maxResid; }
        MPI_Reduce(&(residuals.momentum.x.re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.momentum.refx = to!number(maxResid); }
        MPI_Reduce(&(residuals.momentum.y.re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.momentum.refy = to!number(maxResid); }
        if (GlobalConfig.dimensions == 3) {
            MPI_Reduce(&(residuals.momentum.z.re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals.momentum.refz = to!number(maxResid); }
        }
        MPI_Reduce(&(residuals.total_energy.re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals.total_energy.re = maxResid; }
        foreach (it; 0 .. GlobalConfig.turb_model.nturb) {
            MPI_Reduce(&(residuals.rhoturb[it].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals.rhoturb[it].re = maxResid; }
        }
    }
}

void rewrite_times_file(RestartInfo[] times)
{
    auto fname = format("./config/%s.times", GlobalConfig.base_file_name);
    auto f = File(fname, "w");
    f.writeln("# tindx sim_time dt_global");
    foreach (i, rInfo; times) {
        if ( GlobalConfig.dimensions == 2 ) {
            f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e",
                     i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                     rInfo.globalResidual.re, rInfo.residuals.mass.re,
                     rInfo.residuals.momentum.x.re, rInfo.residuals.momentum.y.re,
                     rInfo.residuals.total_energy.re);
        }
        else {
            f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e %.18e",
                     i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                     rInfo.globalResidual.re, rInfo.residuals.mass.re,
                     rInfo.residuals.momentum.x.re, rInfo.residuals.momentum.y.re, rInfo.residuals.momentum.z.re,
                     rInfo.residuals.total_energy.re);
        }
        foreach(it; 0 .. GlobalConfig.turb_model.nturb){
            f.writef(" %.18e", rInfo.residuals.rhoturb[it].re);
        }
        f.write("\n");
    }
    f.close();
}
