/** steadystate_core.d
 * Core set of functions used in the Newton-Krylov updates for steady-state convergence.
 *
 * Author: Rowan G.
 * Date: 2016-10-09
 */

module steadystate_core;

import core.memory;
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
double[] g0;
double[] g1;
double[] h;
double[] hR;
Matrix!double H0;
Matrix!double H1;
Matrix!double Gamma;
Matrix!double Q0;
Matrix!double Q1;

struct RestartInfo {
    double pseudoSimTime;
    double dt;
    double cfl;
    int step;
    double globalResidual;
    ConservedQuantities residuals;

    this(size_t n)
    {
        residuals = new_ConservedQuantities(n);
    }
}

void extractRestartInfoFromTimesFile(string jobName, ref RestartInfo[] times)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    // remove the conserved mass variable for multi-species gas
    //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
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
            if ( GlobalConfig.gmodel_master.n_species == 1 ) { restartInfo.residuals[cqi.mass] = to!double(tokens[6+MASS]); }
            restartInfo.residuals[cqi.xMom] = to!double(tokens[6+X_MOM]);
            restartInfo.residuals[cqi.yMom] = to!double(tokens[6+Y_MOM]);
            if ( GlobalConfig.dimensions == 3 )
                restartInfo.residuals[cqi.zMom] = to!double(tokens[6+Z_MOM]);
            restartInfo.residuals[cqi.totEnergy] = to!double(tokens[6+TOT_ENERGY]);
            foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
                restartInfo.residuals[cqi.rhoturb+it] = to!double(tokens[6+TKE+it]);
            }
            version(multi_species_gas){
            foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                restartInfo.residuals[cqi.species+sp] = to!double(tokens[6+SPECIES+sp]);
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                restartInfo.residuals[cqi.modes+imode] = to!double(tokens[6+MODES+imode]);
            }
            }
            times ~= restartInfo;
        }
        line = timesFile.readln().strip();
    }
    timesFile.close();
    return;
}

double extractReferenceResidualsFromFile(string jobName, ref ConservedQuantities maxResiduals)
{
/*
    When restarting in the middle of a calculation, we need to read in the reference
    residuals from where they were saved from the original startup. This function
    does the, slightly sketchy, double duty of also setting normRef via its return
    argument.

    @author: Nick Gibbons
*/
    // Make a stack-local copy of conserved quantities info
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;
    size_t SPECIES = GlobalConfig.cqi.species;
    size_t MODES = GlobalConfig.cqi.modes;
    auto cqi = GlobalConfig.cqi;

    // We need to read in the reference residual values from a file.
    string refResidFname = jobName ~ "-ref-residuals.saved";
    auto refResid = File(refResidFname, "r");
    auto line = refResid.readln().strip();
    auto tokens = line.split();
    if (tokens.length==0) throw new Error(format("Error reading %s, no entries found!", refResidFname));

    double normRef = to!double(tokens[0]);
    if ( GlobalConfig.gmodel_master.n_species == 1 ) { maxResiduals[cqi.mass] = to!double(tokens[1+MASS]); }
    maxResiduals[cqi.xMom] = to!double(tokens[1+X_MOM]);
    maxResiduals[cqi.yMom] = to!double(tokens[1+Y_MOM]);
    if ( GlobalConfig.dimensions == 3 )
        maxResiduals[cqi.zMom] = to!double(tokens[1+Z_MOM]);
    maxResiduals[cqi.totEnergy] = to!double(tokens[1+TOT_ENERGY]);
    foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
        maxResiduals[cqi.rhoturb+it] = to!double(tokens[1+TKE+it]);
    }
    version(multi_species_gas){
        if ( GlobalConfig.gmodel_master.n_species > 1 ) {
            foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                maxResiduals[cqi.species+sp] = to!double(tokens[1+SPECIES+sp]);
            }
        }
    }
    version(multi_T_gas){
        foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
            maxResiduals[cqi.modes+imode] = to!double(tokens[1+MODES+imode]);
        }
    }
    return normRef;
}

version(mpi_parallel){

void broadcastResiduals(ref ConservedQuantities residuals){
/*
    Helper function for sending around a conserved quantities object.

    @author: Nick Gibbons
*/
    double[] buffer;

    int size = to!int(residuals.length);
    version(complex_numbers){ size *= 2; }
    buffer.length = size;

    size_t i=0;
    foreach(f; residuals){
        buffer[i] = f.re;
        i++;
        version(complex_numbers){
            buffer[i] = f.im;
            i++;
        }
    }

    MPI_Bcast(buffer.ptr, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    i=0;
    foreach(j; 0 .. residuals.length){
        residuals[j].re = buffer[i];
        i++;
        version(complex_numbers){
            residuals[j].im = buffer[i];
            i++;
        }
    }
    return;
}

void broadcastRestartInfo(ref RestartInfo[] times){
/*
    Helper function for sending around a collection of RestartInfo objects

    Notes: For some reason, this function does need a ref in its argument,
    even though dynamic arrays are supposed to be reference types...

    @author: Nick Gibbons
*/
    // The main issue here is that the non-master processes do not actually know
    // how many of these things there are. First we need sort that out.
    // FIXME: Note that this could cause problems if this function is reused.
    int ntimes = to!int(times.length);
    MPI_Bcast(&ntimes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!GlobalConfig.is_master_task) {
        foreach(i; 0 .. ntimes){
            times ~= RestartInfo(GlobalConfig.cqi.n);
        }
    }

    // Now that we have them allocated, we can fill them out
    foreach(n; 0 .. ntimes){
        MPI_Bcast(&(times[n].pseudoSimTime), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(times[n].dt), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(times[n].cfl), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(times[n].step), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(times[n].globalResidual), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        broadcastResiduals(times[n].residuals);
    }
}

} // end version(mpi_parallel)

void iterate_to_steady_state(int snapshotStart, int maxCPUs, int threadsPerMPITask, bool include_solid_domain = true, bool init_precondition_matrix = true)
{
    auto wallClockStart = Clock.currTime();
    string jobName = GlobalConfig.base_file_name;
    int nsteps = GlobalConfig.sssOptions.nTotalSteps;
    int nIters = 0;
    int nRestarts;
    double linSolResid = 0;
    double relGlobalResidReduction = GlobalConfig.sssOptions.stopOnRelGlobalResid;
    double absGlobalResidReduction = GlobalConfig.sssOptions.stopOnAbsGlobalResid;
    double massBalanceReduction = GlobalConfig.sssOptions.stopOnMassBalance;
    double cfl_max = GlobalConfig.sssOptions.cfl_max;
    double cfl_min = GlobalConfig.sssOptions.cfl_min;
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
    else if (GlobalConfig.interpolation_order > 2) { LHSeval1 = GlobalConfig.interpolation_order; }
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
    // remove the conserved mass variable for multi-species gas
    //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
    immutable size_t MASS = GlobalConfig.cqi.mass;
    immutable size_t X_MOM = GlobalConfig.cqi.xMom;
    immutable size_t Y_MOM = GlobalConfig.cqi.yMom;
    immutable size_t Z_MOM = GlobalConfig.cqi.zMom;
    immutable size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    immutable size_t TKE = GlobalConfig.cqi.rhoturb;
    immutable size_t SPECIES = GlobalConfig.cqi.species;
    immutable size_t MODES = GlobalConfig.cqi.modes;
    immutable size_t nftl = GlobalConfig.n_flow_time_levels;
    immutable size_t ncq  = nConserved; // number of conserved quantities

    ConservedQuantities maxResiduals = new_ConservedQuantities(nConserved);
    ConservedQuantities currResiduals = new_ConservedQuantities(nConserved);
    number mass_balance = 0.0;
    double omega = 1.0;                  // nonlinear update relaxation factor
    double omega_allow_cfl_growth = 0.1; // we freeze the CFL for a relaxation factor below this value
    double omega_min = 0.01;             // minimum allowable relaxation factor to accept an update
    double omega_reduction_factor = 0.7; // factor by which the relaxation factor is reduced during
                                         // the physicality check and line search
    double theta = GlobalConfig.sssOptions.physicalityCheckTheta;
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
        if (GlobalConfig.verbosity_level>1){
            writefln("MPI-task %d : running with %d threads.", GlobalConfig.mpi_rank_for_local_task, extraThreadsInPool+1);
        }
    }
    else {
        writefln("Single process running with %d threads.", extraThreadsInPool+1); // +1 for main thread.
    }
    double normRef = 0.0;
    bool residualsUpToDate = false;
    bool finalStep = false;
    bool usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
    if (usePreconditioner && init_precondition_matrix) {
        evalRHS(0.0, 0);
        // initialize the flow Jacobians used as local precondition matrices for GMRES
        final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
            case PreconditionMatrixType.diagonal:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(-1, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.jacobi:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.sgs:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
            case PreconditionMatrixType.ilu:
                foreach (blk; localFluidBlocks) { blk.initialize_jacobian(0, GlobalConfig.sssOptions.preconditionerSigma); }
                break;
        } // end switch
        //foreach (blk; localFluidBlocks) { blk.verify_jacobian(GlobalConfig.sssOptions.preconditionerSigma); }
    }

    // We need to calculate the initial residual if we are starting from scratch.
    if ( snapshotStart == 0 ) {
        evalRHS(0.0, 0);
        max_residuals(maxResiduals);
        foreach (blk; parallel(localFluidBlocks, 1)) {
            foreach (i; 0 .. blk.ncells*ncq) {
                blk.FU[i] = -blk.celldata.dUdt0[i].re;
            }
        }

        if (GlobalConfig.sssOptions.include_turb_quantities_in_residual == false) {
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.cqi.n_turb;
                int cellCount = 0;
                foreach (i; 0 .. blk.ncells) {
                    foreach(it; 0 .. nturb) { blk.FU[cellCount+TKE+it] = 0.0; }
                    cellCount += nConserved;
                }
            }
        }

        mixin(dot_over_blocks("normRef", "FU", "FU"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(normRef), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        normRef = sqrt(normRef);

        if (GlobalConfig.sssOptions.include_turb_quantities_in_residual == false) {
            foreach (blk; parallel(localFluidBlocks,1)) {
                size_t nturb = blk.myConfig.cqi.n_turb;
                int cellCount = 0;
                foreach (i; 0 .. blk.ncells) {
                    size_t s0 = i*ncq;
                    foreach(it; 0 .. nturb) { blk.FU[cellCount+TKE+it] = -blk.celldata.dUdt0[s0+TKE+it].re; }
                    cellCount += nConserved;
                }
            }
        }
    }

    // if we are restarting a simulation we need read the initial residuals from a file
    RestartInfo[] times;
    if (snapshotStart > 0) {
        if (GlobalConfig.is_master_task){
            extractRestartInfoFromTimesFile(jobName, times);
            normRef = extractReferenceResidualsFromFile(jobName, maxResiduals);
        }
        // Only the master task reads from disk, in MPI we need to send the info around
        version(mpi_parallel){
            MPI_Bcast(&normRef, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            broadcastRestartInfo(times);
            broadcastResiduals(maxResiduals);
        }
        normOld = times[snapshotStart].globalResidual;

        // We also need to determine how many snapshots have already been written
        // We don't count the initial solution as a written snapshot
        nWrittenSnapshots = to!int(times.length) - 1;

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
        if ( GlobalConfig.gmodel_master.n_species == 1 ) { writefln("MASS:           %.12e", maxResiduals[cqi.mass].re); }
        writefln("X-MOMENTUM:     %.12e", maxResiduals[cqi.xMom].re);
        writefln("Y-MOMENTUM:     %.12e", maxResiduals[cqi.yMom].re);
        if ( GlobalConfig.dimensions == 3 )
            writefln("Z-MOMENTUM:     %.12e", maxResiduals[cqi.zMom].re);
        writefln("ENERGY:         %.12e", maxResiduals[cqi.totEnergy].re);
        foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
            string tvname = capitalize(GlobalConfig.turb_model.primitive_variable_name(it));
            writefln("%s:            %.12e",tvname, maxResiduals[cqi.rhoturb+it].re);
        }
        version(multi_species_gas){
            if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                    string spname = capitalize(GlobalConfig.gmodel_master.species_name(sp));
                    writefln("%s:            %.12e",spname, maxResiduals[cqi.species+sp].re);
                }
            }
        }
        version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                string modename = "T_MODES["~to!string(imode)~"]"; //capitalize(GlobalConfig.gmodel_master.energy_mode_name(imode));
                writefln("%s:            %.12e",modename, maxResiduals[cqi.modes+imode].re);
            }
        }
        if (GlobalConfig.sssOptions.include_turb_quantities_in_residual == false) {
            writeln("WARNING: The GLOBAL reference residual does not include turbulence quantities since include_turb_quantities_in_residual is set to false.");
        }
        // store the initial residuals
        string refResidFname = jobName ~ "-ref-residuals.saved";
        auto refResid = File(refResidFname, "w");
        if ( GlobalConfig.dimensions == 2 ) {
            if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                refResid.writef("%.18e %.18e %.18e %.18e %.18e",
                                normRef, maxResiduals[cqi.mass].re, maxResiduals[cqi.xMom].re,
                                maxResiduals[cqi.yMom].re, maxResiduals[cqi.totEnergy].re);
            } else {
                refResid.writef("%.18e %.18e %.18e %.18e",
                                normRef, maxResiduals[cqi.xMom].re,
                                maxResiduals[cqi.yMom].re, maxResiduals[cqi.totEnergy].re);
            }
        } else {
            if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                refResid.writef("%.18e %.18e %.18e %.18e %.18e %.18e",
                                normRef, maxResiduals[cqi.mass].re, maxResiduals[cqi.xMom].re,
                                maxResiduals[cqi.yMom].re, maxResiduals[cqi.zMom].re,
                                maxResiduals[cqi.totEnergy].re);
            } else {
                refResid.writef("%.18e %.18e %.18e %.18e %.18e",
                                normRef, maxResiduals[cqi.xMom].re,
                                maxResiduals[cqi.yMom].re, maxResiduals[cqi.zMom].re,
                                maxResiduals[cqi.totEnergy].re);
            }
        }
        foreach(it; 0 .. GlobalConfig.turb_model.nturb) {
            refResid.writef(" %.18e", maxResiduals[cqi.rhoturb+it].re);
        }
        version(multi_species_gas){
            if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                    refResid.writef(" %.18e", maxResiduals[cqi.species+sp].re);
                }
            }
        }
        version(multi_T_gas){
            foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
                refResid.writef(" %.18e", maxResiduals[cqi.modes+imode].re);
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
        cfl = (nStartUpSteps==0) ? cfl1 : cfl0;
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
            if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                fResid.writefln("#  %02d: mass-abs", 12+2*MASS);
                fResid.writefln("# %02d: mass-rel", 12+2*MASS+1);
            }
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
            if (step > GlobalConfig.shock_detector_freeze_step) {
                GlobalConfig.frozen_shock_detector = true;
            }
        }

        // solve linear system for dU
        version(pir) { point_implicit_relaxation_solve(step, pseudoSimTime, dt, normNew, linSolResid, startStep); }
        else { rpcGMRES_solve(step, pseudoSimTime, dt, eta, sigma, usePreconditioner, normNew, nRestarts, nIters, linSolResid, startStep, LHSeval, RHSeval, pc_matrix_evaluated); }

        // calculate a relaxation factor for the nonlinear update via a physicality check
        omega = 1.0; // start by allowing a full update
        if (GlobalConfig.sssOptions.usePhysicalityCheck) {

            // we first limit the change in the conserved mass
            foreach (blk; localFluidBlocks) {
                auto cqi = blk.myConfig.cqi;
                int cellCount = 0;
                number rel_diff_limit, U, dU;
                foreach (i; 0 .. blk.ncells) {
                    if (cqi.n_species == 1) {
                        U  = blk.celldata.U0[cellCount+MASS];
                        dU = blk.dU[cellCount+MASS];
                    } else { // sum up the species densities
                        U = 0.0;
                        dU = 0.0;
                        foreach(isp; 0 .. cqi.n_species) {
                            U += blk.celldata.U0[cellCount+SPECIES+isp];
                            dU += blk.dU[cellCount+SPECIES+isp];
                        }
                    }
                    rel_diff_limit = fabs(dU/(theta*U));
                    omega = 1.0/(fmax(rel_diff_limit.re,1.0/omega));
                    cellCount += nConserved;
                }
            }

            // communicate minimum relaxation factor based on conserved mass to all processes
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &(omega), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }

            // we now check if the current relaxation factor is sufficient enough to produce realizable primitive variables
            // if it isn't, we reduce the relaxation factor by some prescribed factor and then try again
            foreach (blk; localFluidBlocks) {
                int cellCount = 0;
                foreach (i; 0 .. blk.ncells) {
                    bool failed_decode = false;
                    blk.fs_save.copy_values_from(blk.celldata.flowstates[i]);
                    while (omega >= omega_min) {
                        // check positivity of primitive variables
                        foreach (j; 0 .. nConserved) {
                            blk.celldata.U1[cellCount+j] = blk.celldata.U0[cellCount+j] + omega*blk.dU[cellCount+j];
                        }
                        try {
                            decode_conserved(blk.celldata.positions[i], blk.celldata.U1[cellCount .. cellCount+nConserved],
                                             blk.celldata.flowstates[i], blk.omegaz, i, blk.myConfig);
                        }
                        catch (FlowSolverException e) {
                            failed_decode = true;
                        }

                        // return cell to original state
                        blk.celldata.flowstates[i].copy_values_from(blk.fs_save);

                        if (failed_decode) {
                            omega *= omega_reduction_factor;
                            failed_decode = false;
                        } else {
                            // if we reach here we have a suitable relaxation factor for this cell
                            break;
                        }
                    }
                    cellCount += nConserved;
                }
            }

            // the relaxation factor may have been updated, so communicate minimum omega to all processes again
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &(omega), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }

        } // end physicality check

        // ensure the relaxation factor is sufficient enough to reduce the unsteady residual
        // this is done via a simple backtracking line search algorithm
        bool failed_line_search = false;
        if ( (omega > omega_min) && GlobalConfig.sssOptions.useLineSearch) {
            // residual at current state
            auto RU0 = normNew;
            // unsteady residual at updated state (we will fill this later)
            number RUn = 0.0;

            //  find omega such that the unsteady residual is reduced
            bool reduce_omega = true;
            while (reduce_omega) {

                // 1. compute unsteady term
                foreach (blk; parallel(localFluidBlocks,1)) {
                    int cellCount = 0;
                    foreach (cell; blk.cells) {
                        double dtInv;
                        if (blk.myConfig.with_local_time_stepping) { dtInv = 1.0/cell.dt_local; }
                        else { dtInv = 1.0/dt; }
                        foreach (j; 0 .. nConserved) {
                            blk.FU[cellCount+j] = -dtInv*omega*blk.dU[cellCount+j];
                        }
                        cellCount += nConserved;
                    }
                }

                // 2. compute residual at updated state term
                foreach (blk; parallel(localFluidBlocks,1)) {
                    int cellCount = 0;
                    foreach (cell; blk.cells) {
                        cell.U[1].copy_values_from(cell.U[0]);
                        foreach (j; 0 .. nConserved) {
                            cell.U[1][j] = cell.U[0][j] + omega*blk.dU[cellCount+j];
                        }
                        cell.decode_conserved(0, 1, blk.omegaz);
                        cellCount += nConserved;
                    }
                }
                evalRHS(0.0, 1);
                foreach (blk; parallel(localFluidBlocks,1)) {
                    int cellCount = 0;
                    foreach (cell; blk.cells) {
                        foreach (j; 0 .. nConserved) {
                            blk.FU[cellCount+j] += cell.dUdt[1][j].re;
                        }
                        // return cell to original state
                        cell.decode_conserved(0, 0, blk.omegaz);
                        cellCount += nConserved;
                    }
                }

                // 3. add smoothing source term
                if (GlobalConfig.residual_smoothing) {
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        int cellCount = 0;
                        foreach (cell; blk.cells) {
                            double dtInv;
                            if (blk.myConfig.with_local_time_stepping) { dtInv = 1.0/cell.dt_local; }
                            else { dtInv = 1.0/dt; }
                            foreach (k; 0 .. nConserved) {
                                blk.FU[cellCount+k] += dtInv*blk.DinvR[cellCount+k];
                            }
                            cellCount += nConserved;
                        }
                    }
                }

                // compute norm of unsteady residual
                mixin(dot_over_blocks("RUn", "FU", "FU"));
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &RUn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                }
                RUn = sqrt(RUn);

                // check if unsteady residual is reduced
                if (RUn < RU0 || omega < omega_min) {
                    reduce_omega = false;
                } else {
                    omega *= omega_reduction_factor;
                }
            }
        } // end line search

        if ( (omega < omega_min) && residual_based_cfl_scheduling)  {
            // the update isn't good, reduce the CFL and try again
            cfl = 0.5*cfl;
            dt = determine_dt(cfl);
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
            if ( GlobalConfig.is_master_task ) {
                writefln("WARNING: nonlinear update relaxation factor too small for step= %d", step);
                writefln("         Taking nonlinear step again with the CFL reduced by a factor of 0.5");
            }

            // we don't proceed with the nonlinear update for this step
            continue;
        }

        // If we get here, things are good
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (i; 0 .. blk.ncells*ncq) {
                blk.celldata.U1[i] = blk.celldata.U0[i] + omega*blk.dU[i];
            }
            foreach (cidx; 0 .. blk.ncells) {
                size_t idx = cidx*ncq;
                decode_conserved(blk.celldata.positions[cidx], blk.celldata.U1[idx .. idx+ncq], blk.celldata.flowstates[cidx], blk.omegaz, cidx, blk.myConfig);
            }
            // Put flow state into U[0] ready for next iteration.
            foreach (i; 0 .. blk.ncells*ncq) {
                blk.celldata.U0[i] = blk.celldata.U1[i];
            }
        }

        // after a successful fluid domain update, proceed to perform a solid domain update
        if (include_solid_domain && localSolidBlocks.length > 0) { solid_update(step, pseudoSimTime, cfl, eta, sigma); }

        pseudoSimTime += dt;
        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();

	if (!limiterFreezingCondition && (normNew/normRef <= limiterFreezingResidReduction)) {
	    countsBeforeFreezing++;
	    if (countsBeforeFreezing >= limiterFreezingCount) {
                if (GlobalConfig.frozen_limiter == false) {
                    evalRHS(pseudoSimTime, 0);
                    GlobalConfig.frozen_limiter = true;
                    GlobalConfig.frozen_shock_detector = true;
                }
                limiterFreezingCondition = true;
                writefln("=== limiter freezing condition met at step: %d ===", step);
            }
	}
        // Check on some stopping criteria
        if ( (omega < omega_min) && !residual_based_cfl_scheduling ) {
            if (GlobalConfig.is_master_task) {
                writefln("WARNING: The simulation is stopping because the nonlinear update relaxation factor is below the minimum allowable value.");
            }
            finalStep = true;
        }
        if ( (cfl < cfl_min) && residual_based_cfl_scheduling ) {
            if (GlobalConfig.is_master_task) {
                writefln("WARNING: The simulation is stopping because the CFL (%.3e) is below the minimum allowable CFL value (%.3e)", cfl, cfl_min);
            }
            finalStep = true;
        }
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

        if ((SimState.maxWallClockSeconds > 0) && (wallClockElapsed > SimState.maxWallClockSeconds)) {
            if (GlobalConfig.is_master_task) {
                writefln("Reached maximum wall-clock time with elapsed time %s.", to!string(wallClockElapsed));
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
            if ( fabs(mass_balance).re <= massBalanceReduction) {
                if (GlobalConfig.is_master_task) {
                    writeln("STOPPING: The global mass balance is below target value.");
                    writefln("          current value= %.12e   target value= %.12e", fabs(mass_balance).re, massBalanceReduction);
                }
                finalStep = true;
            }
            // Write out residuals
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            if (GlobalConfig.is_master_task) {
                auto cqi = GlobalConfig.cqi;
                fResid = File(residFname, "a");
                if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                    fResid.writef("%8d  %20.16e  %20.16e %20.16e %20.16e %3d %3d %5d %.8f %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e ",
                                  step, pseudoSimTime, dt, cfl, eta, nRestarts, nIters, fnCount, wallClockElapsed,
                                  normNew, normNew/normRef,
                                  currResiduals[cqi.mass].re, currResiduals[cqi.mass].re/maxResiduals[cqi.mass].re,
                                  currResiduals[cqi.xMom].re, currResiduals[cqi.xMom].re/maxResiduals[cqi.xMom].re,
                                  currResiduals[cqi.yMom].re, currResiduals[cqi.yMom].re/maxResiduals[cqi.yMom].re);
                } else {
                    fResid.writef("%8d  %20.16e  %20.16e %20.16e %20.16e %3d %3d %5d %.8f %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e ",
                                  step, pseudoSimTime, dt, cfl, eta, nRestarts, nIters, fnCount, wallClockElapsed,
                                  normNew, normNew/normRef,
                                  currResiduals[cqi.xMom].re, currResiduals[cqi.xMom].re/maxResiduals[cqi.xMom].re,
                                  currResiduals[cqi.yMom].re, currResiduals[cqi.yMom].re/maxResiduals[cqi.yMom].re);
                }
                if ( GlobalConfig.dimensions == 3 )
                    fResid.writef("%20.16e  %20.16e  ", currResiduals[cqi.zMom].re, currResiduals[cqi.zMom].re/maxResiduals[cqi.zMom].re);
                fResid.writef("%20.16e  %20.16e  ",
                              currResiduals[cqi.totEnergy].re, currResiduals[cqi.totEnergy].re/maxResiduals[cqi.totEnergy].re);
                foreach(it; 0 .. GlobalConfig.turb_model.nturb){
                    fResid.writef("%20.16e  %20.16e  ",
                                  currResiduals[cqi.rhoturb+it].re, currResiduals[cqi.rhoturb+it].re/maxResiduals[cqi.rhoturb+it].re);
                }
                version(multi_species_gas){
                if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                    foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                        fResid.writef("%20.16e  %20.16e  ",
                                      currResiduals[cqi.species+sp].re, currResiduals[cqi.species+sp].re/maxResiduals[cqi.species+sp].re);
                    }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
                    fResid.writef("%20.16e  %20.16e  ",
                                  currResiduals[cqi.modes+imode].re, currResiduals[cqi.modes+imode].re/maxResiduals[cqi.modes+imode].re);
                }
                }
                fResid.writef("%20.16e ", fabs(mass_balance.re));
                fResid.writef("%20.16e ", linSolResid);
                if (GlobalConfig.sssOptions.useLineSearch || GlobalConfig.sssOptions.usePhysicalityCheck) {
                    fResid.writef("%20.16e ", omega);
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
                if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                    formattedWrite(writer, "  mass           %10.6e    %10.6e\n", currResiduals[cqi.mass].re, currResiduals[cqi.mass].re/maxResiduals[cqi.mass].re);
                }
                formattedWrite(writer, "  x-mom          %10.6e    %10.6e\n", currResiduals[cqi.xMom].re, currResiduals[cqi.xMom].re/maxResiduals[cqi.xMom].re);
                formattedWrite(writer, "  y-mom          %10.6e    %10.6e\n", currResiduals[cqi.yMom].re, currResiduals[cqi.yMom].re/maxResiduals[cqi.yMom].re);
                if ( GlobalConfig.dimensions == 3 )
                    formattedWrite(writer, "  z-mom          %10.6e    %10.6e\n", currResiduals[cqi.zMom].re, currResiduals[cqi.zMom].re/maxResiduals[cqi.zMom].re);
                formattedWrite(writer, "  total-energy   %10.6e    %10.6e\n", currResiduals[cqi.totEnergy].re, currResiduals[cqi.totEnergy].re/maxResiduals[cqi.totEnergy].re);
                foreach(it; 0 .. GlobalConfig.turb_model.nturb){
                    auto tvname = GlobalConfig.turb_model.primitive_variable_name(it);
                    formattedWrite(writer, "  %s            %10.6e    %10.6e\n", tvname, currResiduals[cqi.rhoturb+it].re, currResiduals[cqi.rhoturb+it].re/maxResiduals[cqi.rhoturb+it].re);
                }
                version(multi_species_gas){
                if ( GlobalConfig.gmodel_master.n_species > 1 ) {
                    foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                        auto spname = GlobalConfig.gmodel_master.species_name(sp);
                        formattedWrite(writer, "  %s            %10.6e    %10.6e\n", spname, currResiduals[cqi.species+sp].re, currResiduals[cqi.species+sp].re/maxResiduals[cqi.species+sp].re);
                    }
                }
                }
                version(multi_T_gas){
                foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
                    auto modename = "T_MODES["~to!string(imode)~"]"; //GlobalConfig.gmodel_master.energy_mode_name(imode);
                    formattedWrite(writer, "  %s            %10.6e    %10.6e\n", modename, currResiduals[cqi.modes+imode].re, currResiduals[cqi.modes+imode].re/maxResiduals[cqi.modes+imode].re);
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
                    } else {
                        foreach(io; io_list) {
                            string path = "CellData/"~io.tag;
                            if (io.do_save()) ensure_directory_is_present(make_path_name(path, nWrittenSnapshots));
                        }
                    }
                    ensure_directory_is_present(make_path_name!"solid"(nWrittenSnapshots));
                }
                version(mpi_parallel) {
                    MPI_Barrier(MPI_COMM_WORLD);
                }
                foreach (blk; localFluidBlocks) {
                    if (legacy) {
                        auto fileName = make_file_name!"flow"(jobName, blk.id, nWrittenSnapshots, GlobalConfig.flowFileExt);
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
            // Writing to file produces a large amount of temporary storage that needs to be cleaned up.
            // Forcing the Garbage Collector to go off here prevents this freeable memory from
            // accumulating over time, which can upset the queueing systems in HPC Jobs.
            // 24/05/22 (NNG)
            GC.collect();
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
                if (omega >= omega_allow_cfl_growth) {
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
    /*
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
    */
}

void allocate_global_fluid_workspace()
{
    size_t mOuter = to!size_t(GlobalConfig.sssOptions.maxOuterIterations);
    g0.length = mOuter+1;
    g1.length = mOuter+1;
    h.length = mOuter+1;
    hR.length = mOuter+1;
    H0 = new Matrix!double(mOuter+1, mOuter);
    H1 = new Matrix!double(mOuter+1, mOuter);
    Gamma = new Matrix!double(mOuter+1, mOuter+1);
    Q0 = new Matrix!double(mOuter+1, mOuter+1);
    Q1 = new Matrix!double(mOuter+1, mOuter+1);
}

double determine_dt(double cflInit)
{
    double signal;
    double dt = double.max;
    foreach (blk; localFluidBlocks) {
        if (blk.myConfig.sssOptions.inviscidCFL) {
            foreach (i; 0 .. blk.ncells) {
                // calculate the signal using the maximum inviscid wave speed
                // ref. Computational Fluid Dynamics: Principles and Applications, J. Blazek, 2015, pg. 175
                // Note: this approximates cell width dx as the cell volume divided by face area
                signal = 0.0;
                foreach (fidx; blk.celldata.c2f[i]) {
                    Vector3 vel = blk.facedata.flowstates[fidx].vel;
                    Vector3 n   = blk.facedata.normals[fidx];
                    number a = blk.facedata.flowstates[fidx].gas.a;
                    number area = blk.facedata.areas[fidx];

                    number un = fabs(vel.dot(n));
                    number signal_f = (un+a)*area;
                    signal += signal_f.re;
                }
                signal *= (1.0/blk.celldata.volumes[i].re);

                double dt_local = cflInit/signal;
                blk.celldata.dt_local[i] = dt_local;
                dt = fmin(dt, dt_local);
            }
        } else {
            // use the default signal frequency routine from the time-accurate code path
            foreach (i, cell; blk.cells) {
                signal = cell.signal_frequency();

                double dt_local = cflInit/signal;
                cell.dt_local = dt_local;
                blk.celldata.dt_local[i] = dt_local;
                dt = fmin(dt, dt_local);
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
        blk.clear_cell_source_vectors();
    }
    exchange_ghost_cell_boundary_data(pseudoSimTime, 0, ftl);
    foreach (blk; localFluidBlocks) {
        blk.applyPreReconAction(pseudoSimTime, 0, ftl);
    }

    // We don't want to switch between flux calculator application while
    // doing the Frechet derivative, so we'll only search for shock points
    // at ftl = 0, which is when the F(U) evaluation is made.
    if (ftl == 0 && GlobalConfig.do_shock_detect && GlobalConfig.frozen_shock_detector == false) { detect_shocks(0, ftl); }

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
        blk.convective_flux_phase0new(allow_high_order_interpolation);
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(pseudoSimTime, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase1new(allow_high_order_interpolation);
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(pseudoSimTime, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase2new(allow_high_order_interpolation);
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
                blk.average_lsq_cell_derivs_to_faces();
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
            blk.applyPostDiffFluxAction(pseudoSimTime, 0, ftl);
        }
    }

    immutable int step = SimState.step;
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.eval_fluid_source_vectors(blk.omegaz);
        blk.eval_thermochem_source_vector(step);
        blk.eval_udf_source_vectors(pseudoSimTime);
        blk.time_derivatives(0, ftl);
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
    // remove the conserved mass variable for multi-species gas
    //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
    immutable size_t ncq  = nConserved; // number of conserved quantities

    // We perform a Frechet derivative to evaluate J*D^(-1)v
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.clear_fluxes_of_conserved_quantities();
        blk.clear_cell_source_vectors();

        foreach(i; 0 .. ncq*blk.ncells){
            blk.celldata.U1[i] = blk.celldata.U0[i] + sigma*blk.zed[i];
        }
        foreach (i; 0 .. blk.ncells) {
            size_t idx = i*ncq;
            decode_conserved(blk.celldata.positions[i], blk.celldata.U1[idx .. idx+ncq], blk.celldata.flowstates[i], blk.omegaz, i, blk.myConfig);
        }
    }
    evalRHS(pseudoSimTime, 1);
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach(i; 0 .. blk.ncells*ncq){
            blk.zed[i] = (blk.celldata.dUdt1[i].re - blk.FU[i])/sigma;
        }
        foreach (i; 0 .. blk.ncells) {
            size_t idx = i*ncq;
            decode_conserved(blk.celldata.positions[i], blk.celldata.U0[idx .. idx+ncq], blk.celldata.flowstates[i], blk.omegaz, i, blk.myConfig);
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
        // remove the conserved mass variable for multi-species gas
        //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
        immutable size_t ncq  = nConserved; // number of conserved quantities

        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.clear_fluxes_of_conserved_quantities();
            blk.clear_cell_source_vectors();
            
            foreach(i; 0 .. blk.ncells*ncq){
                blk.celldata.U1[i] = blk.celldata.U0[i];
                blk.celldata.U1[i] += complex(0.0, sigma*blk.zed[i].re);
            }

            foreach (i; 0 .. blk.ncells) {
                size_t s1 = i*ncq;
                decode_conserved(blk.celldata.positions[i], blk.celldata.U1[s1 .. s1+ncq], blk.celldata.flowstates[i], blk.omegaz, i, blk.myConfig);
            }
        }
        evalRHS(pseudoSimTime, 1);
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach(i; 0 .. blk.ncells*ncq){
                blk.zed[i] = blk.celldata.dUdt1[i].im/(sigma);
            }

            foreach(i; 0 .. blk.ncells) { blk.celldata.flowstates[i].clear_imaginary_components(); }
            foreach(bc; blk.bc) {
                foreach(ghostcell; bc.ghostcells) { ghostcell.fs.clear_imaginary_components(); }
            }
            foreach(i; 0 .. blk.nfaces) { blk.facedata.flowstates[i].clear_imaginary_components(); }
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

void point_implicit_relaxation_solve(int step, double pseudoSimTime, double dt, ref double residual, ref double linSolResid, int startStep)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    // remove the conserved mass variable for multi-species gas
    //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }

    // we start with a guess of dU = 0
    foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }

    // evaluate R.H.S. residual vector R
    evalRHS(pseudoSimTime, 0);

    // calculate global residual
    foreach (blk; parallel(localFluidBlocks,1)) {
        int cellCount = 0;
        foreach (cell; blk.cells) {
            foreach (j; 0 .. nConserved) { blk.FU[cellCount+j] = cell.dUdt[0][j].re; }
            cellCount += nConserved;
        }
    }
    mixin(dot_over_blocks("residual", "FU", "FU"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(residual), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    residual = sqrt(residual);

    // evaluate Jacobian for use in the implicit operator [A] for matrix-based solvers
    if (step == startStep || step%GlobalConfig.sssOptions.frozenPreconditionerCount == 0) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.evaluate_jacobian();
            blk.flowJacobian.augment_with_dt(blk.celldata.dt_local, dt, blk.ncells, nConserved);
            nm.smla.invert_block_diagonal(blk.flowJacobian.local, blk.flowJacobian.D, blk.flowJacobian.Dinv, blk.cells.length, nConserved);
        }
    }

    // solve linear system [A].dU = R
    final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
    case PreconditionMatrixType.diagonal:
        mixin(diagonal_solve("dU", "FU"));
        break;
    case PreconditionMatrixType.jacobi:
        mixin(jacobi_solve("dU", "FU"));
        break;
    case PreconditionMatrixType.sgs:
        mixin(sgs_solve("dU", "FU"));
        break;
    case PreconditionMatrixType.ilu:
        throw new Error("ILU cannot be used in the point implicit relaxation update.");
    } // end switch

    // compute linear solve residual
    foreach (blk; parallel(localFluidBlocks,1)) {
        // we need to temporarily reinvert the diagonal blocks, since we store A = L + D^(-1) + U
        nm.smla.invert_block_diagonal(blk.flowJacobian.local, blk.flowJacobian.D, blk.flowJacobian.Dinv, blk.cells.length, nConserved);
        nm.smla.multiply(blk.flowJacobian.local, blk.dU, blk.rhs[]);
        nm.smla.invert_block_diagonal(blk.flowJacobian.local, blk.flowJacobian.D, blk.flowJacobian.Dinv, blk.cells.length, nConserved);
        blk.rhs[] = blk.rhs[] - blk.FU[];
    }
    mixin(dot_over_blocks("linSolResid", "rhs", "rhs"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(linSolResid), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    linSolResid = sqrt(linSolResid)/residual; // relative drop in the residual
}

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

    int kmax = GlobalConfig.sssOptions.maxSubIterations;
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

    int kmax = GlobalConfig.sssOptions.maxSubIterations;
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

void rpcGMRES_solve(int step, double pseudoSimTime, double dt, double eta, double sigma, bool usePreconditioner,
                    ref double residual, ref int nRestarts, ref int nIters, ref double linSolResid, int startStep, int LHSeval, int RHSeval, ref bool pc_matrix_evaluated)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    // remove the conserved mass variable for multi-species gas
    //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
    size_t MASS = GlobalConfig.cqi.mass;
    size_t X_MOM = GlobalConfig.cqi.xMom;
    size_t Y_MOM = GlobalConfig.cqi.yMom;
    size_t Z_MOM = GlobalConfig.cqi.zMom;
    size_t TOT_ENERGY = GlobalConfig.cqi.totEnergy;
    size_t TKE = GlobalConfig.cqi.rhoturb;
    size_t SPECIES = GlobalConfig.cqi.species;
    size_t MODES = GlobalConfig.cqi.modes;
    immutable size_t ncq = nConserved;

    int interpOrderSave = GlobalConfig.interpolation_order;
    // We add 1 because the user thinks of "re"starts, so they
    // might legitimately ask for no restarts. We still have
    // to execute at least once.
    int maxRestarts = GlobalConfig.sssOptions.maxRestarts + 1;
    int maxIters = GlobalConfig.sssOptions.maxOuterIterations;
    size_t m = to!size_t(maxIters);
    size_t r;
    size_t iterCount;
    double beta, beta0;
    double outerTol;
    number resid;

    // Variables for max rates of change
    // Use these for equation scaling.
    double minNonDimVal = 1.0; // minimum value used for non-dimensionalisation
                               // when our time rates of change are very small
                               // then we'll avoid non-dimensionalising by
                               // values close to zero.

    // Compute the RHS residual and store dUdt[0] as F(U)
    evalRHS(pseudoSimTime, 0);
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach(i; 0 .. blk.ncells*ncq) {
            blk.FU[i] = blk.celldata.dUdt0[i].re;
        }
    }

    // Compute the approximate Jacobian matrix for preconditioning, if requested
    pc_matrix_evaluated = false;
    if (usePreconditioner && ( (m == nIters && GlobalConfig.sssOptions.useAdaptivePreconditioner) ||
                               (step == startStep) ||
                               (step%GlobalConfig.sssOptions.frozenPreconditionerCount == 0) )) {
        pc_matrix_evaluated = true;
        final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
        case PreconditionMatrixType.diagonal:
            goto case PreconditionMatrixType.sgs;
        case PreconditionMatrixType.jacobi:
            goto case PreconditionMatrixType.sgs;
        case PreconditionMatrixType.sgs:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.evaluate_jacobian();
                blk.flowJacobian.augment_with_dt(blk.celldata.dt_local, dt, blk.ncells, nConserved);
                nm.smla.invert_block_diagonal(blk.flowJacobian.local, blk.flowJacobian.D, blk.flowJacobian.Dinv, blk.cells.length, nConserved);
            }
            break;
        case PreconditionMatrixType.ilu:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.evaluate_jacobian();
                blk.flowJacobian.augment_with_dt(blk.celldata.dt_local, dt, blk.ncells, nConserved);
                nm.smla.decompILU0(blk.flowJacobian.local);
            }
            break;
        }
    }

    // Apply residual smoothing to RHS, if requested
    // ref. A Residual Smoothing Strategy for Accelerating Newton Method Continuation, D. J. Mavriplis, Computers & Fluids, 2021
    if (GlobalConfig.residual_smoothing) {
        // compute approximate solution via dU = D^{-1}*F(U) where we set D = precondition matrix
        final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
        case PreconditionMatrixType.diagonal:
            foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = 0.0; }
            mixin(diagonal_solve("DinvR", "FU"));
            break;
        case PreconditionMatrixType.jacobi:
            foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = 0.0; }
            mixin(jacobi_solve("DinvR", "FU"));
            break;
        case PreconditionMatrixType.sgs:
            foreach (blk; parallel(localFluidBlocks,1)) { blk.DinvR[] = 0.0; }
            mixin(sgs_solve("DinvR", "FU"));
            break;
        case PreconditionMatrixType.ilu:
            foreach (blk; parallel(localFluidBlocks,1)) {
                blk.DinvR[] = blk.FU[];
                nm.smla.solve(blk.flowJacobian.local, blk.DinvR);
            }
            break;
        } // end switch

        // add smoothing source term to RHS
        foreach (blk; parallel(localFluidBlocks,1)) {
            int cellCount = 0;
            foreach (cell; blk.cells) {
                double dtInv;
                if (blk.myConfig.with_local_time_stepping) { dtInv = 1.0/cell.dt_local; }
                else { dtInv = 1.0/dt; }
                foreach (k; 0 .. nConserved) {
                    blk.FU[cellCount+k] += dtInv*blk.DinvR[cellCount+k];
                }
                cellCount += nConserved;
            }
        }
    }

    // Determine the max rates in F(U) for scaling the linear system
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.maxRate.clear();
        // determine max rates
        foreach (i; 0 .. blk.ncells) {
            foreach(j; 0 .. ncq){
                size_t idx = i*ncq;
                blk.maxRate[j] = fmax(blk.maxRate[j], blk.FU[idx+j]);
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
        if ( blk.myConfig.gmodel.n_species == 1 ) { maxMass = fmax(maxMass, blk.maxRate[cqi.mass]); }
        maxMomX = fmax(maxMomX, blk.maxRate[cqi.xMom]);
        maxMomY = fmax(maxMomY, blk.maxRate[cqi.yMom]);
        if ( blk.myConfig.dimensions == 3 )
            maxMomZ = fmax(maxMomZ, blk.maxRate[cqi.zMom]);
        maxEnergy = fmax(maxEnergy, blk.maxRate[cqi.totEnergy]);
        foreach(it; 0 .. blk.myConfig.turb_model.nturb){
            maxTurb[it] = fmax(maxTurb[it], blk.maxRate[cqi.rhoturb+it]);
        }
        version(multi_species_gas){
        if ( blk.myConfig.gmodel.n_species > 1 ) {
            foreach(sp; 0 .. blk.myConfig.gmodel.n_species){
                maxSpecies[sp] = fmax(maxSpecies[sp], blk.maxRate[cqi.species+sp]);
            }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. blk.myConfig.gmodel.n_modes){
            maxModes[imode] = fmax(maxModes[imode], blk.maxRate[cqi.modes+imode]);
        }
        }
    }
    // In distributed memory, reduce the max values and ensure everyone has a copy
    version(mpi_parallel) {
        if ( GlobalConfig.gmodel_master.n_species == 1 ) { MPI_Allreduce(MPI_IN_PLACE, &(maxMass.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); }
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
    if ( GlobalConfig.gmodel_master.n_species == 1 ) { maxMass = fmax(maxMass, minNonDimVal); }
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
            if ( blk.myConfig.gmodel.n_species == 1 ) { blk.maxRate[cqi.mass] = maxMass; }
            blk.maxRate[cqi.xMom] = maxMomX;
            blk.maxRate[cqi.yMom] = maxMomY;
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate[cqi.zMom] = maxMomZ;
            blk.maxRate[cqi.totEnergy] = maxEnergy;
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                blk.maxRate[cqi.rhoturb+it] = maxTurb[it];
            }
            version(multi_species_gas){
            if ( blk.myConfig.gmodel.n_species > 1 ) {
                foreach(sp; 0 .. blk.myConfig.gmodel.n_species){
                    blk.maxRate[cqi.species+sp] = maxSpecies[sp];
                }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. blk.myConfig.gmodel.n_modes){
                blk.maxRate[cqi.modes+imode] = maxModes[imode];
            }
            }
        }
        else { // just scale by 1
            if ( blk.myConfig.gmodel.n_species == 1 ) { blk.maxRate[cqi.mass] = 1.0; }
            blk.maxRate[cqi.xMom] = 1.0;
            blk.maxRate[cqi.yMom] = 1.0;
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate[cqi.zMom] = 1.0;
            blk.maxRate[cqi.totEnergy] = 1.0;
            foreach(it; 0 .. blk.myConfig.turb_model.nturb){
                blk.maxRate[cqi.rhoturb+it] = 1.0;
            }
            version(multi_species_gas){
            if ( blk.myConfig.gmodel.n_species ) {
                foreach(sp; 0 .. blk.myConfig.gmodel.n_species){
                    blk.maxRate[cqi.species+sp] = 1.0;
                }
            }
            }
            version(multi_T_gas){
            foreach(imode; 0 .. blk.myConfig.gmodel.n_modes){
                blk.maxRate[cqi.modes+imode] = 1.0;
            }
            }

        }
    }

    // Compute the unscaled L2 norm for reporting the residual of the non-linear system of equations
    double unscaledNorm2;
    // we optionally remove the turbulent conserved quantities from the global residual since their
    // magnitude can be several orders larger than the flow quantities, and thus their convergence
    // history can adversely dominate the CFL growth algorithm.
    if (GlobalConfig.sssOptions.include_turb_quantities_in_residual == false) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            auto cqi = blk.myConfig.cqi;
            int cellCount = 0;
            foreach (i; 0 .. blk.ncells) {
                foreach(it; 0 .. cqi.n_turb) { blk.FU[cellCount+TKE+it] = 0.0; }
                cellCount += nConserved;
            }
        }
    }

    // L2 norm calculation
    mixin(dot_over_blocks("unscaledNorm2", "FU", "FU"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &unscaledNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    unscaledNorm2 = sqrt(unscaledNorm2);

    // if we removed the turbulent conserved quantities in the global residual, we have to add them back in
    if (GlobalConfig.sssOptions.include_turb_quantities_in_residual == false) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            auto cqi = blk.myConfig.cqi;
            int cellCount = 0;
            foreach (i; 0 .. blk.ncells) {
                foreach(it; 0 .. cqi.n_turb) { blk.FU[cellCount+TKE+it] = blk.celldata.dUdt0[cellCount+cqi.rhoturb+it].re; }
                cellCount += nConserved;
            }
        }
    }

    // the remainder of the routine closely follows the structure of rpcGMRES found in nm/smla.d
    // we set the initial guess to zero
    foreach (blk; parallel(localFluidBlocks,1)) { blk.x0[] = 0.0; }

    // Start outer-loop of restarted GMRES
    for ( r = 0; r < maxRestarts; r++ ) {

        // 0. Initialise some arrays and matrices that have already been allocated
        g0[] = 0.0;
        g1[] = 0.0;
        H0.zeros();
        H1.zeros();
        Gamma.eye();

        // 1. Evaluate r0 = b - A.x0, beta, v1

        // evaluate A.x0 using a Frechet derivative (note that the zed[] array is baked into the evalJacobianVecProd routine).
        foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = blk.x0[]; }

        // Prepare 'w' with (I/dt)(P^-1)v term;
        if (GlobalConfig.with_local_time_stepping) {
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (i; 0 .. blk.ncells) {
                    double dtInv = 1.0/blk.celldata.dt_local[i];
                    size_t offset = i*nConserved;
                    foreach (k; 0 .. nConserved) {
                        blk.w[offset+k] = dtInv*blk.zed[offset+k];
                    }
                }
            }
        } else { // Global timestepping version:
            double dtInv = 1.0/dt;
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (idx; 0 .. blk.nvars) {
                    blk.w[idx] = dtInv*blk.zed[idx];
                }
            }
        }

        // Evaluate Jz and place in z
        evalJacobianVecProd(pseudoSimTime, sigma, LHSeval, RHSeval);

        // Now we can complete calculation of r0
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars)  blk.r0[k] = blk.FU[k] - (blk.w[k] - blk.zed[k]);
        }

        // apply the system scaling to r0
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach(i; 0 .. blk.ncells){
                size_t idx = i*ncq;
                foreach(j; 0 .. ncq){
                    blk.r0[idx+j] *= 1.0/blk.maxRate[j].re;
                }
            }
        }
        // then compute v = r0/||r0|| and set first residual entry
        mixin(dot_over_blocks("beta", "r0", "r0"));
        version(mpi_parallel) {
            // NOTE: this dot product has been observed to be sensitive to the order of operations,
            //       the use of MPI_Allreduce means that the same convergence behaviour can not be expected
            //       for a different mapping of blocks over the MPI tasks and/or a shared memory calculation
            //       2022-09-30 (KAD).
            MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        beta = sqrt(beta);
        g0[0] = beta;
        //writefln("beta= %e", beta);
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.r0[k]/beta;
                blk.VT[k] = blk.v[k];
            }
        }

        // Compute outer tolerance on first restart and store initial residual
        if (r == 0) {
            outerTol = eta*beta;
            beta0 = beta;
        }

        // 2. Do 'm' iterations of update
        foreach (j; 0 .. m) {
            iterCount = j+1;

            // Undo the linear system scaling for Jacobian-vector evaluation
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach(i; 0 .. blk.ncells){
                    size_t idx = i*ncq;
                    foreach(j; 0 .. ncq){
                        blk.v[idx+j] *= blk.maxRate[j].re;
                    }
                }
            }

            // Apply preconditioning step
            if (usePreconditioner && step >= GlobalConfig.sssOptions.startPreconditioning) {
                final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
                case PreconditionMatrixType.diagonal:
                    foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = 0.0; }
                    mixin(diagonal_solve("zed", "v"));
                    break;
                case PreconditionMatrixType.jacobi:
                    foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = 0.0; }
                    mixin(jacobi_solve("zed", "v"));
                    break;
                case PreconditionMatrixType.sgs:
                    foreach (blk; parallel(localFluidBlocks,1)) { blk.zed[] = 0.0; }
                    mixin(sgs_solve("zed", "v"));
                    break;
                case PreconditionMatrixType.ilu:
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        blk.zed[] = blk.v[];
                        nm.smla.solve(blk.flowJacobian.local, blk.flowJacobian.local_diags, blk.zed);
                    }
                    break;
                } // end switch
            } else {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    blk.zed[] = blk.v[];
                }
            }

            // evaluate Jz using either a real or complex valued Frechet derivative...
            // for the real-valued solver we compute an ideal perturbation parameter for the Frechet derivative as per equation 11/12 from
            // Knoll, D.A. and McHugh, P.R., Newton-Krylov Methods Applied to a System of Convection-Diffusion-Reaction Equations, 1994
            // note: calculate this parameter WITHOUT scaling applied
            if (!GlobalConfig.sssOptions.useComplexMatVecEval) {
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
                    foreach(i; 0 .. blk.nvars) {
                        sume += eps0*abs(blk.celldata.U0[i]) + eps0;
                    }
                    N += blk.nvars;
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
            }

            // Prepare 'w' with (I/dt)(P^-1)v term;
            if (GlobalConfig.with_local_time_stepping) {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (i; 0 .. blk.ncells) {
                        double dtInv = 1.0/blk.celldata.dt_local[i];
                        size_t offset = i*nConserved;
                        foreach (k; 0 .. nConserved) {
                            blk.w[offset+k] = dtInv*blk.zed[offset+k];
                        }
                    }
                }
            } else { // Global timestepping version:
                double dtInv = 1.0/dt;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (idx; 0 .. blk.nvars) {
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

            // apply the linear system scaling
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach(i; 0 .. blk.ncells){
                    size_t idx = i*ncq;
                    foreach(j; 0 .. ncq){
                        blk.w[idx+j] *= 1.0/blk.maxRate[j].re;
                    }
                }
            }

            // The remainder of the algorithm looks a lot like any standard GMRES implementation (for example, see smla.d)
            foreach (i; 0 .. j+1) {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    // Extract column 'i'
                    size_t offset = i*blk.nvars;
                    foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.VT[offset + k];
                }
                double H0_ij;
                mixin(dot_over_blocks("H0_ij", "w", "v"));
                version(mpi_parallel) {
                    // NOTE: this dot product has been observed to be sensitive to the order of operations,
                    //       the use of MPI_Allreduce means that the same convergence behaviour can not be expected
                    //       for a different mapping of blocks over the MPI tasks and/or a shared memory calculation
                    //       2022-09-30 (KAD).
                    MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                }
                H0[i,j] = H0_ij;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k];
                }
            }
            double H0_jp1j;
            mixin(dot_over_blocks("H0_jp1j", "w", "w"));
            version(mpi_parallel) {
                // NOTE: this dot product has been observed to be sensitive to the order of operations,
                //       the use of MPI_Allreduce means that the same convergence behaviour can not be expected
                //       for a different mapping of blocks over the MPI tasks and/or a shared memory calculation
                //       2022-09-30 (KAD).
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
                nm.bbla.dot!double(Gamma, j+2, j+2, Q0, j+2, Q1);
            }

            // Get residual
            resid = fabs(g1[j+1]);
            //writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
            nIters = to!int(iterCount);
            linSolResid = (resid/beta0).re;
            if ( resid <= outerTol ) {
                m = j+1;
                // DEBUG:
                //      writefln("OUTER: TOL ACHIEVED restart-count= %d iteration-count= %d, resid= %e", r, m, resid);
                //      writefln("RANK %d: tolerance achieved on iteration: %d", GlobalConfig.mpi_rank_for_local_task, m);
                break;
            }

            // Prepare for next iteration
            copy(H1, H0);
            g0[] = g1[];
            copy(Q1, Q0);
        }

        if (iterCount == maxIters)
            m = maxIters;

        // At end H := R up to row m
        //        g := gm up to row m
        upperSolve!double(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.transpose_and_dot!double(blk.VT, blk.nvars, m, blk.nvars, blk.g1, blk.zed);
        }

        // Undo the linear system scaling to recover the unscaled solution vector
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach(i; 0 .. blk.ncells){
                size_t idx = i*ncq;
                foreach(j; 0 .. ncq){
                    blk.zed[idx+j] *= blk.maxRate[j].re;
                }
            }
        }

        // Apply preconditioning step
        if (usePreconditioner && step >= GlobalConfig.sssOptions.startPreconditioning) {
            final switch (GlobalConfig.sssOptions.preconditionMatrixType) {
                case PreconditionMatrixType.diagonal:
                    foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }
                    mixin(diagonal_solve("dU", "zed"));
                    break;
                case PreconditionMatrixType.jacobi:
                    foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }
                    mixin(jacobi_solve("dU", "zed"));
                    break;
                case PreconditionMatrixType.sgs:
                    foreach (blk; parallel(localFluidBlocks,1)) { blk.dU[] = 0.0; }
                    mixin(sgs_solve("dU", "zed"));
                    break;
                case PreconditionMatrixType.ilu:
                    foreach(blk; parallel(localFluidBlocks,1)) {
                        blk.dU[] = blk.zed[];
                        nm.smla.solve(blk.flowJacobian.local, blk.dU);
                    }
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

        // Else, prepare for a restart by setting the inital
        // guess to the current best estimate of the solution
        foreach (blk; parallel(localFluidBlocks,1)) { blk.x0[] = blk.dU[]; }

        /*
          We were originally peforming the computation of r0 as per Fraysse etal (2005),
          however, we have observed that the resid value at the end of the last GMRES iteration
          and the beta value computed using the r0 as per Fraysse differed quite substantially,
          sometimes by an order of magnitude. It hasn't been determined whether this is due to
          a bug in the implementation of the trick or a fundamental deficiency of the method.
          We have thus reverted to computing the r0 vector explicitly, i.e. r0 = b - Ax0 where Ax0
          is computed using a Frechet derivative, in doing so the resid and beta values have a far
          better agreement. [KAD 2023-07-18]

          // old Fraysse method...
          foreach (blk; localFluidBlocks) copy(Q1, blk.Q1);
          // Set all values in g0 to 0.0 except for final (m+1) value
          foreach (i; 0 .. m) g0[i] = 0.0;
          foreach (blk; localFluidBlocks) blk.g0[] = g0[];
          foreach (blk; parallel(localFluidBlocks,1)) {
          nm.bbla.dot(blk.Q1, m, m+1, blk.g0, blk.g1);
          }
          foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.transpose_and_dot(blk.VT, blk.nvars, m+1, blk.nvars, blk.g1, blk.r0);
            }
        */
    }

    residual = unscaledNorm2;
    nRestarts = to!int(r);
}

void max_residuals(ref ConservedQuantities residuals)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = GlobalConfig.cqi.n;
    // remove the conserved mass variable for multi-species gas
    //if (GlobalConfig.cqi.n_species > 1) { nConserved -= 1; }
    immutable size_t ncq = nConserved;

    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach(j; 0 .. ncq) {
            blk.residuals[j] = fabs(blk.celldata.dUdt0[0+j]);
        }

        foreach(i; 0 .. blk.ncells){
            size_t idx = i*ncq;
            foreach(j; 0 .. ncq){
                blk.residuals[j] = fmax(blk.residuals[j], fabs(blk.celldata.dUdt0[idx+j]));
            }
        }
    }
    residuals.copy_values_from(localFluidBlocks[0].residuals);
    foreach (blk; localFluidBlocks) {
        foreach(j; 0 .. ncq){
            residuals[j] = fmax(residuals[j], blk.residuals[j]);
        }
    }
    version(mpi_parallel) {
        // In MPI context, only the master task (rank 0) has collated the residuals
        double maxResid;
        auto cqi = GlobalConfig.cqi;
        if ( GlobalConfig.gmodel_master.n_species == 1 ) {
            MPI_Reduce(&(residuals[cqi.mass].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals[cqi.mass].re = maxResid; }
        }
        MPI_Reduce(&(residuals[cqi.xMom].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals[cqi.xMom] = to!number(maxResid); }
        MPI_Reduce(&(residuals[cqi.yMom].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals[cqi.yMom] = to!number(maxResid); }
        if (GlobalConfig.dimensions == 3) {
            MPI_Reduce(&(residuals[cqi.zMom].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals[cqi.zMom] = to!number(maxResid); }
        }
        MPI_Reduce(&(residuals[cqi.totEnergy].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) { residuals[cqi.totEnergy].re = maxResid; }
        foreach (it; 0 .. GlobalConfig.turb_model.nturb) {
            MPI_Reduce(&(residuals[cqi.rhoturb+it].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals[cqi.rhoturb+it].re = maxResid; }
        }
        version(multi_species_gas){
        if ( GlobalConfig.gmodel_master.n_species > 1 ) {
            foreach (sp; 0 .. GlobalConfig.gmodel_master.n_species) {
                MPI_Reduce(&(residuals[cqi.species+sp].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                if (GlobalConfig.is_master_task) { residuals[cqi.species+sp].re = maxResid; }
            }
        }
        }
        version(multi_T_gas){
        foreach (imode; 0 .. GlobalConfig.gmodel_master.n_modes) {
            MPI_Reduce(&(residuals[cqi.modes+imode].re), &maxResid, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) { residuals[cqi.modes+imode].re = maxResid; }
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
            if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e",
                         i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                         rInfo.globalResidual.re, rInfo.residuals[cqi.mass].re,
                         rInfo.residuals[cqi.xMom].re, rInfo.residuals[cqi.yMom].re,
                         rInfo.residuals[cqi.totEnergy].re);
            } else {
                f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e",
                         i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                         rInfo.globalResidual.re,
                         rInfo.residuals[cqi.xMom].re, rInfo.residuals[cqi.yMom].re,
                         rInfo.residuals[cqi.totEnergy].re);
            }
        }
        else {
            if ( GlobalConfig.gmodel_master.n_species == 1 ) {
                f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e %.18e",
                         i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                         rInfo.globalResidual.re, rInfo.residuals[cqi.mass].re,
                         rInfo.residuals[cqi.xMom].re, rInfo.residuals[cqi.yMom].re, rInfo.residuals[cqi.zMom].re,
                         rInfo.residuals[cqi.totEnergy].re);
            } else {
                f.writef("%04d %.18e %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e",
                         i, rInfo.pseudoSimTime.re, rInfo.dt.re, rInfo.cfl.re, rInfo.step,
                         rInfo.globalResidual.re,
                         rInfo.residuals[cqi.xMom].re, rInfo.residuals[cqi.yMom].re, rInfo.residuals[cqi.zMom].re,
                         rInfo.residuals[cqi.totEnergy].re);
            }
        }
        foreach(it; 0 .. GlobalConfig.turb_model.nturb){
            f.writef(" %.18e", rInfo.residuals[cqi.rhoturb+it].re);
        }
        version(multi_species_gas){
        if ( GlobalConfig.gmodel_master.n_species > 1 ) {
            foreach(sp; 0 .. GlobalConfig.gmodel_master.n_species){
                f.writef(" %.18e", rInfo.residuals[cqi.species+sp].re);
            }
        }
        }
        version(multi_T_gas){
        foreach(imode; 0 .. GlobalConfig.gmodel_master.n_modes){
            f.writef(" %.18e", rInfo.residuals[cqi.modes+imode].re);
        }
        }
        f.write("\n");
    }
    f.close();
}
