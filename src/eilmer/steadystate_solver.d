/** steadystate_solver.d
 * Newton-Krylov updates for steady-state convergence.
 *
 * Author: Rowan G.
 * Date: 2016-10-09
 *
 * Note: This is the first attempt at 'production' code.
 * Some test implementations began on 2016-03-29.
 *
 * History:
 *   2016-10-17 : Implement right-preconditioned GMRES and
 *                add a preconditioner.
 *   2016-10-28 : Add right-preconditioned flexible GMRES so that 
 *                we can use a variable preconditioning step.
 *   2016-11-02 : Add restarted GMRES method.
 */

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

import nm.smla;
import nm.bbla;

import block;
import sblock;
import globaldata;
import globalconfig;
import simcore;
import fvcore;
import fileutil;
import user_defined_source_terms;

immutable double dU_essentially_zero = 1.0e-16;

int fnCount = 0;

struct ResidValues {
    double mass;
    double xMomentum;
    double yMomentum;
    double zMomentum;
    double energy;
}

void main(string[] args)
{
    //    GlobalConfig.use_steady_state_solver = true;

    writeln("Eilmer compressible-flow simulation code -- steady state solver.");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                               Comment:\n";
    msg       ~= "e4sss    [--job=<string>]            file names built from this string\n";
    msg       ~= "         [--verbosity=<int>]         defaults to 0\n";
    msg       ~= "\n";
    msg       ~= "         [--restart]                 start from latest available iteration\n";
    msg       ~= "         [--max-cpus=<int>]          defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine\n";
    msg       ~= "         [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "\n";
    msg       ~= "         [--help]                    writes this message\n";
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    int verbosityLevel = 0;
    int tindxStart = 0;
    bool restart = false;
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;
    try {
	getopt(args,
	       "job", &jobName,
	       "verbosity", &verbosityLevel,
	       "restart", &restart,
	       "max-cpus", &maxCPUs,
	       "max-wall-clock", &maxWallClock,
	       "help", &helpWanted
	       );
    } catch (Exception e) {
	writeln("Problem parsing command-line options.");
	writeln("Arguments not processed: ");
	args = args[1 .. $]; // Dispose of program name in first argument
	foreach (arg; args) writeln("   arg: ", arg);
	write(msg);
	exit(1);
    }
    if (helpWanted) {
	write(msg);
	exit(0);
    }

    if (jobName.length == 0) {
	writeln("Need to specify a job name.");
	write(msg);
	exit(1);
    }

    GlobalConfig.base_file_name = jobName;
    GlobalConfig.verbosity_level = verbosityLevel;
    maxCPUs = min(max(maxCPUs, 1), totalCPUs); // don't ask for more than available

    if (restart) {
	string errMsg = "Restart option not yet implemented.";
	throw new Error(errMsg);
	/*
	auto times_dict = readTimesFile(jobName);
	auto tindx_list = times_dict.keys;
	sort(tindx_list);
	tindxStart = tindx_list[$-1];
	*/
    }

    if (verbosityLevel > 0) {
	writeln("Begin simulation with command-line arguments.");
	writeln("  jobName: ", jobName);
	writeln("  restart: ", restart);
	writeln("  maxWallClock: ", maxWallClock);
	writeln("  verbosityLevel: ", verbosityLevel);
	writeln("  maxCPUs: ", maxCPUs);
    }
	
    init_simulation(tindxStart, maxCPUs, maxWallClock);
    // Additional initialisation
    foreach (blk; gasBlocks) {
	blk.allocate_GMRES_workspace();
    }

    /* Check that items are implemented. */
    bool goodToProceed = true;
    if ( GlobalConfig.dimensions == 3 ) {
	writeln("Steady-state solver not implemented for 3D calculations.");
	goodToProceed = false;
    }
    //if ( GlobalConfig.viscous == true ) {
    //	writeln("Steady-state solver not implemented for viscous calculations.");
    //	goodToProceed = false;
    //    }
    if ( gasBlocks.length > 1 ) {
	writeln("Steady-state solver not implemented nBlocks > 1.");
	goodToProceed = false;
    }

    if ( !goodToProceed ) {
	writeln("One or more options are not yet available for the steady-state solver.");
	writeln("Bailing out!");
	exit(1);
    }

    iterate_to_steady_state();
    writeln("Done simulation.");
}


void iterate_to_steady_state()
{
    bool withPTC = true;
    string jobName = GlobalConfig.base_file_name;
    int nsteps = GlobalConfig.sssOptions.nNewtonSteps;
    int nOuterIterations;
    int maxNumberAttempts = GlobalConfig.sssOptions.maxNumberAttempts;
    int nConserved = GlobalConfig.sssOptions.nConserved;
    double cflInit = GlobalConfig.sssOptions.cflInit;
    SBlock sblk = cast(SBlock) gasBlocks[0];
    int interpOrderSave = sblk.get_interpolation_order();

    double dt = determine_initial_dt(cflInit);
    double dtTrial;
    bool failedAttempt;
    double pseudoSimTime = 0.0;
    double normOld, normNew;
    int snapshotsFreq = GlobalConfig.sssOptions.snapshotsFrequency;
    int snapshotsCount = GlobalConfig.sssOptions.snapshotsCount;
    int writtenSnapshotsCount = 0;

    int nPreSteps = GlobalConfig.sssOptions.nPreSteps;
    int nLowOrderSteps = GlobalConfig.sssOptions.nLowOrderSteps;
    double tau = GlobalConfig.sssOptions.tau;
    bool inexactNewtonPhase = false;

    double[string][] times;
    times ~= [ "pst":0.0, "dt":dt ];

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
    double normMax = 0.0;
    ResidValues maxResid, currResid;

    bool withPreconditioning = false;

    writeln("Begin pre-iterations to establish sensible max residual.");
    foreach ( preStep; -nPreSteps .. 0 ) {
	sblk.set_interpolation_order(1);
	foreach (attempt; 0 .. maxNumberAttempts) {
	    FGMRES_solve(pseudoSimTime, dt, withPreconditioning, withPTC, normOld, nOuterIterations);
	    foreach (blk; gasBlocks) {
		int cellCount = 0;
		foreach (cell; blk.cells) {
		    cell.U[1].copy_values_from(cell.U[0]);
		    cell.U[1].mass = cell.U[0].mass + blk.dU[cellCount+0];
		    cell.U[1].momentum.refx = cell.U[0].momentum.x + blk.dU[cellCount+1];
		    cell.U[1].momentum.refy = cell.U[0].momentum.y + blk.dU[cellCount+2];
		    cell.U[1].total_energy = cell.U[0].total_energy + blk.dU[cellCount+3];
		    try {
			cell.decode_conserved(0, 1, 0.0);
		    }
		    catch (FlowSolverException e) {
			writefln("Failed attempt %d: dt= %e", attempt, dt);
			failedAttempt = true;
			dt = 0.1*dt;
			break;
		    }
		    cellCount += nConserved;
		}
	    }

	    if ( failedAttempt )
		continue;

	    // If we get here, things are good. Put flow state into U[0]
	    // ready for next iteration.
	    foreach (blk; gasBlocks) {
		foreach (cell; blk.cells) {
		    swap(cell.U[0], cell.U[1]);
		}
	    }
	    // If we got here, we can break the attempts loop.
	    break;
	}
	if ( failedAttempt ) {
	    writefln("Pre-step failed: %d", step);
	    writeln("Bailing out!");
	    exit(1);
	}
	
	writefln("PRE-STEP  %d  ::  dt= %.3e  global-norm= %.12e", preStep, dt, normOld); 
	if ( normOld > normMax ) {
	    normMax = normOld;
	    max_residuals(gasBlocks[0], maxResid);
	}
    }

    writeln("Pre-step phase complete.");
    
    if ( nPreSteps <= 0 ) {
	// Take initial residual as max residual
	evalRHS(gasBlocks[0], 0.0, 0);
	max_residuals(gasBlocks[0], maxResid);
    }

    writeln("Maximum resisudals are established as:");
    writefln("MASS:           %.12e", maxResid.mass);
    writefln("X-MOMENTUM:     %.12e", maxResid.xMomentum);
    writefln("Y-MOMENTUM:     %.12e", maxResid.yMomentum);
    writefln("ENERGY:         %.12e", maxResid.energy);

    // Open file for writing residuals
    auto residFname = "e4sss.diagnostics.dat";
    auto fResid = File(residFname, "w");
    fResid.write("# step   pseudo-time    dt   nOuterIterations  nFnCalls mass-abs   mass-rel    x-mom-abs  x-mom-rel   y-mom-abs  y-mom-rel   energy-abs  energy-rel\n");
    fResid.close();

    // Begin Newton steps
    foreach (step; 1 .. nsteps+1) {
	if ( step <= nLowOrderSteps ) {
	    sblk.set_interpolation_order(1);
	    withPreconditioning = false;
	}
	else {
	    sblk.set_interpolation_order(interpOrderSave);
	    withPreconditioning = GlobalConfig.sssOptions.usePreconditioning;
	}

	foreach (attempt; 0 .. maxNumberAttempts) {
	    failedAttempt = false;
	    FGMRES_solve(pseudoSimTime, dt, withPreconditioning, withPTC, normNew, nOuterIterations);
	    foreach (blk; gasBlocks) {
		int cellCount = 0;
		foreach (cell; blk.cells) {
		    cell.U[1].copy_values_from(cell.U[0]);
		    cell.U[1].mass = cell.U[0].mass + blk.dU[cellCount+0];
		    cell.U[1].momentum.refx = cell.U[0].momentum.x + blk.dU[cellCount+1];
		    cell.U[1].momentum.refy = cell.U[0].momentum.y + blk.dU[cellCount+2];
		    cell.U[1].total_energy = cell.U[0].total_energy + blk.dU[cellCount+3];
		    try {
			cell.decode_conserved(0, 1, 0.0);
		    }
		    catch (FlowSolverException e) {
			writefln("Failed attempt %d: dt= %e", attempt, dt);
			failedAttempt = true;
			dt = 0.1*dt;
			break;
		    }
		    cellCount += nConserved;
		}
	    }

	    if ( failedAttempt )
		continue;

	    // If we get here, things are good. Put flow state into U[0]
	    // ready for next iteration.
	    foreach (blk; gasBlocks) {
		foreach (cell; blk.cells) {
		    swap(cell.U[0], cell.U[1]);
		}
	    }
	    // If we got here, we can break the attempts loop.
	    break;
	}
	if ( failedAttempt ) {
	    writefln("Step failed: %d", step);
	    writeln("Bailing out!");
	    exit(1);
	}

	pseudoSimTime += dt;
	
	// Now do some output and diagnostics work
	// Write out residuals
	max_residuals(gasBlocks[0], currResid);
	fResid = File(residFname, "a");
	fResid.writef("%8d  %20.16e  %20.16e  %3d %5d %20.16e  %20.16e  %20.16e %20.16e  %20.16e  %20.16e  %20.16e  %20.16e\n",
		      step, pseudoSimTime, dt, nOuterIterations, fnCount,
		      currResid.mass, currResid.mass/maxResid.mass,
		      currResid.xMomentum, currResid.xMomentum/maxResid.xMomentum,
		      currResid.yMomentum, currResid.yMomentum/maxResid.yMomentum,
		      currResid.energy, currResid.energy/maxResid.energy);
	fResid.close();

	if ( (step % GlobalConfig.print_count) == 0 ) {
	    auto writer = appender!string();
	    formattedWrite(writer, "Iteration= %7d  pseudo-time=%10.3e dt=%10.3e global-residual= %10.6e\n", step, pseudoSimTime, dt, normNew);
	    formattedWrite(writer, "Residuals:  mass  x-mom  y-mom  energy\n");
	    formattedWrite(writer, "absolute   %10.6e %10.6e %10.6e %10.6e\n",
			   currResid.mass, currResid.xMomentum, currResid.yMomentum, currResid.energy);
	    formattedWrite(writer, "relative   %10.6e %10.6e %10.6e %10.6e\n",
			   currResid.mass/maxResid.mass,
			   currResid.xMomentum/maxResid.xMomentum,
			   currResid.yMomentum/maxResid.yMomentum,
			   currResid.energy/maxResid.energy);
	    writeln(writer.data);
	}

	// Write out the flow field, if required
	if ( (step % snapshotsFreq) == 0 ) {
	    writefln("Writing flow solution at pseudo-time= %e", pseudoSimTime);
	    writtenSnapshotsCount++;
	    if ( writtenSnapshotsCount <= snapshotsCount ) {
		ensure_directory_is_present(make_path_name!"flow"(writtenSnapshotsCount));
		foreach ( iblk, blk; gasBlocks ) {
		    auto fileName = make_file_name!"flow"(jobName, to!int(iblk), writtenSnapshotsCount);
		    blk.write_solution(fileName, pseudoSimTime);
		}
		times ~= [ "pst": pseudoSimTime, "dt":dt ];
		rewrite_times_file(times);
	    }
	    else {
		// We need to shuffle all of the snapshots...
		foreach ( iSnap; 2 .. snapshotsCount+1) {
		    foreach ( iblk; 0 .. gasBlocks.length ) {
			auto fromName = make_file_name!"flow"(jobName, to!int(iblk), iSnap);
			auto toName = make_file_name!"flow"(jobName, to!int(iblk), iSnap-1);
			rename(fromName, toName);
		    }
		}
		// ... and add the new snapshot.
		foreach ( iblk, blk; gasBlocks ) {
		    auto fileName = make_file_name!"flow"(jobName, to!int(iblk), snapshotsCount);	
		    blk.write_solution(fileName, pseudoSimTime);
		}
		remove(times, 1);
		times[$-1] = [ "pst": pseudoSimTime, "dt":dt ];
		rewrite_times_file(times);
	    }
	}

	if ( !inexactNewtonPhase && currResid.mass/maxResid.mass < tau ) {
	    // Switch to inexactNewtonPhase
	    inexactNewtonPhase = true;
	}

	// Choose a new timestep.
	if ( inexactNewtonPhase ) {
	    if ( step <= nLowOrderSteps ) {
		// Let's assume we're still letting the shock settle
		// when doing low order steps, so we use a power of 0.75
		dtTrial = dt*pow(normOld/normNew, 0.75);
	    }
	    else {
		// We use a power of 1.0
		dtTrial = dt*(normOld/normNew);
	    }
	    dtTrial = min(dtTrial, 2.0*dt);
	    dtTrial = max(dtTrial, 0.1*dt);
	    dt = dtTrial;
	}

	normOld = normNew;
    }

}

double determine_initial_dt(double cflInit)
{
    double signal, dt_local, dt;
    bool first = true;

    foreach (blk; gasBlocks) {
	foreach (cell; blk.cells) {
	    signal = cell.signal_frequency();
	    dt_local = cflInit / signal;
	    if ( first )
		dt = dt_local;
	    else
		dt = fmin(dt, dt_local);
	}
    }
    return dt;
}

void evalRHS(Block blk, double pseudoSimTime, int ftl)
{
    fnCount++;

    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    blk.applyPreReconAction(pseudoSimTime, 0, ftl);
    blk.convective_flux_phase0();
    blk.convective_flux_phase1();
    blk.applyPostConvFluxAction(pseudoSimTime, 0, ftl);
    if (GlobalConfig.viscous) {
	blk.applyPreSpatialDerivAction(pseudoSimTime, 0, ftl);
	blk.flow_property_derivatives(0); 
	blk.estimate_turbulence_viscosity();
	blk.viscous_flux();
	blk.applyPostDiffFluxAction(pseudoSimTime, 0, ftl);
    }

    foreach (i, cell; blk.cells) {
	cell.add_inviscid_source_vector(0, 0.0);
	if (blk.myConfig.udf_source_terms) {
	    addUDFSourceTermsToCell(blk.myL, cell, 0, 
				    pseudoSimTime, blk.myConfig.gmodel);
	}
	cell.time_derivatives(0, ftl, false);
    }
}

void evalJacobianVecProd(Block blk, double pseudoSimTime, double[] v)
{
    double sigma = GlobalConfig.sssOptions.sigma;
    int nConserved = GlobalConfig.sssOptions.nConserved;
    SBlock sblk = cast(SBlock) gasBlocks[0];
    int interpOrder = sblk.get_interpolation_order();
    // We perform a Frechet derivative to evaluate J*v
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    int cellCount = 0;
    foreach (cell; blk.cells) {
	cell.U[1].copy_values_from(cell.U[0]);
	cell.U[1].mass += sigma*blk.S[cellCount+0]*v[cellCount+0];
	cell.U[1].momentum.refx += sigma*blk.S[cellCount+1]*v[cellCount+1];
	cell.U[1].momentum.refy += sigma*blk.S[cellCount+2]*v[cellCount+2];
	cell.U[1].total_energy += sigma*blk.S[cellCount+3]*v[cellCount+3];
	cell.decode_conserved(0, 1, 0.0);
	cellCount += nConserved;
    }
    evalRHS(blk, pseudoSimTime, 1);
    cellCount = 0;
    foreach (cell; blk.cells) {
	v[cellCount+0] = (-cell.dUdt[1].mass - blk.FU[cellCount+0])/(sigma*blk.S[cellCount+0]);
	v[cellCount+1] = (-cell.dUdt[1].momentum.x - blk.FU[cellCount+1])/(sigma*blk.S[cellCount+1]);
	v[cellCount+2] = (-cell.dUdt[1].momentum.y - blk.FU[cellCount+2])/(sigma*blk.S[cellCount+2]);
	v[cellCount+3] = (-cell.dUdt[1].total_energy - blk.FU[cellCount+3])/(sigma*blk.S[cellCount+3]);
	cellCount += nConserved;
    }
}


void FGMRES_solve(double pseudoSimTime, double dt, bool withPreconditioning, bool withPTC, ref double residual, ref int nRestarts)
{
    double resid;
    int nConserved = GlobalConfig.sssOptions.nConserved;
    // Presently, just do one block
    Block blk = gasBlocks[0];
    SBlock sblk = cast(SBlock) gasBlocks[0];
    int interpOrder = sblk.get_interpolation_order();
    int maxIters = GlobalConfig.sssOptions.maxOuterIterations;
    // We add 1 because the user thinks of "re"starts, so they
    // might legitimately ask for no restarts. We still have
    // to execute at least once.
    int maxRestarts = GlobalConfig.sssOptions.maxRestarts + 1; 
    size_t m = to!size_t(maxIters);
    size_t n = blk.v_outer.length;
    size_t r;
    size_t iterCount;

    // Variables for max rates of change
    // Use these for equation scaling.
    double maxMass = 0.0;
    double maxMom = 0.0;
    double maxEnergy = 0.0;
    double minNonDimVal = 1.0; // minimum value used for non-dimensionalisation
                               // when our time rates of change are very small
                               // then we'll avoid non-dimensionalising by 
                               // values close to zero.

    // 1. Evaluate r0, beta, v1
    evalRHS(blk, pseudoSimTime, 0);
    // Store dUdt[0] as F(U)
    int cellCount = 0;
    foreach (cell; blk.cells) {
	blk.FU[cellCount+0] = -cell.dUdt[0].mass;
	blk.FU[cellCount+1] = -cell.dUdt[0].momentum.x;
	blk.FU[cellCount+2] = -cell.dUdt[0].momentum.y;
	blk.FU[cellCount+3] = -cell.dUdt[0].total_energy;
	cellCount += nConserved;
	maxMass = fmax(maxMass, cell.dUdt[0].mass);
	maxMom = max(maxMom, cell.dUdt[0].momentum.x, cell.dUdt[0].momentum.y);
	maxEnergy = fmax(maxEnergy, cell.dUdt[0].total_energy);
    }
    double unscaledNorm2 = norm2(blk.FU);
    // Place some guards when time-rate-of-changes are very small.
    maxMass = fmax(maxMass, minNonDimVal);
    maxMom = fmax(maxMom, minNonDimVal);
    maxEnergy = fmax(maxEnergy, minNonDimVal);
    // Now set-up scaling matrix (which is diagonal, so just store as a vector)
    cellCount = 0;
    foreach (cell; blk.cells) {
	blk.S[cellCount+0] = maxMass;
	blk.S[cellCount+1] = maxMom;
	blk.S[cellCount+2] = maxMom;
	blk.S[cellCount+3] = maxEnergy;
	cellCount += nConserved;
    }

    // Initialise some arrays and matrices that have already been allocated
    blk.g0_outer[] = 0.0;
    blk.g1_outer[] = 0.0;
    blk.H0_outer.zeros();
    blk.H1_outer.zeros();

    // We'll scale r0 against these max rates of change.
    // r0 = b - A*x0
    // Taking x0 = [0] (as is common) gives r0 = b = FU
    blk.x0[] = 0.0;
    cellCount = 0;
    foreach (cell; blk.cells) {
	blk.r0[cellCount+0] = -(1./blk.S[cellCount+0])*blk.FU[cellCount+0];
	blk.r0[cellCount+1] = -(1./blk.S[cellCount+1])*blk.FU[cellCount+1];
	blk.r0[cellCount+2] = -(1./blk.S[cellCount+2])*blk.FU[cellCount+2];
	blk.r0[cellCount+3] = -(1./blk.S[cellCount+3])*blk.FU[cellCount+3];
	cellCount += nConserved;
    }
    // Then compute v = r0/||r0||
    auto beta = norm2(blk.r0);
    // DEBUG:  writefln("OUTER: beta= %e", beta);
    blk.g0_outer[0] = beta;
    foreach (k; 0 .. n) {
	blk.v_outer[k] = blk.r0[k]/beta;
	blk.V_outer[k,0] = blk.v_outer[k];
    }

    // Compute tolerance
    auto outerTol = GlobalConfig.sssOptions.eta*beta;
    // DEBUG:  writefln("OUTER: outerTol= %e", outerTol);

    // 2. Start outer-loop of restarted GMRES
    for ( r = 0; r < maxRestarts; r++ ) {
	// 2a. Begin iterations
	foreach (j; 0 .. m) {
	    iterCount = j+1;
	    if ( withPreconditioning ) {
		GMRES_solve(pseudoSimTime, dt, withPTC);
	    }
	    else {
		blk.z_outer[] = blk.v_outer[];
	    }
	    // Save z vector in Z.
	    foreach (k; 0 .. n) blk.Z_outer[k,j] = blk.z_outer[k];

	    // Prepare 'w' with (I/dt)v term;
	    double dtInv = 1.0/dt;
	    foreach (k; 0 .. n) {
		blk.w_outer[k] = dtInv*blk.z_outer[k];
	    }
	
	    // Evaluate Jz and place result in z_outer
	    evalJacobianVecProd(blk, pseudoSimTime, blk.z_outer);
	    // Now we can complete calculation of w
	    foreach (k; 0 .. n)  blk.w_outer[k] += blk.z_outer[k];
	    // The remainder of the algorithm looks a lot like any standard
	    // GMRES implementation (for example, see smla.d)
	    foreach (i; 0 .. j+1) {
		foreach (k; 0 .. n ) blk.v_outer[k] = blk.V_outer[k,i]; // Extract column 'i'
		blk.H0_outer[i,j] = dot(blk.w_outer, blk.v_outer);
		foreach (k; 0 .. n) blk.w_outer[k] -= blk.H0_outer[i,j]*blk.v_outer[k]; 
	    }
	    blk.H0_outer[j+1,j] = norm2(blk.w_outer);
	
	    foreach (k; 0 .. n) {
		blk.v_outer[k] = blk.w_outer[k]/blk.H0_outer[j+1,j];
		blk.V_outer[k,j+1] = blk.v_outer[k];
	    }

	    // Build rotated Hessenberg progressively
	    if ( j != 0 ) {
		// Extract final column in H
		foreach (i; 0 .. j+1) blk.h_outer[i] = blk.H0_outer[i,j];
		// Rotate column by previous rotations (stored in Q0)
		nm.bbla.dot(blk.Q0_outer, j+1, j+1, blk.h_outer, blk.hR_outer);
		// Place column back in H
		foreach (i; 0 .. j+1) blk.H0_outer[i,j] = blk.hR_outer[i];
	    }
	    // Now form new Gamma
	    blk.Gamma_outer.eye();
	    auto denom = sqrt(blk.H0_outer[j,j]*blk.H0_outer[j,j] + blk.H0_outer[j+1,j]*blk.H0_outer[j+1,j]);
	    auto s_j = blk.H0_outer[j+1,j]/denom; 
	    auto c_j = blk.H0_outer[j,j]/denom;
	    blk.Gamma_outer[j,j] = c_j; blk.Gamma_outer[j,j+1] = s_j;
	    blk.Gamma_outer[j+1,j] = -s_j; blk.Gamma_outer[j+1,j+1] = c_j;
	    // Apply rotations
	    nm.bbla.dot(blk.Gamma_outer, j+2, j+2, blk.H0_outer, j+1, blk.H1_outer);
	    nm.bbla.dot(blk.Gamma_outer, j+2, j+2, blk.g0_outer, blk.g1_outer);
	    // Accumulate Gamma rotations in Q.
	    if ( j == 0 ) {
		copy(blk.Gamma_outer, blk.Q1_outer);
	    }
	    else {
		nm.bbla.dot(blk.Gamma_outer, j+2, j+2, blk.Q0_outer, j+2, blk.Q1_outer);
	    }

	    // Prepare for next step
	    copy(blk.H1_outer, blk.H0_outer);
	    blk.g0_outer[] = blk.g1_outer[];
	    copy(blk.Q1_outer, blk.Q0_outer);

	    // Get residual
	    resid = fabs(blk.g1_outer[j+1]);
	    // DEBUG:
	    //	    writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
	    if ( resid <= outerTol ) {
		m = j+1;
		// DEBUG:
		//	writefln("OUTER: TOL ACHIEVED restart-count= %d iteration-count= %d, resid= %e", r, m, resid);
		break;
	    }

	}
	if ( iterCount == maxIters )
	    m = maxIters;

	// At end H := R up to row m
	//        g := gm up to row m
	upperSolve(blk.H1_outer, to!int(m), blk.g1_outer);
	nm.bbla.dot(blk.Z_outer, n, m, blk.g1_outer, blk.dU);
	foreach (k; 0 .. n) blk.dU[k] += blk.x0[k];

	if ( resid <= outerTol || r+1 == maxRestarts ) {
	    // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
	    // DEBUG:  writefln("Breaking restart loop.");
	    break;
	}

	// Else, we prepare for restart by computing r0 and setting x0
	blk.x0[] = blk.dU[];
	foreach (k; 0 .. n) blk.Z_outer[k,m] = blk.V_outer[k,m]; // store final Arnoldi vector in Z (unpreconditioned)
	nm.bbla.dot(blk.Z_outer, n, m+1, blk.Q1_outer, m+1, blk.V_outer);
	foreach (i; 0 .. m) blk.g0_outer[i] = 0.0;
	nm.bbla.dot(blk.V_outer, n, m+1, blk.g0_outer, blk.r0);

	beta = norm2(blk.r0);
	// DEBUG: writefln("OUTER: ON RESTART beta= %e", beta);
	foreach (k; 0 .. n) {
	    blk.v_outer[k] = blk.r0[k]/beta;
	    blk.V_outer[k,0] = blk.v_outer[k];
	}
	// Re-initialise some vectors and matrices for restart
	blk.g0_outer[] = 0.0;
	blk.g1_outer[] = 0.0;
	blk.H0_outer.zeros();
	blk.H1_outer.zeros();
	// And set first residual entry
	blk.g0_outer[0] = beta;
    }
    // Rescale
    cellCount = 0;
    foreach (cell; blk.cells) {
	blk.dU[cellCount+0] *= blk.S[cellCount+0];
	blk.dU[cellCount+1] *= blk.S[cellCount+1];
	blk.dU[cellCount+2] *= blk.S[cellCount+2];
	blk.dU[cellCount+3] *= blk.S[cellCount+3];
	cellCount += nConserved;
    }
    
    residual = unscaledNorm2;
    nRestarts = to!int(m);
}

void GMRES_solve(double pseudoSimTime, double dt, bool withPTC)
{
    // We always perform 'nInnerIterations' for the preconditioning step
    // That is, we do NOT stop on some tolerance in this solve.
    double resid;
    int nConserved = GlobalConfig.sssOptions.nConserved;
    // Presently, just do one block
    Block blk = gasBlocks[0];
    int maxIters = GlobalConfig.sssOptions.nInnerIterations;
    size_t m = to!size_t(maxIters);
    size_t n = blk.v_inner.length;
    size_t iterCount;

    // 1. Evaluate r0, beta, v1
    blk.g0_inner[] = 0.0;
    blk.g1_inner[] = 0.0;
    blk.H0_inner.zeros();
    blk.H1_inner.zeros();
    blk.Gamma_inner.eye();

    // r0 = v_outer - A x0, with x0 = [0]
    // Then compute v = r0/||r0||
    auto beta = norm2(blk.v_outer);
    blk.g0_inner[0] = beta;

    foreach (k; 0 .. n) {
	blk.v_inner[k] = blk.v_outer[k]/beta;
	blk.V_inner[k,0] = blk.v_inner[k];
    }

    // 2. Begin iterations
    foreach (j; 0 .. m) {
	iterCount = j+1;
	
	// Prepare 'w' with (I/dt)v term;
	double dtInv = 1.0/dt;
	foreach (k; 0 .. n) {
	    blk.w_inner[k] = dtInv*blk.v_inner[k];
	}
	// Evaluate Jv and place in v
	evalJacobianVecProd(blk, pseudoSimTime, blk.v_inner);
	// Complete calculation of w.
	foreach (k; 0 .. n) blk.w_inner[k] += blk.v_inner[k];

	// The remainder of the algorithm looks a lot like any standard
	// GMRES implementation (for example, see smla.d)
	foreach (i; 0 .. j+1) {
	    foreach (k; 0 .. n ) blk.v_inner[k] = blk.V_inner[k,i]; // Extract column 'i'
	    blk.H0_inner[i,j] = dot(blk.w_inner, blk.v_inner);
	    foreach (k; 0 .. n) blk.w_inner[k] -= blk.H0_inner[i,j]*blk.v_inner[k]; 
	}
	blk.H0_inner[j+1,j] = norm2(blk.w_inner);
	
	foreach (k; 0 .. n) {
	    blk.v_inner[k] = blk.w_inner[k]/blk.H0_inner[j+1,j];
	    blk.V_inner[k,j+1] = blk.v_inner[k];
	}

	// Build rotated Hessenberg progressively
	if ( j != 0 ) {
	    // Extract final column in H
	    foreach (i; 0 .. j+1) blk.h_inner[i] = blk.H0_inner[i,j];
	    // Rotate column by previous rotations (stored in Q0)
	    nm.bbla.dot(blk.Q0_inner, j+1, j+1, blk.h_inner, blk.hR_inner);
	    // Place column back in H
	    foreach (i; 0 .. j+1) blk.H0_inner[i,j] = blk.hR_inner[i];
	}
	// Now form new Gamma
	blk.Gamma_inner.eye();
	auto denom = sqrt(blk.H0_inner[j,j]*blk.H0_inner[j,j] + blk.H0_inner[j+1,j]*blk.H0_inner[j+1,j]);
	auto s_j = blk.H0_inner[j+1,j]/denom; 
	auto c_j = blk.H0_inner[j,j]/denom;
	blk.Gamma_inner[j,j] = c_j; blk.Gamma_inner[j,j+1] = s_j;
	blk.Gamma_inner[j+1,j] = -s_j; blk.Gamma_inner[j+1,j+1] = c_j;
	// Apply rotations
	nm.bbla.dot(blk.Gamma_inner, j+2, j+2, blk.H0_inner, j+1, blk.H1_inner);
	nm.bbla.dot(blk.Gamma_inner, j+2, j+2, blk.g0_inner, blk.g1_inner);
	// Accumulate Gamma rotations in Q.
	if ( j == 0 ) {
	    copy(blk.Gamma_inner, blk.Q1_inner);
	}
	else {
	    nm.bbla.dot(blk.Gamma_inner, j+2, j+2, blk.Q0_inner, j+2, blk.Q1_inner);
	}

	// Prepare for next step
	copy(blk.H1_inner, blk.H0_inner);
	blk.g0_inner[] = blk.g1_inner[];
	copy(blk.Q1_inner, blk.Q0_inner);
    }

    m = maxIters;
    // At end H := R up to row m
    //        g := gm up to row m
    upperSolve(blk.H1_inner, to!int(m), blk.g1_inner);
    nm.bbla.dot(blk.V_inner, n, m, blk.g1_inner, blk.z_outer);

}

void max_residuals(Block blk, ref ResidValues rv)
{
    int nc = GlobalConfig.sssOptions.nConserved;
    double massMax = blk.cells[0].dUdt[0].mass;
    double xMomMax = blk.cells[0].dUdt[0].momentum.x;
    double yMomMax = blk.cells[0].dUdt[0].momentum.y;
    double energyMax = blk.cells[0].dUdt[0].total_energy;

    double massLocal, xMomLocal, yMomLocal, energyLocal;
    foreach (cell; blk.cells) {
	massLocal = cell.dUdt[0].mass;
	xMomLocal = cell.dUdt[0].momentum.x;
	yMomLocal = cell.dUdt[0].momentum.y;
	energyLocal = cell.dUdt[0].total_energy;

	massMax = fmax(massMax, massLocal);
	xMomMax = fmax(xMomMax, xMomLocal);
	yMomMax = fmax(yMomMax, yMomLocal);
	energyMax = fmax(energyMax, energyLocal);
    }

    rv.mass = massMax;
    rv.xMomentum = xMomMax;
    rv.yMomentum = yMomMax;
    rv.energy = energyMax;
}

void rewrite_times_file(double[string][] times)
{
    auto fname = format("%s.times", GlobalConfig.base_file_name);
    auto f = File(fname, "w");
    f.writeln("# tindx sim_time dt_global");
    foreach (i, t; times) {
	f.writefln("%04d %.18e %.18e", i, t["pst"], t["dt"]);
    }
    f.close();
}
