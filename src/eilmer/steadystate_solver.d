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
    string jobName = GlobalConfig.base_file_name;
    int nsteps = GlobalConfig.sssOptions.nOuterIterations;
    int nInnerIterations;
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

    int nPreIterations = GlobalConfig.sssOptions.nPreIterations;
    int nLowOrderSteps = GlobalConfig.sssOptions.nLowOrderIterations;
    double tau = GlobalConfig.sssOptions.tau;
    bool inexactNewtonPhase = false;

    double[string][] times;
    times ~= [ "pst":0.0, "dt":dt ];

    // The initial residual is usually a poor choice for basing decisions about how the
    // residual is dropping, particularly when a constant initial condition is given.
    // A constant initial condition gives a zero residual everywhere in the interior
    // of the domain. Therefore, the only residual can come from boundary condition influences.
    // What we do here is use some early iterations (I've called them "pre"-iterations) to
    // allow the boundary conditions to exert some influence into the interior.
    // We'll find the max residual in that start-up process and use that as our initial value
    // for residual. Our relative residuals will then be based on that.
    // We can use a fixed timestep (typically small) and low-order reconstruction in this
    // pre-iteration stage.
    double normMax = 0.0;
    ResidValues maxResid, currResid;


    writeln("Begin pre-iterations to establish sensible max residual.");
    foreach ( preStep; -nPreIterations .. 0 ) {
	sblk.set_interpolation_order(1);
	foreach (attempt; 0 .. maxNumberAttempts) {
	    GMRES_solve(pseudoSimTime, dt, normOld, nInnerIterations);
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
	
	writefln("pre-iteration %d:  dt= %.3e  global-norm= %.12e", preStep, dt, normOld); 
	if ( normOld > normMax ) {
	    normMax = normOld;
	    max_residuals(gasBlocks[0], maxResid);
	}
    }

    writeln("Pre-iteration phase complete.");
    writeln("Maximum resisudals are established as:");
    writefln("MASS:           %.12e", maxResid.mass);
    writefln("X-MOMENTUM:     %.12e", maxResid.xMomentum);
    writefln("Y-MOMENTUM:     %.12e", maxResid.yMomentum);
    writefln("ENERGY:         %.12e", maxResid.energy);

    // Open file for writing residuals
    auto fname = "e4sss.diagnostics.dat";
    auto fResid = File(fname, "w");
    fResid.write("# iteration   pseudo-time    dt     mass-abs   mass-rel    x-mom-abs  x-mom-rel   y-mom-abs  y-mom-rel   energy-abs  energy-rel\n");

    // Begin Newton iterations
    foreach (step; 1 .. nsteps+1) {
	if ( step <= nLowOrderSteps ) {
	    sblk.set_interpolation_order(1);
	}
	else {
	    sblk.set_interpolation_order(interpOrderSave);
	}


	foreach (attempt; 0 .. maxNumberAttempts) {
	    failedAttempt = false;
	    GMRES_solve(pseudoSimTime, dt, normNew, nInnerIterations);
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
	fResid.writef("%8d  %20.16e  %20.16e  %3d %20.16e  %20.16e  %20.16e %20.16e  %20.16e  %20.16e  %20.16e  %20.16e\n",
		      step, pseudoSimTime, dt, nInnerIterations, 
		      currResid.mass, currResid.mass/maxResid.mass,
		      currResid.xMomentum, currResid.xMomentum/maxResid.xMomentum,
		      currResid.yMomentum, currResid.yMomentum/maxResid.yMomentum,
		      currResid.energy, currResid.energy/maxResid.energy);

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
		    //writeln("fileName= ", fileName);
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
			//writeln("Shuffle: ", fromName, " --> ", toName);
		    }
		}
		// ... and add the new snapshot.
		foreach ( iblk, blk; gasBlocks ) {
		    auto fileName = make_file_name!"flow"(jobName, to!int(iblk), snapshotsCount);	
		    //writeln("fileName= ", fileName);
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
	    dtTrial = dt*pow(normOld/normNew, 0.75);
	    //writefln("DEBUG: normNew=%e normOld=%e dt= %e  dtTrial= %e", normNew, normOld, dt, dtTrial);
	    dtTrial = min(dtTrial, 2.0*dt);
	    dtTrial = max(dtTrial, 0.1*dt);
	    dt = dtTrial;
	    //writefln("DEBUG: dt chosen= %e", dt);
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
    int nConserved = GlobalConfig.sssOptions.nConserved;
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

    foreach (cell; blk.cells) {
	cell.add_inviscid_source_vector(0, 0.0);
	if (blk.myConfig.udf_source_terms) {
	    addUDFSourceTermsToCell(blk.myL, cell, 0, 
				    pseudoSimTime, blk.myConfig.gmodel);
	}
	cell.time_derivatives(0, ftl, false);
    }
}

void evalJacobianVecProd(Block blk, double pseudoSimTime)
{
    double sigma = GlobalConfig.sssOptions.sigma;
    int nConserved = GlobalConfig.sssOptions.nConserved;
    // We perform a Frechet derivative to evaluate J*v
    blk.clear_fluxes_of_conserved_quantities();
    foreach (cell; blk.cells) cell.clear_source_vector();
    int cellCount = 0;
    foreach (cell; blk.cells) {
	cell.U[1].copy_values_from(cell.U[0]);
	cell.U[1].mass += sigma*blk.S[cellCount+0]*blk.v[cellCount+0];
	cell.U[1].momentum.refx += sigma*blk.S[cellCount+1]*blk.v[cellCount+1];
	cell.U[1].momentum.refy += sigma*blk.S[cellCount+2]*blk.v[cellCount+2];
	cell.U[1].total_energy += sigma*blk.S[cellCount+3]*blk.v[cellCount+3];
	cell.decode_conserved(0, 1, 0.0);
	cellCount += nConserved;
    }
    evalRHS(blk, pseudoSimTime, 1);
    cellCount = 0;
    foreach (cell; blk.cells) {
	blk.v[cellCount+0] = (cell.dUdt[1].mass - blk.FU[cellCount+0])/(sigma*blk.S[cellCount+0]);
	blk.v[cellCount+1] = (cell.dUdt[1].momentum.x - blk.FU[cellCount+1])/(sigma*blk.S[cellCount+1]);
	blk.v[cellCount+2] = (cell.dUdt[1].momentum.y - blk.FU[cellCount+2])/(sigma*blk.S[cellCount+2]);
	blk.v[cellCount+3] = (cell.dUdt[1].total_energy - blk.FU[cellCount+3])/(sigma*blk.S[cellCount+3]);
	cellCount += nConserved;
    }
}

void GMRES_solve(double pseudoSimTime, double dt, ref double residual, ref int nInnerIterations)
{
    double resid;
    int nConserved = GlobalConfig.sssOptions.nConserved;
    // Presently, just do one block
    Block blk = gasBlocks[0];
    int maxIters = GlobalConfig.sssOptions.maxInnerIterations;
    size_t m = to!size_t(maxIters);
    size_t n = blk.v.length;
    size_t iterCount;
    double sigma = 1.0e-8;
    // Initialise some arrays and matrices that have already been allocated
    blk.g0[] = 0.0;
    blk.g1[] = 0.0;
    blk.H0.zeros();
    blk.H1.zeros();
    blk.Gamma.eye();

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
	blk.FU[cellCount+0] = cell.dUdt[0].mass;
	blk.FU[cellCount+1] = cell.dUdt[0].momentum.x;
	blk.FU[cellCount+2] = cell.dUdt[0].momentum.y;
	blk.FU[cellCount+3] = cell.dUdt[0].total_energy;
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
    //DEBUG: writefln("maxMass= %e  maxMom= %e maxEnergy= %e", maxMass, maxMom, maxEnergy);
    // Now set-up scaling matrix (which is diagonal, so just store as a vector)
    cellCount = 0;
    foreach (cell; blk.cells) {
	blk.S[cellCount+0] = maxMass;
	blk.S[cellCount+1] = maxMom;
	blk.S[cellCount+2] = maxMom;
	blk.S[cellCount+3] = maxEnergy;
	cellCount += nConserved;
    }

    // We'll scale r0 against these max rates of change.
    // r0 = b - A*x0
    // Taking x0 = [0] (as is common) gives r0 = b = FU
    cellCount = 0;
    foreach (cell; blk.cells) {
	blk.r0[cellCount+0] = (1./blk.S[cellCount+0])*blk.FU[cellCount+0];
	blk.r0[cellCount+1] = (1./blk.S[cellCount+1])*blk.FU[cellCount+1];
	blk.r0[cellCount+2] = (1./blk.S[cellCount+2])*blk.FU[cellCount+2];
	blk.r0[cellCount+3] = (1./blk.S[cellCount+3])*blk.FU[cellCount+3];
	cellCount += nConserved;
    }
    // Then compute v = r0/||r0||
    auto beta = norm2(blk.r0);
    // DEBUG:   writefln("beta= %e", beta);
    blk.g0[0] = beta;
    foreach (k; 0 .. n) {
	blk.v[k] = blk.r0[k]/beta;
	blk.V[k,0] = blk.v[k];
    }

    // Compute tolerance
    auto tol = GlobalConfig.sssOptions.eta*beta;
    //DEBUG: writefln("tolerance= %e", tol);

    // 2. Begin iterations
    foreach (j; 0 .. m) {
	iterCount = j+1;
	// Prepare 'w' with (I/dt)v term;
	double dtInv = 1.0/dt;
	foreach (k; 0 .. n) {
	    blk.w[k] = dtInv*blk.v[k];
	}
	// Evaluate Jv and place in v
	evalJacobianVecProd(blk, pseudoSimTime);
	// Complete calculation of w.
	foreach (k; 0 .. n) {
	    blk.w[k] -= blk.v[k];
	}
	// The remainder of the algorithm looks a lot like any standard
	// GMRES implementation (for example, see smla.d)
	foreach (i; 0 .. j+1) {
	    foreach (k; 0 .. n ) blk.v[k] = blk.V[k,i]; // Extract column 'i'
	    blk.H0[i,j] = dot(blk.w, blk.v);
	    foreach (k; 0 .. n) blk.w[k] -= blk.H0[i,j]*blk.v[k]; 
	}
	blk.H0[j+1,j] = norm2(blk.w);
	
	foreach (k; 0 .. n) {
	    blk.v[k] = blk.w[k]/blk.H0[j+1,j];
	    blk.V[k,j+1] = blk.v[k];
	}

	// Build rotated Hessenberg progressively
	if ( j != 0 ) {
	    // Extract final column in H
	    foreach (i; 0 .. j+1) blk.h[i] = blk.H0[i,j];
	    // Rotate column by previous rotations (stored in Q0)
	    nm.bbla.dot(blk.Q0, j+1, j+1, blk.h, blk.hR);
	    // Place column back in H
	    foreach (i; 0 .. j+1) blk.H0[i,j] = blk.hR[i];
	}
	// Now form new Gamma
	blk.Gamma.eye();
	auto denom = sqrt(blk.H0[j,j]*blk.H0[j,j] + blk.H0[j+1,j]*blk.H0[j+1,j]);
	auto s_j = blk.H0[j+1,j]/denom; 
	auto c_j = blk.H0[j,j]/denom;
	blk.Gamma[j,j] = c_j; blk.Gamma[j,j+1] = s_j;
	blk.Gamma[j+1,j] = -s_j; blk.Gamma[j+1,j+1] = c_j;
	// Apply rotations
	nm.bbla.dot(blk.Gamma, j+2, j+2, blk.H0, j+1, blk.H1);
	nm.bbla.dot(blk.Gamma, j+2, j+2, blk.g0, blk.g1);
	// Accumulate Gamma rotations in Q.
	if ( j == 0 ) {
	    copy(blk.Gamma, blk.Q1);
	}
	else {
	    nm.bbla.dot(blk.Gamma, j+2, j+2, blk.Q0, j+2, blk.Q1);
	}
	// Get residual
	resid = fabs(blk.g1[j+1]);
	// DEBUG: writefln("iteration= %d, resid= %e", j, resid);
	if ( resid <= tol ) {
	    m = j+1;
	    // DEBUG: writefln("iteration-count= %d, resid= %e", m, resid);
	    break;
	}

	// Prepare for next step
	copy(blk.H1, blk.H0);
	blk.g0[] = blk.g1[];
	copy(blk.Q1, blk.Q0);
    }
    if ( iterCount == maxIters )
	m = maxIters;

    // At end H := R up to row m
    //        g := gm up to row m
    upperSolve(blk.H1, to!int(m), blk.g1);
    nm.bbla.dot(blk.V, n, m, blk.g1, blk.dU);

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
    nInnerIterations = to!int(m);
}

struct ResidValues {
    double mass;
    double xMomentum;
    double yMomentum;
    double zMomentum;
    double energy;
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
