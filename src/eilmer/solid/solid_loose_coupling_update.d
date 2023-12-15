/** solid_loose_coupling_update.d
 *
 * Core set of functions used in the Newton-Krylov updates for steady-state convergence.
 *
 * Author: Kyle Damm
 * Date: 01-04-2021
 */

module solid_loose_coupling_update;

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
import ntypes.complex;
import nm.number;

import globalconfig;
import globaldata;
import solid_full_face_copy;
import solid_gas_full_face_copy;
import simcore_exchange;
import solid_udf_source_terms;
import solidprops;
import simcore_io;
import simcore;
import ssolidblock;

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

@nogc
double determine_dt(double cfl_value)
{
    double dt = double.max;
    double dt_local = double.max;
    bool first = true;
    foreach (sblk; localSolidBlocks) {
        dt_local = sblk.determine_time_step_size(cfl_value);
        if (first) {
            dt = dt_local;
            first = false;
        } else {
            dt = fmin(dt, dt_local);
        }
    }
    return dt;
} // end determine_dt

void integrate_solid_in_time_explicit(double target_time)
{
    auto wallClockStart = Clock.currTime();
    double wallClockElapsed;

    auto super_time_steps = GlobalConfig.sdluOptions.superTimeSteps;
    GlobalConfig.max_time = target_time;
    SimState.s_RKL = super_time_steps;
    SimState.time = 0.0;
    foreach (blk; parallel(localFluidBlocks,1)) { blk.active = false; }
    integrate_in_time(target_time);
    foreach (blk; parallel(localFluidBlocks,1)) { blk.active = true; }

    wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();
    if (GlobalConfig.is_master_task) {
        writef("WALL-CLOCK (s): %.8f \n", wallClockElapsed);
    }
}

void integrate_solid_in_time_implicit(double target_time, bool init_precondition_matrix)
{
    auto wallClockStart = Clock.currTime();
    double wallClockElapsed;

    // set some time integration parameters
    bool dual_time_stepping = false;
    double physicalSimTime = 0.0;
    int physical_step = 0;
    double target_physical_time, dt_physical;
    int temporal_order, max_physical_steps;
    if (GlobalConfig.sdluOptions.implicitTimeIntegrationMode == 0) {
        // steady-state operation via a backward Euler method
        dual_time_stepping = false;
        temporal_order = 0;
        dt_physical = GlobalConfig.dt_init;
        target_physical_time = dt_physical;
        max_physical_steps = 1; // only perform one nonlinear solve
    } else if (GlobalConfig.sdluOptions.implicitTimeIntegrationMode == 1 ||
               GlobalConfig.sdluOptions.implicitTimeIntegrationMode == 2) {
        // time-accurate operation via a backward-difference formula (BDF)
        dual_time_stepping = true;
        temporal_order = 1; // higher-order BDF schemes are not self-starting
        dt_physical = GlobalConfig.dt_init;
        target_physical_time = target_time;
        max_physical_steps = to!int(ceil(target_physical_time/dt_physical)); // maximum number of expected steps to reach target_physical_time
    } else {
        throw new Error("Invalid implicit_time_integration_mode set in user input script, please select either 0 (for steady-state), 1 (for BDF1), or 2 (for BDF2)");
    }
    int startStep = 0;
    double residual = 0.0;
    double normNew, normRef;
    auto cfl = GlobalConfig.sdluOptions.cfl;
    auto nSteps = GlobalConfig.sdluOptions.maxNewtonIterations;
    auto eta = GlobalConfig.sdluOptions.GMRESSolveTolerance;
    auto sigma = GlobalConfig.sdluOptions.perturbationSize;
    auto tol = GlobalConfig.sdluOptions.NewtonSolveTolerance;

    // fill out all entries in the conserved quantity vector with the initial state
    foreach (sblk; parallel(localSolidBlocks,1)) {
        foreach (i, scell; sblk.cells) {
            foreach (ftl; 0..sblk.myConfig.n_flow_time_levels) {
                scell.e[ftl] = scell.e[0];
            }
        }
    }

    // initialize the precondition matrix
    if (init_precondition_matrix) {
        int spatial_order = 0;
        evalRHS(0.0, 0); // ensure the ghost cells are filled with good data
        foreach (sblk; parallel(localSolidBlocks,1)) {
            sblk.initialize_jacobian(spatial_order, sigma);
        }
    }

    // we calculate the first reference residual here
    evalRHS(0.0, 0);
    foreach (sblk; parallel(localSolidBlocks,1)) {
        foreach (i, scell; sblk.cells) {
            if (temporal_order == 0) {
                sblk.Fe[i] = scell.dedt[0].re;
            } else if (temporal_order == 1) {
                // add BDF1 unsteady term TODO: could we need to handle this better?
                sblk.Fe[i] = scell.dedt[0].re - (1.0/dt_physical)*scell.e[0].re + (1.0/dt_physical)*scell.e[1].re;
            } else { // temporal_order = 2
                // add BDF2 unsteady term TODO: could we need to handle this better?
                sblk.Fe[i] = scell.dedt[0].re - (1.5/dt_physical)*scell.e[0].re + (2.0/dt_physical)*scell.e[1].re - (0.5/dt_physical)*scell.e[2].re;
            }
        }
    }
    mixin(dot_over_blocks("normRef", "Fe", "Fe"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &normRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    normRef = sqrt(normRef);

    // start of main time-stepping loop
    while (physicalSimTime < target_physical_time && physical_step < max_physical_steps) {

        if (GlobalConfig.is_master_task) {
            writeln("TIME: ", physicalSimTime, "    TIME-STEP: ", dt_physical);
        }

        // we need to evaluate time-dependent boundary conditions and source terms
        // at the future state, so we update the time at the start of a step
        physicalSimTime += dt_physical;
        physical_step += 1;

        // determine the allowable timestep
        double dt = determine_dt(cfl);
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        }

        // Begin Newton steps ...
        foreach (step; startStep .. nSteps+1) {
            SimState.step = step;

            // implicit solid update
            rpcGMRES_solve(step, physicalSimTime, dt, eta, sigma, dual_time_stepping, temporal_order, dt_physical, normNew);
            foreach (sblk; parallel(localSolidBlocks,1)) {
                int ftl = to!int(sblk.myConfig.n_flow_time_levels-1);
                foreach (i, scell; sblk.cells) {
                    scell.e[ftl] = scell.e[0] + sblk.de[i];
                    scell.T = updateTemperature(scell.sp, scell.e[ftl]);
                }
            }

            // Put new solid state into e[0] ready for next iteration.
            foreach (sblk; parallel(localSolidBlocks,1)) {
                int ftl = to!int(sblk.myConfig.n_flow_time_levels-1);
                foreach (scell; sblk.cells) {
                    swap(scell.e[0],scell.e[ftl]);
                }
            }

            if (GlobalConfig.is_master_task) {
                writeln("RELATIVE RESIDUAL: ", normNew/normRef);
            }

            if (normNew/normRef < tol) { break; }
        }

        startStep = SimState.step+1;

        // after a BDF1 step we can now switch to BDF2 if requested
        if (physical_step == 1) temporal_order = GlobalConfig.sdluOptions.implicitTimeIntegrationMode;

        // shuffle conserved quantities:
        if (temporal_order == 1) {
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (scell; sblk.cells) {
                    // U_n+1 => U_n
                    scell.e[1] = scell.e[0];
                }
            }
        } else { // temporal_order == 2
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (scell; sblk.cells) {
                    // U_n => U_n-1
                    scell.e[2] = scell.e[1];
                    // U_n+1 => U_n
                    scell.e[1] = scell.e[0];
                }
            }
        }

        // we calculate the new reference residual here
        evalRHS(physicalSimTime, 0);
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (i, scell; sblk.cells) {
                if (temporal_order == 0) {
                    sblk.Fe[i] = scell.dedt[0].re;
                } else if (temporal_order == 1) {
                    // add BDF1 unsteady term TODO: could we need to handle this better?
                    sblk.Fe[i] = scell.dedt[0].re - (1.0/dt_physical)*scell.e[0].re + (1.0/dt_physical)*scell.e[1].re;
                } else { // temporal_order = 2
                    // add BDF2 unsteady term TODO: could we need to handle this better?
                    sblk.Fe[i] = scell.dedt[0].re - (1.5/dt_physical)*scell.e[0].re + (2.0/dt_physical)*scell.e[1].re - (0.5/dt_physical)*scell.e[2].re;
                }
            }
        }
        mixin(dot_over_blocks("normRef", "Fe", "Fe"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &normRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        normRef = sqrt(normRef);

        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();
        if (GlobalConfig.is_master_task) {
            writef("WALL-CLOCK (s): %.8f \n", wallClockElapsed);
        }
    }
}

void allocate_global_solid_workspace()
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

void evalRHS(double sim_time, int ftl)
{
    fnCount++;

    exchange_ghost_cell_solid_boundary_data();

    foreach (sblk; localSolidBlocks) {
        sblk.applyPreSpatialDerivActionAtBndryFaces(sim_time, ftl);
    }

    foreach (sblk; localSolidBlocks) {
        sblk.applyPreSpatialDerivActionAtBndryCells(sim_time, ftl);
    }

    foreach (sblk; parallel(localSolidBlocks, 1)) {
        sblk.averageTemperatures();
        sblk.clearSources();
        sblk.computeSpatialDerivatives(ftl);
    }

    exchange_ghost_cell_solid_boundary_data();
    exchange_ghost_cell_gas_solid_boundary_data();

    foreach (sblk; parallel(localSolidBlocks, 1)) {
        sblk.averageTGradients();
        sblk.computeFluxes();
    }

    foreach (sblk; localSolidBlocks) {
        sblk.applyPostFluxAction(sim_time, ftl);
    }

    foreach (sblk; parallel(localSolidBlocks, 1)) {
        foreach (scell; sblk.cells) {
            if (GlobalConfig.udfSolidSourceTerms) {
                addUDFSourceTermsToSolidCell(sblk.myL, scell, sim_time, sblk);
            }
            scell.timeDerivatives(ftl, GlobalConfig.dimensions);
        }
    }

} // end evalRHS

void evalMatVecProd(double pseudoSimTime, double sigma)
{
    version(complex_numbers) {
        int perturbed_ftl = to!int(GlobalConfig.n_flow_time_levels-1);
        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach (sblk; parallel(localSolidBlocks,1)) {
            sblk.clearSources();
            foreach (i, scell; sblk.cells) {
                scell.e[perturbed_ftl] = scell.e[0];
                scell.e[perturbed_ftl] += complex(0.0, sigma*sblk.zed[i].re);
                scell.T = updateTemperature(scell.sp, scell.e[perturbed_ftl]);
            }
        }

        evalRHS(pseudoSimTime, perturbed_ftl);

        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (i, scell; sblk.cells) {
                sblk.zed[i] = scell.dedt[perturbed_ftl].im/(sigma);
            }
            // we must explicitly remove the imaginary components from the cell and interface flowstates
            foreach(i, scell; sblk.cells) {
                scell.clear_imaginary_components();
            }

        }
    } else {
        throw new Error("Oops. Steady-State Solver setting: useComplexMatVecEval is not compatible with real-number version of the code.");
    }
}

string dot_over_blocks(string dot, string A, string B)
{
    return `
foreach (sblk; parallel(localSolidBlocks,1)) {
   sblk.dotAcc = 0.0;
   foreach (k; 0 .. sblk.nvars) {
      sblk.dotAcc += sblk.`~A~`[k].re*sblk.`~B~`[k].re;
   }
}
`~dot~` = 0.0;
foreach (sblk; localSolidBlocks) `~dot~` += sblk.dotAcc;`;

}

void rpcGMRES_solve(int step, double pseudoSimTime, double dt, double eta, double sigma,
                    bool dual_time_stepping, int temporal_order, double dt_physical, ref double residual)
{

    int maxIters = GlobalConfig.sssOptions.maxOuterIterations;
    // We add 1 because the user thinks of "re"starts, so they
    // might legitimately ask for no restarts. We still have
    // to execute at least once.
    int maxRestarts = GlobalConfig.sssOptions.maxRestarts + 1;
    size_t m = to!size_t(maxIters);
    size_t r;
    size_t iterCount;
    double beta, beta0;
    double outerTol;
    number resid;
    double linSolResid;

    // Variables for max rates of change
    // Use these for equation scaling.
    double minNonDimVal = 1.0; // minimum value used for non-dimensionalisation
                               // when our time rates of change are very small
                               // then we'll avoid non-dimensionalising by
                               // values close to zero.

    // Compute the RHS residual and store dedt[0] as F(e)
    evalRHS(pseudoSimTime, 0);

    foreach (sblk; parallel(localSolidBlocks,1)) {
        foreach (i, scell; sblk.cells) {
            if (temporal_order == 0) {
                sblk.Fe[i] = scell.dedt[0].re;
            } else if (temporal_order == 1) {
                // add BDF1 unsteady term TODO: could we need to handle this better?
                sblk.Fe[i] = scell.dedt[0].re - (1.0/dt_physical)*scell.e[0].re + (1.0/dt_physical)*scell.e[1].re;
            } else { // temporal_order = 2
                // add BDF2 unsteady term TODO: could we need to handle this better?
                sblk.Fe[i] = scell.dedt[0].re - (1.5/dt_physical)*scell.e[0].re + (2.0/dt_physical)*scell.e[1].re - (0.5/dt_physical)*scell.e[2].re;
            }
        }
    }

    // Compute the approximate Jacobian matrix for preconditioning
    foreach (sblk; parallel(localSolidBlocks,1)) {
        sblk.evaluate_jacobian();
        foreach ( ref entry; sblk.jacobian.local.aa) { entry *= -1; }
        foreach (cell; sblk.cells) {
            double dtInv = 1.0/dt;
            if (dual_time_stepping) {
                if (temporal_order == 1) {
                    dtInv = dtInv + 1.0/dt_physical;
                } else {
                    double dtInv_physical  = 1.5/dt_physical;
                    dtInv = dtInv + dtInv_physical;
                }
            }
            sblk.jacobian.local[cell.id,cell.id] = sblk.jacobian.local[cell.id,cell.id] + dtInv;
        }
        nm.smla.decompILU0(sblk.jacobian.local);
    }

    // Determine the max rates in F(e) for scaling the linear system
    foreach (sblk; parallel(localSolidBlocks,1)) {
        sblk.maxRate = 0.0;
        foreach (i, cell; sblk.cells) {
            sblk.maxRate = fmax(sblk.maxRate, fabs(sblk.Fe[i]));
        }
    }
    double maxRate = 0.0;
    foreach (sblk; localSolidBlocks) {
        maxRate = fmax(maxRate, sblk.maxRate);
    }

    // In distributed memory, reduce the max values and ensure everyone has a copy
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(maxRate.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    // Place some guards when time-rate-of-changes are very small
    maxRate = fmax(maxRate, minNonDimVal);

    // Get a copy of the maxRate out to each block
    bool useScaling = true;
    foreach (sblk; parallel(localSolidBlocks,1)) {
        if (useScaling) {
            sblk.maxRate = maxRate;
        }
        else { // just scale by 1
            sblk.maxRate = 1.0;
        }
    }

    // Compute the unscaled L2 norm for reporting the residual of the non-linear system of equations
    double unscaledNorm2;
    mixin(dot_over_blocks("unscaledNorm2", "Fe", "Fe"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &unscaledNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    unscaledNorm2 = sqrt(unscaledNorm2);

    // the remainder of the routine closely follows the structure of rpcGMRES found in nm/smla.d
    // we set the initial guess to zero
    foreach (sblk; parallel(localSolidBlocks,1)) { sblk.x0[] = 0.0; }

    // Start outer-loop of restarted GMRES
    for ( r = 0; r < maxRestarts; r++ ) {

        // Initialise some arrays and matrices that have already been allocated
        g0[] = 0.0;
        g1[] = 0.0;
        H0.zeros();
        H1.zeros();
        Gamma.eye();

        // 1. Evaluate r0 = b - A.x0, beta, v1

        // evaluate A.x0 using a Frechet derivative (note that the zed[] array is baked into the evalJacobianVecProd routine).
        foreach (sblk; parallel(localSolidBlocks,1)) { sblk.zed[] = sblk.x0[]; }

        // Prepare 'w' with (I/dt)(P^-1)v term;
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (i, cell; sblk.cells) {
                double dtInv = 1.0/dt;
                if (dual_time_stepping) {
                    if (temporal_order == 1) {
                        dtInv = dtInv + 1.0/dt_physical;
                    } else {
                        dtInv = dtInv + 1.5/dt_physical;
                    }
                }
                sblk.w[i] = dtInv*sblk.zed[i];
            }
        }

        // Evaluate Jz and place in z
        evalMatVecProd(pseudoSimTime, sigma);

        // Now we can complete calculation of r0
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (k; 0 .. sblk.nvars) sblk.r0[k] = sblk.Fe[k] - (sblk.w[k] - sblk.zed[k]);
        }

        // apply the system scaling to r0
        double eps = 1.0e-50;
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (i, cell; sblk.cells) {
                sblk.r0[i] = (1.0/sblk.maxRate)*sblk.Fe[i];
                // we have observed that sometimes the nonlinear solver can find a solution such that F(e) is exactly 0.0,
                // this causes division by zero throughout the GMRES algorithm, so we add on a very small finite number
                // to r0. This appears to alleviate the issue without degrading convergence.
                sblk.r0[i] += eps;
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
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (k; 0 .. sblk.nvars) {
                sblk.v[k] = sblk.r0[k]/beta;
                sblk.V[k,0] = sblk.v[k];
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
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (i, cell; sblk.cells) {
                    sblk.v[i] *= (sblk.maxRate);
                }
            }

            // apply preconditioning here
            //foreach (sblk; parallel(localSolidBlocks,1)) {
            //    sblk.zed[] = sblk.v[];
            //}
            foreach (sblk; parallel(localSolidBlocks,1)) {
                sblk.zed[] = sblk.v[];
                nm.smla.solve(sblk.jacobian.local, sblk.zed);
            }

            // Prepare 'w' with (I/dt)(P^-1)v term;
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (i, cell; sblk.cells) {
                    double dtInv = 1.0/(dt);
                    if (dual_time_stepping) {
                        if (temporal_order == 1) {
                            dtInv = dtInv + 1.0/dt_physical;
                        } else {
                            dtInv = dtInv + 1.5/dt_physical;
                        }
                    }
                    sblk.w[i] = dtInv*sblk.zed[i];
                }
            }

            // Evaluate Jz and place in z
            evalMatVecProd(pseudoSimTime, sigma);

            // Now we can complete calculation of w
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (k; 0 .. sblk.nvars)  sblk.w[k] = sblk.w[k] - sblk.zed[k];
            }

            // apply the linear system scaling
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (i, cell; sblk.cells) {
                    sblk.w[i] *= (1.0/sblk.maxRate);
                }
            }

            // The remainder of the algorithm looks a lot like any standard GMRES implementation (for example, see smla.d)
            foreach (i; 0 .. j+1) {
                foreach (sblk; parallel(localSolidBlocks,1)) {
                    // Extract column 'i'
                    foreach (k; 0 .. sblk.nvars) sblk.v[k] = sblk.V[k,i];
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
                foreach (sblk; parallel(localSolidBlocks,1)) {
                    foreach (k; 0 .. sblk.nvars) sblk.w[k] -= H0_ij*sblk.v[k];
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

            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (k; 0 .. sblk.nvars) {
                    sblk.v[k] = sblk.w[k]/H0_jp1j;
                    sblk.V[k,j+1] = sblk.v[k];
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

            // Prepare for next step
            copy(H1, H0);
            g0[] = g1[];
            copy(Q1, Q0);

            // Get residual
            resid = fabs(g1[j+1]);
            //nIters = to!int(iterCount);
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
        foreach (sblk; localSolidBlocks) sblk.g1[] = g1[];
        foreach (sblk; parallel(localSolidBlocks,1)) {
            nm.bbla.dot!double(sblk.V, sblk.nvars, m, sblk.g1, sblk.zed);
        }

        // Undo the linear system scaling to recover the unscaled solution vector
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (i, cell; sblk.cells) {
                sblk.zed[i] *= (sblk.maxRate);
            }
        }

        // Apply preconditioning step
        //foreach(sblk; parallel(localSolidBlocks,1)) {
        //    sblk.de[] = sblk.zed[];
        //}
        foreach(sblk; parallel(localSolidBlocks,1)) {
            sblk.de[] = sblk.zed[];
            nm.smla.solve(sblk.jacobian.local, sblk.de);
        }

        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (k; 0 .. sblk.nvars) sblk.de[k] += sblk.x0[k];
        }

        if ( resid <= outerTol || r+1 == maxRestarts ) {
            // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
            // DEBUG:  writefln("Breaking restart loop.");
            // DEBUG:  writefln("RANK %d: breaking restart loop, resid= %e, r+1= %d", GlobalConfig.mpi_rank_for_local_task, resid, r+1);
            break;
        }

        // Else, prepare for a restart by setting the inital
        // guess to the current best estimate of the solution
        foreach (sblk; parallel(localSolidBlocks,1)) { sblk.x0[] = sblk.de[]; }

    }

    if (GlobalConfig.is_master_task) {
        // TODO: handle this output in a more production-ready way
        auto outFile = File("e4-nk.solid.diagnostics.dat", "a");
        outFile.writef("%d    %.16e    %d    %d    %.16e    %.16e \n", step, dt, r, iterCount, linSolResid, unscaledNorm2);
        outFile.close();
    }
    residual = unscaledNorm2;
}

void verify_jacobian() {
    // This routine performs a verification of the numerical Jacobian implementation.
    // It does this by first constructing a numerical Jacobian and a test vector and
    // then comparing the output of an explicit matrix-vector product of the Jacobian
    // and test vector to the output of a Frechet derivative using the test vector.

    // check we are only operating with a single solid block (we expect 1 FluidBlock and 1 SolidBlock)
    if (GlobalConfig.nBlocks > 2) {
        throw new Error("ERROR: the solid domain Jacobian verification routine currently only operates on a single solid block.");
    }
    SSolidBlock sblk = localSolidBlocks[0];

    // perturbation parameter used in numerical Jacobian construction and Frechet derivative evaluation
    double eps = 1.0e-50;

    // we perform a residual evaluation to ensure the ghost cells are filled with good data
    evalRHS(0.0, 0);

    // calculate the numerical Jacobian
    int spatial_order = 2;
    sblk.initialize_jacobian(spatial_order, eps);
    // write out some intermediate information about the cell and face stencils
    foreach (cell; sblk.cells) {
        writef("cell: %d    cell_list: %d    face_list: %d \n", cell.id, cell.cell_list.length, cell.face_list.length);
    }
    sblk.evaluate_jacobian();

    // create an arbitrary unit vector
    double[] vec;
    vec.length = sblk.cells.length;
    foreach ( i, ref val; vec) { val = i+1; }

    // normalise the vector
    double norm = 0.0;
    foreach( i; 0..vec.length) { norm += vec[i]*vec[i]; }
    norm = sqrt(norm);
    foreach( ref val; vec) { val = val/norm; }

    // prepare the result vectors
    double[] sol1;
    sol1.length = vec.length;
    double[] sol2;
    sol2.length = vec.length;

    // explicit multiplication of J*vec
    nm.smla.multiply(sblk.jacobian.local, vec, sol1);

    // Frechet derivative of J*vec
    evalRHS(0.0, 0);
    evalJacobianVecProd(eps, vec, sol2);

    // write out results
    string fileName = "jacobian_test.output";
    auto outFile = File(fileName, "w");
    foreach( i; 0..vec.length ) {
        outFile.writef("%d    %.16e    %.16e    %.16f    %.16f \n", i, fabs((sol1[i]-sol2[i])/sol1[i]).re, fabs(sol1[i]-sol2[i]).re, sol1[i].re, sol2[i].re);
    }

    // stop the program at this point
    import core.stdc.stdlib : exit;
    exit(0);
} // end verify_jacobian

void evalJacobianVecProd(double eps, double[] vec, ref double[] sol) {
    SSolidBlock sblk = localSolidBlocks[0];
    version(complex_numbers) {
        // We perform a Frechet derivative to evaluate J*D^(-1)v
        sblk.clearSources();
        foreach (i, cell; sblk.cells) {
            cell.e[2] = cell.e[0];
            cell.e[2] += complex(0.0, eps*vec[i].re);
            cell.T = updateTemperature(cell.sp, cell.e[2]);
        }
        evalRHS(0.0, 2);
        foreach (i, cell; sblk.cells) {
            sol[i] = cell.dedt[2].im/eps;
        }
        // we must explicitly remove the imaginary components from the cell and interface flowstates
        foreach(i, cell; sblk.cells) {
            cell.clear_imaginary_components();
        }
    } else {
        throw new Error("ERROR: the Jacobian verification test is currently only implemented for the complex-number code path.");
    }
}
