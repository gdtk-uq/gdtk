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
import nm.complex;
import nm.number;

import globalconfig;
import globaldata;
import solid_full_face_copy;
import solid_gas_full_face_copy;
import simcore_exchange;
import solid_udf_source_terms;
import solidprops;

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

@nogc
double determine_dt(double cfl_value)
{
    double dt_local, dt;
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

void solid_update(int step, double pseudoSimTime, double cfl, double eta, double sigma)
{
    // determine the allowable timestep
    double dt = determine_dt(cfl);
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }

    // transfer data between solid and fluid domains before performing the solid domain update
    exchange_ghost_cell_gas_solid_boundary_data();

    // implicit solid update
    rpcGMRES_solve(step, pseudoSimTime, dt, eta, sigma);
    foreach (sblk; parallel(localSolidBlocks,1)) {
        foreach (i, scell; sblk.cells) {
            scell.e[1] = scell.e[0] + sblk.de[i];
            scell.T = updateTemperature(scell.sp, scell.e[1]);
        }
    }

    // Put new solid state into e[0] ready for next iteration.
    foreach (sblk; parallel(localSolidBlocks,1)) {
        foreach (scell; sblk.cells) {
            scell.e[0] = scell.e[1];
        }
    }

    // transfer data between solid and fluid domains after updating the solid domain
    exchange_ghost_cell_gas_solid_boundary_data();
    return;
}

void allocate_global_solid_workspace()
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


void evalRHS(double pseudoSimTime, int ftl)
{
    fnCount++;

    exchange_ghost_cell_solid_boundary_data();

    foreach (sblk; localSolidBlocks) {
        sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl);
    }

    foreach (sblk; localSolidBlocks) {
        sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl);
    }

    foreach (sblk; parallel(localSolidBlocks, 1)) {
        sblk.averageTemperatures();
        sblk.clearSources();
        sblk.computeSpatialDerivatives(ftl);
    }

    exchange_ghost_cell_solid_boundary_data();

    foreach (sblk; parallel(localSolidBlocks, 1)) {
        sblk.computeFluxes();
    }

    foreach (sblk; localSolidBlocks) {
        sblk.applyPostFluxAction(SimState.time, ftl);
    }

    foreach (sblk; parallel(localSolidBlocks, 1)) {
        foreach (scell; sblk.activeCells) {
            if (GlobalConfig.udfSolidSourceTerms) {
                addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
            }
            scell.timeDerivatives(ftl, GlobalConfig.dimensions);
        }
    }

} // end evalRHS


void evalMatVecProd(double pseudoSimTime, double sigma)
{
    version(complex_numbers) {

        // We perform a Frechet derivative to evaluate J*D^(-1)v
        foreach (sblk; parallel(localSolidBlocks,1)) {
            sblk.clearSources();
            foreach (i, scell; sblk.cells) {
                scell.e[1] = scell.e[0];
                scell.e[1] += complex(0.0, sigma*sblk.zed[i].re);
                scell.T = updateTemperature(scell.sp, scell.e[1]);
            }
        }

        evalRHS(pseudoSimTime, 1);

        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (i, scell; sblk.cells) {
                sblk.zed[i] = scell.dedt[1].im/(sigma);
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

void rpcGMRES_solve(int step, double pseudoSimTime, double dt, double eta, double sigma)
{
    number resid;
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

    // Store dedt[0] as F(e)
    foreach (sblk; parallel(localSolidBlocks,1)) {
        int cellCount = 0;
        sblk.maxRate = 0.0;
        foreach (i, scell; sblk.cells) {
            sblk.Fe[cellCount] = scell.dedt[0];
            cellCount += 1;
            sblk.maxRate = fmax(sblk.maxRate, fabs(scell.dedt[0]));
        }
    }

    number maxRate = 0.0;
    foreach (sblk; localSolidBlocks) {
        maxRate = fmax(maxRate, sblk.maxRate);
    }

    // In distributed memory, reduce the max values and ensure everyone has a copy
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(maxRate.re), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    // Place some guards when time-rate-of-changes are very small.
    maxRate = fmax(maxRate, minNonDimVal);

    // Get a copy of the maxes out to each block
    bool useScaling = true;
    foreach (sblk; parallel(localSolidBlocks,1)) {
        if (useScaling) {
            sblk.maxRate = maxRate;
        }
        else { // just scale by 1
            sblk.maxRate = 1.0;
        }
    }

    double unscaledNorm2;
    mixin(dot_over_blocks("unscaledNorm2", "Fe", "Fe"));
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &unscaledNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    unscaledNorm2 = sqrt(unscaledNorm2);

    // Initialise some arrays and matrices that have already been allocated
    g0[] = to!number(0.0);
    g1[] = to!number(0.0);
    H0.zeros();
    H1.zeros();

    // We'll scale r0 against these max rates of change.
    // r0 = b - A*x0
    // Taking x0 = [0] (as is common) gives r0 = b = FU
    // apply scaling
    foreach (sblk; parallel(localSolidBlocks,1)) {
        sblk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (scell; sblk.cells) {
            sblk.r0[cellCount] = (1./sblk.maxRate)*sblk.Fe[cellCount];
            cellCount += 1;
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
    g0[0] = beta;
    foreach (sblk; parallel(localSolidBlocks,1)) {
        foreach (k; 0 .. sblk.nvars) {
            sblk.v[k] = sblk.r0[k]/beta;
            sblk.V[k,0] = sblk.v[k];
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
            foreach (sblk; parallel(localSolidBlocks,1)) {
                int cellCount = 0;
                foreach (cell; sblk.cells) {
                    sblk.v[cellCount] *= (sblk.maxRate);
                    cellCount += 1;
                }
            }

            // TODO: apply preconditioning here
            foreach (sblk; parallel(localSolidBlocks,1)) {
                sblk.zed[] = sblk.v[];
            }

            // Prepare 'w' with (I/dt)(P^-1)v term;
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (i, cell; sblk.cells) {
                    foreach (k; 0..1) {
                        ulong idx = i*1 + k;
                        number dtInv;
                        dtInv = 1.0/(dt);
                        sblk.w[idx] = dtInv*sblk.zed[idx];
                    }
                }
            }

            // Evaluate Jz and place in z
            evalMatVecProd(pseudoSimTime, sigma);

            // Now we can complete calculation of w
            foreach (sblk; parallel(localSolidBlocks,1)) {
                foreach (k; 0 .. sblk.nvars)  sblk.w[k] = sblk.w[k] - sblk.zed[k];
            }

            // apply scaling
            foreach (sblk; parallel(localSolidBlocks,1)) {
                int cellCount = 0;
                foreach (cell; sblk.cells) {
                    sblk.w[cellCount] *= (1./sblk.maxRate);
                    cellCount += 1;
                }
            }

            // The remainder of the algorithm looks a lot like any standard
            // GMRES implementation (for example, see smla.d)
            foreach (i; 0 .. j+1) {
                foreach (sblk; parallel(localSolidBlocks,1)) {
                    // Extract column 'i'
                    foreach (k; 0 .. sblk.nvars ) sblk.v[k] = sblk.V[k,i];
                }
                number H0_ij;
                mixin(dot_over_blocks("H0_ij", "w", "v"));
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(H0_ij.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
                }
                H0[i,j] = H0_ij;
                foreach (sblk; parallel(localSolidBlocks,1)) {
                    foreach (k; 0 .. sblk.nvars) sblk.w[k] -= H0_ij*sblk.v[k];
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
        foreach (sblk; localSolidBlocks) sblk.g1[] = g1[];
        foreach (sblk; parallel(localSolidBlocks,1)) {
            nm.bbla.dot!number(sblk.V, sblk.nvars, m, sblk.g1, sblk.zed);
        }

        // apply scaling
        foreach (sblk; parallel(localSolidBlocks,1)) {
            int cellCount = 0;
            foreach (cell; sblk.cells) {
                sblk.zed[cellCount] *= (sblk.maxRate);
                cellCount += 1;
            }
        }

        // TODO: apply preconditioning here
        foreach(sblk; parallel(localSolidBlocks,1)) {
            sblk.de[] = sblk.zed[];
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

        // Else, we prepare for restart by setting x0 and computing r0
        // Computation of r0 as per Fraysee etal (2005)
        foreach (sblk; parallel(localSolidBlocks,1)) {
            sblk.x0[] = sblk.de[];
        }

        foreach (sblk; localSolidBlocks) copy(Q1, sblk.Q1);
        // Set all values in g0 to 0.0 except for final (m+1) value
        foreach (i; 0 .. m) g0[i] = 0.0;
        foreach (sblk; localSolidBlocks) sblk.g0[] = g0[];
        foreach (sblk; parallel(localSolidBlocks,1)) {
            nm.bbla.dot(sblk.Q1, m, m+1, sblk.g0, sblk.g1);
        }
        foreach (sblk; parallel(localSolidBlocks,1)) {
            nm.bbla.dot(sblk.V, sblk.nvars, m+1, sblk.g1, sblk.r0);
        }

        // apply scaling
        foreach (sblk; parallel(localSolidBlocks,1)) {
            int cellCount = 0;
            foreach (cell; sblk.cells) {
                sblk.r0[cellCount] *= (1.0/sblk.maxRate);
                cellCount += 1;
            }
        }

        mixin(dot_over_blocks("beta", "r0", "r0"));
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &(beta.re), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            version(complex_numbers) { MPI_Allreduce(MPI_IN_PLACE, &(beta.im), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); }
        }
        beta = sqrt(beta);

        // DEBUG: writefln("OUTER: ON RESTART beta= %e", beta);
        foreach (sblk; parallel(localSolidBlocks,1)) {
            foreach (k; 0 .. sblk.nvars) {
                sblk.v[k] = sblk.r0[k]/beta;
                sblk.V[k,0] = sblk.v[k];
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

    // TODO: handle this output in a more production-ready way
    auto outFile = File("e4-nk.solid.diagnostics.dat", "a");
    outFile.writef("%d %.16e %.16e \n", step, dt, unscaledNorm2);
    outFile.close();
}

