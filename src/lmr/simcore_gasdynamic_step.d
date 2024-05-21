// simcore_gasdynamic_step.d
// The core gasdynamic increment functions extracted from simcore.
// 2021-04-14
//

module simcore_gasdynamic_step;

import std.math;
import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.parallelism;
import ntypes.complex;
import nm.number;
import nm.bbla;
import nm.schedule;

import geom;
import geom.misc.kdtree;
import gas;
import conservedquantities;
import globalconfig;
import globaldata;
import flowstate;
import fluidblock;
import sfluidblock;
import ufluidblock;
import ssolidblock;
import solidfvinterface;
import solid_full_face_copy;
import solid_gas_full_face_copy;
import bc.ghost_cell_effect.gas_solid_full_face_copy;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import grid_motion;
import grid_motion_udf;
import grid_motion_shock_fitting;
version(mpi_parallel) {
    import mpi;
}
version(FSI) { import fsi; }
import simcore_exchange;


// To avoid race conditions, there are a couple of locations where
// each block will put its result into the following arrays,
// then we will reduce across the arrays.
shared static double[] local_dt_allow;
shared static double[] local_dt_allow_parab;
shared static double[] local_cfl_max;
shared static int[] local_invalid_cell_count;


void determine_time_step_size()
{
    // Set the size of the time step to be the minimum allowed for any active block.
    // We will check it occasionally, if we have not elected to keep fixed time steps.
    bool do_dt_check_now =
        ((SimState.step % GlobalConfig.cfl_count) == 0) ||
        (SimState.dt_override > 0.0);
    version(mpi_parallel) {
        // If one task is doing a time-step check, all tasks have to.
        int myFlag = to!int(do_dt_check_now);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        do_dt_check_now = to!bool(myFlag);
    }
    if (do_dt_check_now) {
        // Adjust the time step...
        // First, check what each block thinks should be the allowable step size.
        // Also, if we have done some steps, check the CFL limits.
        double cfl_value = GlobalConfig.cfl_schedule.interpolate_value(SimState.time);
        cfl_value *= GlobalConfig.cfl_scale_factor;
        foreach (i, myblk; parallel(localFluidBlocksBySize,1)) {
            // Note 'i' is not necessarily the block id but
            // that is not important here, just need a unique spot to poke into local_dt_allow.
            if (myblk.active) {
                double[3] results = myblk.determine_time_step_size(SimState.dt_global, cfl_value, (SimState.step > 0));
                local_dt_allow[i] = results[0];
                local_cfl_max[i] = results[1];
                local_dt_allow_parab[i] = results[2];
            }
        }
        // Second, reduce this estimate across all local blocks.
        // The starting values are sure to be replaced.
        SimState.dt_allow = double.max;
	SimState.dt_allow_parab = double.max;
        SimState.cfl_max = 0.0;
        foreach (i, myblk; localFluidBlocks) { // serial loop
            if (myblk.active) {
                SimState.dt_allow = min(SimState.dt_allow, local_dt_allow[i]);
                SimState.dt_allow_parab = min(SimState.dt_allow_parab, local_dt_allow_parab[i]);
                SimState.cfl_max = max(SimState.cfl_max, local_cfl_max[i]);
            }
        }
        version(mpi_parallel) {
            double my_dt_allow = SimState.dt_allow;
            double my_dt_allow_parab = SimState.dt_allow_parab;
            double my_cfl_max = SimState.cfl_max;
            MPI_Allreduce(MPI_IN_PLACE, &my_dt_allow, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &my_dt_allow_parab, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &my_cfl_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            SimState.dt_allow = my_dt_allow;
            SimState.dt_allow_parab = my_dt_allow_parab;
            SimState.cfl_max = my_cfl_max;
        }
        if (GlobalConfig.with_super_time_stepping) {
            if (SimState.step == 0) {
                // When starting out, we may override the computed value.
                // This might be handy for situations where the computed estimate
                // is likely to be not small enough for numerical stability.
                SimState.dt_allow = fmin(GlobalConfig.dt_init, SimState.dt_allow);
                SimState.dt_allow_parab = fmin(GlobalConfig.dt_init, SimState.dt_allow_parab);
            }
            // Now, change the actual time step, as needed.
            if (SimState.dt_allow <= SimState.dt_global) {
                // If we need to reduce the time step, do it immediately.
                SimState.dt_global = SimState.dt_allow;
                SimState.dt_global_parab = SimState.dt_allow_parab;
            } else {
                // Make the transitions to larger time steps gentle.
                SimState.dt_global = min(SimState.dt_global*1.5, SimState.dt_allow);
                SimState.dt_global_parab = min(SimState.dt_global_parab*1.5, SimState.dt_allow_parab);
                // The user may supply, explicitly, a maximum time-step size.
                SimState.dt_global = min(SimState.dt_global, GlobalConfig.dt_max);
            }
        } else if (GlobalConfig.with_local_time_stepping) {
            SimState.dt_global = SimState.dt_allow;
        } else { // do some global time-stepping checks
            if (SimState.step == 0) {
                // When starting out, we may override the computed value.
                // This might be handy for situations where the computed estimate
                // is likely to be not small enough for numerical stability.
                SimState.dt_allow = fmin(GlobalConfig.dt_init, SimState.dt_allow);
            }
            if (SimState.dt_override > 0.0) {
                // The user-defined supervisory function atTimestepStart may have set
                // dt_override because it knows something about the simulation conditions
                // that is not handled well by our generic CFL check.
                SimState.dt_allow = fmin(SimState.dt_override, SimState.dt_allow);
            }
            // Now, change the actual time step, as needed.
            if (SimState.dt_allow <= SimState.dt_global) {
                // If we need to reduce the time step, do it immediately.
                SimState.dt_global = SimState.dt_allow;
            } else {
                // Make the transitions to larger time steps gentle.
                SimState.dt_global = min(SimState.dt_global*1.5, SimState.dt_allow);
                // The user may supply, explicitly, a maximum time-step size.
                SimState.dt_global = min(SimState.dt_global, GlobalConfig.dt_max);
            }
        }
    } // end if do_cfl_check_now
} // end determine_time_step_size()


void sts_gasdynamic_explicit_increment_with_fixed_grid()
{
    if (GlobalConfig.reacting && GlobalConfig.chemistry_update == ChemistryUpdateMode.integral) {
        throw new Error("Not implemented: explicit gasdynamic update with integral chemistry update.");
    }
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    //
    // compute number of super steps
    bool euler_step = false;
    double dt_global;
    double alpha;
    double s_RKL;
    int S;

    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.steady_fluid_transient_solid) {
        if (GlobalConfig.fixed_time_step) {
            // if we have a fixed time step then we need to specify the number of super-steps
            S = SimState.s_RKL;
            if (S == 1) { euler_step = true; }
            dt_global = GlobalConfig.dt_init;
        } else {
            // We determine the allowable timestep here to overwrite the timestep calculated for the fluid domain
            // TODO: think about a more appropriate place to calculate the timestep. KAD 2022-11-08
            double cfl_value = GlobalConfig.cfl_schedule.interpolate_value(SimState.time);
            dt_global = double.max;
            bool first = true;
            foreach (i, sblk; parallel(localSolidBlocks, 1)) {
                double dt_local = sblk.determine_time_step_size(cfl_value);
                if (first) {
                    local_dt_allow[i] = dt_local;
                    first = false;
                } else {
                    local_dt_allow[i] = fmin(local_dt_allow[i], dt_local);
                }
            }
            foreach (i,sblk; localSolidBlocks) {
                dt_global = fmin(dt_global, local_dt_allow[i]);
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }

            // set the super-time-step
            SimState.dt_global_parab = dt_global;
            S = SimState.s_RKL;
            double dt_super;
            if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                dt_super = dt_global * (S*S+S)/(2.0); // RKL1
            } else {
                dt_super = dt_global * (S*S+S-2.0)/(4.0); // RKL2
            }

            // check for a couple of edge cases when using the user specified S is not appropriate...
            double dt_remainder = SimState.target_time - SimState.time;
            if (dt_global > dt_remainder) {
                // the explicit time-step is larger than the remaining simulation time,
                // so we just set the super-time-step to the remaining time and perform an Euler step
                S = 1;
                dt_super = dt_remainder;
            }
            else if (dt_remainder < dt_super) {
                // the remaining time is larger than the explicit time-step but smaller than the super-time-step,
                // we borrow from the mixed equation formulation to find a suitable S
                alpha = dt_remainder/dt_global;
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                    s_RKL = 0.5*(-1.0+sqrt(1+8.0*alpha)); // RKL1
                } else {
                    s_RKL = 0.5*(-1.0+sqrt(9+16.0*alpha));  // RKL2
                }
                S = to!int(floor(s_RKL));
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                    dt_super = dt_global * (S*S+S)/(2.0); // RKL1
                } else {
                    dt_super = dt_global * (S*S+S-2.0)/(4.0); // RKL2
                }
            }

            // place the super-time-step into the dt_global variables
            if (S == 1) { euler_step = true; }
            dt_global = dt_super;
            SimState.dt_global = dt_global;
        }
    } else {
        if (GlobalConfig.fixed_time_step) {

            // if we have a fixed time step then we won't have calculated either a hyperbolic or parabolic dt
            // we need to specify the number of super-steps in this case
            S = 7;
            dt_global = GlobalConfig.dt_init;

        } else {

            // otherwise we will calculate the suitable number of super-steps as per:
            //    A stabilized Runge–Kutta–Legendre method for explicit super-time-stepping of parabolic and mixed equations
            //    Journal of Computational Physics,
            //    Meyer et al. (2014)

            double dt_hyp = SimState.dt_global;
            double dt_parab = SimState.dt_global_parab;
            alpha = (dt_hyp)/(dt_parab);

            if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                s_RKL = 0.5*(-1.0+sqrt(1+8.0*alpha)); // RKL1
            } else {
                s_RKL = 0.5*(-1.0+sqrt(9+16.0*alpha));  // RKL2
            }

            // it is recommended to round S down to the nearest odd integer for stability
            S = to!int(floor(s_RKL));
            if ( fmod(S, 2) == 0 && s_RKL != 1 ) { S = S - 1; }

            // When dt_parab is approxmately equal to or greater than dt_hyper (i.e. S <= 1), then we will just do a simple Euler step
            if (S <= 1) {
                S = 1;
                euler_step = true;
                if ((SimState.step % GlobalConfig.print_count) == 0 && GlobalConfig.is_master_task) {
                    writeln("WARNING: dtPARAB (", SimState.dt_global_parab, ") ~= dtHYPER (", SimState.dt_global, ") .... taking an Euler step.");
                }
            }

            // Due to rounding S we should alter the time-step to be consistent
            if (euler_step) {
                dt_global = SimState.dt_global;
            } else {
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                    dt_global = SimState.dt_global_parab * (S*S+S)/(2.0); // RKL1
                } else {
                    dt_global = SimState.dt_global_parab * (S*S+S-2.0)/(4.0); // RKL2
                }
            }
        }
    }
    SimState.s_RKL = S;
    // --------------------------------------------------
    // j = 1
    // --------------------------------------------------

    // First-stage of gas-dynamic update.
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level remains at zero for the non-moving grid

    // Preparation for the inviscid gas-dynamic flow update.
    if (GlobalConfig.udf_source_terms) {
        foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
            if (!blk.active) continue;
            int local_ftl = ftl;
            int local_gtl = gtl;
            double local_sim_time = SimState.time;
            foreach (cell; blk.cells) {
                size_t i_cell = cell.id; size_t j_cell = 0; size_t k_cell = 0;
                if (blk.grid_type == Grid_t.structured_grid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                    i_cell = ijk_indices[0]; j_cell = ijk_indices[1]; k_cell = ijk_indices[2];
                }
                getUDFSourceTermsForCell(blk.myL, cell, local_gtl, local_sim_time,
                                         blk.myConfig, blk.id, i_cell, j_cell, k_cell);
            }
        }
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) {
	    blk.clear_fluxes_of_conserved_quantities();
	    foreach (cell; blk.cells) {
		cell.clear_source_vector();
		cell.data_is_bad = false;
	    }
	}
    }
    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
    exchange_ghost_cell_gas_solid_boundary_data();
    if (GlobalConfig.apply_bcs_in_parallel) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	}
    } else {
	foreach (blk; localFluidBlocksBySize) {
	    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	}
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) { blk.set_face_flowstates_to_averages_from_cells(); }
    }
    // And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
    foreach (sblk; parallel(localSolidBlocks, 1)) {
	if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
    }
    // We've put this detector step here because it needs the ghost-cell data
    // to be current, as it should be just after a call to apply_convective_bc().
    if (GlobalConfig.do_shock_detect) detect_shocks(gtl, ftl);

    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
    }

    // for unstructured blocks we need to transfer the convective gradients before the flux calc
    if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
        exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
    }

    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) { blk.convective_flux_phase2(allow_high_order_interpolation, gtl); }
    }

    if (GlobalConfig.apply_bcs_in_parallel) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	}
    } else {
	foreach (blk; localFluidBlocksBySize) {
	    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	}
    }
    if (GlobalConfig.viscous) {
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) {
                    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                }
            }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.flow_property_spatial_derivatives(gtl);
            }
        }
        // for unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
        // we need to transfer the viscous gradients before the flux calc
        exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, gtl, ftl);
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                // we need to average cell-centered spatial (/viscous) gradients to get approximations of the gradients
                // at the cell interfaces before the viscous flux calculation.
                if (blk.myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                    foreach(f; blk.faces) {
                        f.average_cell_deriv_values(0);
                    }
                }
            }
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            if (blk.active) {
                blk.estimate_turbulence_viscosity();
            }
        }
        // we exchange boundary data at this point to ensure the
        // ghost cells along block-block boundaries have the most
        // recent mu_t and k_t values.
        exchange_ghost_cell_turbulent_viscosity();
        exchange_ghost_cell_gas_solid_boundary_data();
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.average_turbulent_transprops_to_faces();
                blk.viscous_flux();
            }
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
            }
        }
    } // end if viscous
    foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
	if (!blk.active) continue;
	int blklocal_ftl = ftl;
	int blklocal_gtl = gtl;
	double dt = dt_global;
	double blklocal_sim_time = SimState.time;
	foreach (cell; blk.cells) {
	    cell.add_inviscid_source_vector(blklocal_gtl, blk.omegaz);
	    if (blk.myConfig.viscous) { cell.add_viscous_source_vector(); }
	    if (blk.myConfig.udf_source_terms) { cell.add_udf_source_vector(); }
	    cell.time_derivatives(blklocal_gtl, blklocal_ftl);
	}
	if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(blklocal_ftl); }
        auto cqi = blk.myConfig.cqi;
        // Coefficients for the j=1 stage update.
        double s = SimState.s_RKL;
        double muj_tilde;
        if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 || euler_step) {
            // RKL1 (j=1)
            muj_tilde = 2.0/(s*s+s);
        } else {
            // RKL2 (j=1)
            muj_tilde = 4.0/(3.0*(s*s+s-2.0));
        }
        // The update itself.
	foreach (cell; blk.cells) {
            auto U0 = cell.U[0]; auto U1 = cell.U[1]; auto dUdt0 = cell.dUdt[0];
            foreach (k; 0 .. cqi.n) {
                U1[k] = U0[k] + muj_tilde*dt*dUdt0[k];
            }
	    if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 || euler_step) {
                // RKL1 (j=1)
                // Nothing more to do.
            } else {
                // RKL2 (j=1)
                // Make a copy of the initial conserved quantities.
                // It will be used in the j > 1 stages.
                cell.U[3].copy_values_from(cell.U[0]);
            }
            cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
	}
        try {
            local_invalid_cell_count[i] = blk.count_invalid_cells(blklocal_gtl, blklocal_ftl+1);
        }  catch (Exception e) {
            debug {
                writefln("Exception thrown in count_invalid_cells: %s", e.msg);
                writefln("  mpi_rank=%d", GlobalConfig.mpi_rank_for_local_task);
            }
            local_invalid_cell_count[i] = to!int(blk.cells.length);
        }
    } // end foreach blk
    //
    int flagTooManyBadCells = 0;
    foreach (i, blk; localFluidBlocksBySize) { // serial loop
	if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
	    flagTooManyBadCells = 1;
	    writefln("Following first-stage gasdynamic update: %d bad cells in block[%d].",
		     local_invalid_cell_count[i], i);
	}
    }
    version(mpi_parallel) {
	MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (flagTooManyBadCells > 0) {
	throw new FlowSolverException("Too many bad cells; go home.");
    }
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ||
        GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.steady_fluid_transient_solid) {
	// Next do solid domain update IMMEDIATELY after at same flow time level
        exchange_ghost_cell_solid_boundary_data(); // we need up to date temperatures for the spatial derivative calc.
	foreach (sblk; parallel(localSolidBlocks, 1)) {
	    if (!sblk.active) continue;
	    sblk.averageTemperatures();
	    sblk.averageProperties();
	    sblk.clearSources();
	    sblk.computeSpatialDerivatives(ftl);
        }
        exchange_ghost_cell_solid_boundary_data(); // we need to transfer the spatial derivatives for the flux eval.
        exchange_ghost_cell_gas_solid_boundary_data();
        foreach (sblk; parallel(localSolidBlocks, 1)) {
            if (!sblk.active) continue;
            sblk.averageTGradients();
            sblk.computeFluxes();
	}
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (sblk; parallel(localSolidBlocks, 1)) {
		if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
	    }
	} else {
	    foreach (sblk; localSolidBlocks) {
		    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
	    }
	}
	// We need to synchronise before updating
	foreach (sblk; parallel(localSolidBlocks, 1)) {
	    foreach (scell; sblk.cells) {
		if (GlobalConfig.udfSolidSourceTerms) {
		    addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
		}
		scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                if (euler_step) {
                    scell.eulerUpdate(dt_global);
                } else {
                    if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                        scell.stage1RKL1Update(dt_global, 1, SimState.s_RKL); // RKL1 (j=1)
                    } else {
                        scell.stage1RKL2Update(dt_global, 1, SimState.s_RKL); // RKL2 (j=1)
                    }
                }
                scell.ss.e = scell.e[ftl+1];
                sblk.stm.updateTemperature(scell.ss);
                scell.T = scell.ss.T;
	    } // end foreach scell
	} // end foreach sblk
    } // end if tight solid domain coupling.

    if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 || euler_step) {
        SimState.time = t0 + (2.0/(S*S+S))*dt_global; // RKL1
    } else {
        SimState.time = t0 + (4.0/(3.0*(S*S+S-2.0)))*dt_global; // RKL2
    }
    // --------------------------------------------------
    // 2 <= j <= S
    // --------------------------------------------------
    //writeln("max: ", SimState.s_RKL+1);
    foreach (j; 2..SimState.s_RKL+1) {
        if (euler_step) { continue;}
        ftl = 1;

        // Preparation for the inviscid gas-dynamic flow update.
        if (GlobalConfig.udf_source_terms && GlobalConfig.eval_udf_source_terms_at_each_stage) {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                int blklocal_ftl = ftl;
                int blklocal_gtl = gtl;
                double blklocal_sim_time = SimState.time;
                foreach (cell; blk.cells) {
                    size_t i_cell = cell.id; size_t j_cell = 0; size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                        i_cell = ijk_indices[0]; j_cell = ijk_indices[1]; k_cell = ijk_indices[2];
                    }
                    getUDFSourceTermsForCell(blk.myL, cell, blklocal_gtl, blklocal_sim_time,
                                             blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                }
            }
        }
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {
		blk.clear_fluxes_of_conserved_quantities();
		foreach (cell; blk.cells) {
		    cell.clear_source_vector();
		    cell.data_is_bad = false;
		}
	    }
	}
        exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
        exchange_ghost_cell_gas_solid_boundary_data();
        if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	    }
	}
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.set_face_flowstates_to_averages_from_cells(); }
        }
	// And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
	foreach (sblk; parallel(localSolidBlocks, 1)) {
	    if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
	}
	// We've put this detector step here because it needs the ghost-cell data
	// to be current, as it should be just after a call to apply_convective_bc().
        if (GlobalConfig.do_shock_detect) detect_shocks(gtl, ftl);

	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
	}

        // for unstructured blocks we need to transfer the convective gradients before the flux calc
        if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
            exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
        }

	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
	}

        // for unstructured blocks we need to transfer the convective gradients before the flux calc
        if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
            exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
        }

	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.convective_flux_phase2(allow_high_order_interpolation, gtl); }
	}

	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	    }
	}
        if (GlobalConfig.viscous) {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            }
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.flow_property_spatial_derivatives(gtl);
                }
            }
            // for unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
            // we need to transfer the viscous gradients before the flux calc
            exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, gtl, ftl);
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    // we need to average cell-centered spatial (/viscous) gradients to get approximations of the gradients
                    // at the cell interfaces before the viscous flux calculation.
                    if (blk.myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                        foreach(f; blk.faces) {
                            f.average_cell_deriv_values(0);
                        }
                    }
                }
            }
            foreach (blk; parallel(localFluidBlocks,1)) {
                if (blk.active) {
                    blk.estimate_turbulence_viscosity();
                }
            }
            // we exchange boundary data at this point to ensure the
            // ghost cells along block-block boundaries have the most
            // recent mu_t and k_t values.
            exchange_ghost_cell_turbulent_viscosity();
            exchange_ghost_cell_gas_solid_boundary_data();
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.average_turbulent_transprops_to_faces();
                    blk.viscous_flux();
                }
            }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            }
        } // end if viscous
	foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
	    if (!blk.active) continue;
	    int blklocal_ftl = ftl;
	    int blklocal_gtl = gtl;
	    double dt = dt_global;
	    foreach (cell; blk.cells) {
		cell.add_inviscid_source_vector(blklocal_gtl, blk.omegaz);
		if (blk.myConfig.viscous) { cell.add_viscous_source_vector(); }
		if (blk.myConfig.udf_source_terms) { cell.add_udf_source_vector(); }
		cell.time_derivatives(blklocal_gtl, blklocal_ftl);
	    }
	    if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(blklocal_ftl); }
            //
            // Coefficients for the update, stage j.
            double s = SimState.s_RKL;
            double ajm1; double bj; double bjm1; double bjm2;
            double muj; double vuj; double muj_tilde; double gam_tilde;
            if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                // RKL1
                muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);
                muj = (2.0*j-1)/j;
                vuj = (1.0-j)/j;
            } else {
                // RKL2
                if (j == 2) {
                    bj = 1.0/3.0;
                    bjm1 = 1.0/3.0;
                    bjm2 = 1.0/3.0;
                } else if (j == 3) {
                    bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
                    bjm1 = 1.0/3.0;
                    bjm2 = 1.0/3.0;
                } else if (j == 4) {
                    bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
                    bjm1 = ((j-1.0)*(j-1.0)+(j-1.0)-2.0)/(2.0*(j-1.0)*((j-1.0)+1.0));
                    bjm2 = 1.0/3.0;
                } else {
                    bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
                    bjm1 = ((j-1.0)*(j-1.0)+(j-1.0)-2.0)/(2.0*(j-1.0)*((j-1.0)+1.0));
                    bjm2 = ((j-2.0)*(j-2.0)+(j-2.0)-2.0)/(2.0*(j-2.0)*((j-2.0)+1.0));
                }
                ajm1 = 1.0-bjm1;
                muj_tilde = (4.0*(2.0*j-1.0))/(j*(s*s+s-2.0)) * (bj/bjm1);
                gam_tilde = -ajm1*muj_tilde;
                muj = (2*j-1.0)/(j) * (bj/bjm1);
                vuj = -(j-1.0)/(j) * (bj/bjm2);
            }
            //
            // The update itself.
            auto cqi = blk.myConfig.cqi;
	    foreach (cell; blk.cells) {
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                    // RKL1
                    auto U0 = cell.U[0]; auto U1 = cell.U[1]; auto U2 = cell.U[2];
                    auto dUdt0 = cell.dUdt[1];
                    foreach (k; 0 .. cqi.n) {
                        U2[k] = muj*U1[k] + vuj*U0[k] + muj_tilde*dt*dUdt0[k];
                    }
                } else {
                    // RKL2
                    auto U0 = cell.U[0]; auto U1 = cell.U[1]; auto U2 = cell.U[2]; auto U3 = cell.U[3];
                    auto dUdt0 = cell.dUdt[1]; auto dUdtO = cell.dUdt[0]; // Note the subtly-different names.
                    foreach (k; 0 .. cqi.n) {
                        U2[k] = muj*U1[k] + vuj*U0[k] + (1.0-muj-vuj)*U3[k] +
                            muj_tilde*dt*dUdt0[k] + gam_tilde*dt*dUdtO[k];
                    }
                }
		cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
	    }
            try {
                local_invalid_cell_count[i] = blk.count_invalid_cells(blklocal_gtl, blklocal_ftl+1);
            }  catch (Exception e) {
                debug {
                    writefln("Exception thrown in count_invalid_cells: %s", e.msg);
                    writefln("  mpi_rank=%d", GlobalConfig.mpi_rank_for_local_task);
                }
                local_invalid_cell_count[i] = to!int(blk.cells.length);
            }
	} // end foreach blk
	//
	flagTooManyBadCells = 0;
	foreach (i, blk; localFluidBlocksBySize) { // serial loop
	    if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
		flagTooManyBadCells = 1;
		writefln("Following first-stage gasdynamic update: %d bad cells in block[%d].",
			 local_invalid_cell_count[i], i);
	    }
	}
	version(mpi_parallel) {
	    MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	}
	if (flagTooManyBadCells > 0) {
	    throw new FlowSolverException("Too many bad cells; go home.");
	}
	//
	if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ||
            GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.steady_fluid_transient_solid) {
        // Next do solid domain update IMMEDIATELY after at same flow time level
        exchange_ghost_cell_solid_boundary_data(); // we need up to date temperatures for the spatial derivative calc.
        foreach (sblk; parallel(localSolidBlocks, 1)) {
            if (!sblk.active) continue;
            sblk.averageTemperatures();
            sblk.averageProperties();
            sblk.clearSources();
            sblk.computeSpatialDerivatives(ftl);
        }
        exchange_ghost_cell_solid_boundary_data(); // we need to transfer the spatial derivatives for the flux eval.
        foreach (sblk; parallel(localSolidBlocks, 1)) {
            if (!sblk.active) continue;
            sblk.averageTGradients();
            sblk.computeFluxes();
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (sblk; parallel(localSolidBlocks, 1)) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        } else {
            foreach (sblk; localSolidBlocks) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        }
        // We need to synchronise before updating
        foreach (sblk; parallel(localSolidBlocks, 1)) {
            foreach (scell; sblk.cells) {
                if (GlobalConfig.udfSolidSourceTerms) {
                    addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
                }
                scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                    scell.stage2RKL1Update(dt_global, j, SimState.s_RKL); // RKL1 (j=1)
                } else {
                    scell.stage2RKL2Update(dt_global, j, SimState.s_RKL); // RKL2 (j=1)
                }
                scell.ss.e = scell.e[ftl+1];
                sblk.stm.updateTemperature(scell.ss);
                scell.T = scell.ss.T;
            } // end foreach scell
        } // end foreach sblk
    } // end if tight solid domain coupling.

	// shuffle time-levels for next iteration (U1 goes to U0 & U2 goes to U1)
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {
		foreach (cell; blk.cells) {
                    cell.U[0].copy_values_from(cell.U[1]);
                    cell.U[1].copy_values_from(cell.U[2]);
                }
	    }
	} // end foreach blk
        if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ||
            GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.steady_fluid_transient_solid) {
            foreach (sblk; parallel(localSolidBlocks,1)) {
                if (sblk.active) {
                    foreach (scell; sblk.cells) {
                        scell.e[0] = scell.e[1];
                        scell.e[1] = scell.e[2];
                    }
                }
            } // end foreach blk
        }
        if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
            SimState.time = t0 + ((j*j+j)/(S*S+S))*dt_global; // RKL1
        } else {
            SimState.time = t0 + ((j*j+j-2.0)/(S*S+S-2.0))*dt_global; // RKL2
        }
    } // end foreach (J; 2..S+1)

    size_t end_indx;
    if (euler_step) { end_indx = 1; }
    else { end_indx = 2; }

    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            foreach (cell; blk.cells) { cell.U[0].copy_values_from(cell.U[end_indx]); }
        }
    } // end foreach blk
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ||
        GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.steady_fluid_transient_solid) {
        foreach (sblk; parallel(localSolidBlocks, 1)) {
            if (sblk.active) {
                foreach (scell; sblk.cells) { scell.e[0] = scell.e[end_indx]; }
            }
        } // end foreach sblk
    }
    // Finally, update the globally known simulation time for the whole step.
    SimState.time = t0 + dt_global;
} // end sts_gasdynamic_explicit_increment_with_fixed_grid()


void gasdynamic_explicit_increment_with_fixed_grid()
{
    if (GlobalConfig.reacting && GlobalConfig.chemistry_update == ChemistryUpdateMode.integral) {
        throw new Error("Not implemented: explicit gasdynamic update with integral chemistry update.");
    }
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    //
    // Set the time-step coefficients for the stages of the update scheme.
    // These coefficients are used to increment the time variable for each stage of update.
    // Later, each stage will set its own coefficients for the dependent variables.
    shared double c2 = 1.0; // default for predictor-corrector update
    shared double c3 = 1.0; // default for predictor-corrector update
    shared double c4 = 1.0; // Only used for classic_rk4
    final switch ( GlobalConfig.gasdynamic_update_scheme ) {
    case GasdynamicUpdate.euler:
    case GasdynamicUpdate.backward_euler:
    case GasdynamicUpdate.implicit_rk1:
    case GasdynamicUpdate.pc: c2 = 1.0; c3 = 1.0; break;
    case GasdynamicUpdate.midpoint: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.classic_rk3: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.tvd_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.denman_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.rkl1:
    case GasdynamicUpdate.rkl2: assert(false, "invalid option");
    case GasdynamicUpdate.moving_grid_1_stage:
    case GasdynamicUpdate.moving_grid_2_stage: assert(false, "invalid option");
    case GasdynamicUpdate.classic_rk4: c2 = 0.5; c3 = 0.5; c4 = 1.0;
    }
    //
    int attempt_number = 0;
    int step_failed = 0; // Use int because we want to reduce across MPI ranks.
    int n_stages = to!int(number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme));
    int flagTooManyBadCells = 0;
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            foreach (cell; blk.cells) { cell.data_is_bad = false; }
        }
    }
    do {
        ++attempt_number;
        step_failed = 0;
        foreach (stage; 1 .. n_stages+1) {
            shared int ftl = stage-1; // time-level within the overall convective-update
            shared int gtl = 0; // grid time-level remains at zero for the non-moving grid
            // Could simplify the following by changing c to an array.
            switch (stage) {
            case 1: SimState.time = t0; break;
            case 2: SimState.time = t0 + c2 * SimState.dt_global; break;
            case 3: SimState.time = t0 + c3 * SimState.dt_global; break;
            case 4: SimState.time = t0 + c4 * SimState.dt_global; break;
            default: throw new Error("Invalid state for explicit update.");
            }
            // Phase 01 LOCAL
            if (GlobalConfig.udf_source_terms &&
                ((stage == 1) || GlobalConfig.eval_udf_source_terms_at_each_stage)) {
                foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                    if (!blk.active) continue;
                    int blklocal_ftl = ftl;
                    int blklocal_gtl = gtl;
                    double blklocal_sim_time = SimState.time;
                    foreach (cell; blk.cells) {
                        size_t i_cell = cell.id; size_t j_cell = 0; size_t k_cell = 0;
                        if (blk.grid_type == Grid_t.structured_grid) {
                            auto sblk = cast(SFluidBlock) blk;
                            assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                            auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                            i_cell = ijk_indices[0]; j_cell = ijk_indices[1]; k_cell = ijk_indices[2];
                        }
                        getUDFSourceTermsForCell(blk.myL, cell, blklocal_gtl, blklocal_sim_time,
                                                 blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                    }
                }
            }
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.clear_fluxes_of_conserved_quantities();
                    foreach (cell; blk.cells) { cell.clear_source_vector(); }
                }
            }
            // Phase 02 (maybe) MPI
            // We are relying on exchanging boundary data as a pre-reconstruction activity.
            exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
            exchange_ghost_cell_gas_solid_boundary_data();
            // Phase 03 LOCAL
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                }
            }
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.set_face_flowstates_to_averages_from_cells(); }
            }
            // Phase 04 MPI
            // We've put this detector step here because it needs the ghost-cell data
            // to be current, as it should be just after a call to apply_convective_bc().
            if ((GlobalConfig.do_shock_detect) && (stage == 1) &&
                ((!GlobalConfig.frozen_shock_detector) ||
                 (GlobalConfig.shock_detector_freeze_step > SimState.step))) {
                detect_shocks(gtl, ftl);
            }
            // Phase 05a LOCAL
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
            }
            // Phase 05b MPI
            // for unstructured blocks we need to transfer the convective gradients before the flux calc
            if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
            }
            // Phase 06a LOCAL
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
            }
            // Phase 06b MPI
            // for unstructured blocks we need to transfer the convective gradients before the flux calc
            if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
            }
            // Phase 07 LOCAL
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.convective_flux_phase2(allow_high_order_interpolation, gtl); }
            }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            }
            if (GlobalConfig.viscous) {
                // Phase 08 LOCAL
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
                            blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                            blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                        }
                    }
                } else {
                    foreach (blk; localFluidBlocksBySize) {
                        if (blk.active) {
                            blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                            blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                        }
                    }
                }
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.flow_property_spatial_derivatives(gtl);
                    }
                }
                // Phase 09 MPI
                // For unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
                // we need to transfer the viscous gradients before the flux calc
                exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, gtl, ftl);
                //
                // Phase 10 LOCAL
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        // we need to average cell-centered spatial (/viscous) gradients
                        // to get approximations of the gradients
                        // at the cell interfaces before the viscous flux calculation.
                        if (blk.myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                            foreach(f; blk.faces) {
                                f.average_cell_deriv_values(0);
                            }
                        }
                    }
                }
                foreach (blk; parallel(localFluidBlocks,1)) {
                    if (blk.active) {
                        blk.estimate_turbulence_viscosity();
                    }
                }
                // Phase 11 MPI
                // We exchange boundary data at this point to ensure the
                // ghost cells along block-block boundaries have the most
                // recent mu_t and k_t values.
                exchange_ghost_cell_turbulent_viscosity();
                //
                // Phase 12 LOCAL
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.average_turbulent_transprops_to_faces();
                        blk.viscous_flux();
                    }
                }
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                    }
                } else {
                    foreach (blk; localFluidBlocksBySize) {
                        if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                    }
                }
            } // end if viscous
            //
            // Phase 13 LOCAL
            try {
                foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                    if (!blk.active) continue;
                    auto cqi = blk.myConfig.cqi;
                    int blklocal_ftl = ftl;
                    int blklocal_gtl = gtl;
                    bool blklocal_with_local_time_stepping = with_local_time_stepping;
                    double blklocal_dt_global = SimState.dt_global;
                    foreach (cell; blk.cells) {
                        cell.add_inviscid_source_vector(blklocal_gtl, blk.omegaz);
                        if (blk.myConfig.viscous) { cell.add_viscous_source_vector(); }
                        if (blk.myConfig.udf_source_terms) { cell.add_udf_source_vector(); }
                        cell.time_derivatives(blklocal_gtl, blklocal_ftl);
                    }
                    if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(blklocal_ftl); }
                    switch (stage) {
                    case 1:
                        double gamma_1 = 1.0; // for normal Predictor-Corrector or Euler update.
                        final switch (blk.myConfig.gasdynamic_update_scheme) {
                        case GasdynamicUpdate.euler:
                        case GasdynamicUpdate.backward_euler:
                        case GasdynamicUpdate.implicit_rk1:
                        case GasdynamicUpdate.moving_grid_1_stage:
                        case GasdynamicUpdate.moving_grid_2_stage:
                        case GasdynamicUpdate.pc: gamma_1 = 1.0; break;
                        case GasdynamicUpdate.midpoint: gamma_1 = 0.5; break;
                        case GasdynamicUpdate.classic_rk3: gamma_1 = 0.5; break;
                        case GasdynamicUpdate.tvd_rk3: gamma_1 = 1.0; break;
                        case GasdynamicUpdate.denman_rk3: gamma_1 = 8.0/15.0; break;
                        case GasdynamicUpdate.rkl1:
                        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
                        case GasdynamicUpdate.classic_rk4: gamma_1 = 1.0 / 2.0; break;
                        }
                        foreach (cell; blk.cells) {
                            auto dt = (blklocal_with_local_time_stepping) ? cell.dt_local : blklocal_dt_global;
                            auto dUdt0 = cell.dUdt[0];
                            auto U0 = cell.U[0];
                            auto U1 = cell.U[1];
                            foreach (j; 0 .. cqi.n) {
                                U1[j] = U0[j] + dt*gamma_1*dUdt0[j];
                            }
                            version(turbulence) {
                                foreach(j; 0 .. cqi.n_turb){
                                    U1[cqi.rhoturb+j] = fmax(U1[cqi.rhoturb+j], U0[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                                }
                                // ...assuming a minimum value of 1.0 for omega
                                // It may occur (near steps in the wall) that a large flux of romega
                                // through one of the cell interfaces causes romega within the cell
                                // to drop rapidly.
                                // The large values of omega come from Menter's near-wall correction that may be
                                // applied outside the control of this finite-volume core code.
                                // These large values of omega will be convected along the wall and,
                                // if they are convected past a corner with a strong expansion,
                                // there will be an unreasonably-large flux out of the cell.
                            }
                            version(MHD) {
                                if (blk.myConfig.MHD && blk.myConfig.divergence_cleaning && !blk.myConfig.MHD_static_field) {
                                    U1[cqi.psi] *= cell.divergence_damping_factor(dt, blk.myConfig.c_h, blk.myConfig.divB_damping_length);
                                }
                            }
                            cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
                        } // end foreach cell
                        break;
                    case 2:
                        double gamma_1 = 0.5; // Presume predictor-corrector.
                        double gamma_2 = 0.5;
                        final switch (blk.myConfig.gasdynamic_update_scheme) {
                        case GasdynamicUpdate.euler:
                        case GasdynamicUpdate.backward_euler:
                        case GasdynamicUpdate.implicit_rk1:
                        case GasdynamicUpdate.moving_grid_1_stage:
                            assert(false, "invalid option for multi-stage update.");
                        case GasdynamicUpdate.moving_grid_2_stage:
                        case GasdynamicUpdate.pc: gamma_1 = 0.5, gamma_2 = 0.5; break;
                        case GasdynamicUpdate.midpoint: gamma_1 = 0.0; gamma_2 = 1.0; break;
                        case GasdynamicUpdate.classic_rk3: gamma_1 = -1.0; gamma_2 = 2.0; break;
                        case GasdynamicUpdate.tvd_rk3: gamma_1 = 0.25; gamma_2 = 0.25; break;
                        case GasdynamicUpdate.denman_rk3: gamma_1 = -17.0/60.0; gamma_2 = 5.0/12.0; break;
                        case GasdynamicUpdate.rkl1:
                        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
                        case GasdynamicUpdate.classic_rk4: gamma_1 = 0.0; gamma_2 = 1.0 / 2.0; break;
                        }
                        foreach (cell; blk.cells) {
                            auto dt = (blklocal_with_local_time_stepping) ? cell.dt_local : blklocal_dt_global;
                            auto dUdt0 = cell.dUdt[0];
                            auto dUdt1 = cell.dUdt[1];
                            auto U_old = cell.U[0];
                            if (blk.myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) { U_old = cell.U[1]; }
                            auto U2 = cell.U[2];
                            foreach (j; 0 .. cqi.n) {
                                U2[j] = U_old[j] + dt*(gamma_1*dUdt0[j] + gamma_2*dUdt1[j]);
                            }
                            version(turbulence) {
                                foreach(j; 0 .. cqi.n_turb){
                                    U2[cqi.rhoturb+j] = fmax(U2[cqi.rhoturb+j], U_old[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                                }
                            }
                            version(MHD) {
                                if (blk.myConfig.MHD && blk.myConfig.divergence_cleaning && !blk.myConfig.MHD_static_field) {
                                    U2[cqi.psi] *= cell.divergence_damping_factor(dt, blk.myConfig.c_h, blk.myConfig.divB_damping_length);
                                }
                            }
                            cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
                        } // end foreach cell
                        break;
                    case 3:
                        double gamma_1 = 1.0/6.0; // presume TVD_RK3 scheme.
                        double gamma_2 = 1.0/6.0;
                        double gamma_3 = 4.0/6.0;
                        final switch (blk.myConfig.gasdynamic_update_scheme) {
                        case GasdynamicUpdate.euler:
                        case GasdynamicUpdate.backward_euler:
                        case GasdynamicUpdate.implicit_rk1:
                        case GasdynamicUpdate.moving_grid_1_stage:
                        case GasdynamicUpdate.moving_grid_2_stage:
                        case GasdynamicUpdate.pc:
                        case GasdynamicUpdate.midpoint:
                            assert(false, "invalid option for 3-stage update.");
                        case GasdynamicUpdate.classic_rk3:
                            gamma_1 = 1.0/6.0; gamma_2 = 4.0/6.0; gamma_3 = 1.0/6.0;
                            break;
                        case GasdynamicUpdate.tvd_rk3:
                            gamma_1 = 1.0/6.0; gamma_2 = 1.0/6.0; gamma_3 = 4.0/6.0;
                            break;
                        case GasdynamicUpdate.denman_rk3:
                            // FIX-ME: Check that we have Andrew Denman's scheme ported correctly.
                            gamma_1 = 0.0; gamma_2 = -5.0/12.0; gamma_3 = 3.0/4.0;
                            break;
                        case GasdynamicUpdate.rkl1:
                        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
                        case GasdynamicUpdate.classic_rk4:
                            gamma_1 = 0.0; gamma_2 = 0.0; gamma_3 = 1.0;
                            break;
                        }
                        foreach (cell; blk.cells) {
                            auto dt = (blklocal_with_local_time_stepping) ? cell.dt_local : blklocal_dt_global;
                            auto dUdt0 = cell.dUdt[0];
                            auto dUdt1 = cell.dUdt[1];
                            auto dUdt2 = cell.dUdt[2];
                            auto U_old = cell.U[0];
                            if (blk.myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) { U_old = cell.U[2]; }
                            auto U3 = cell.U[3];
                            foreach (j; 0 .. cqi.n) {
                                U3[j] = U_old[j] + dt * (gamma_1*dUdt0[j] + gamma_2*dUdt1[j] + gamma_3*dUdt2[j]);
                            }
                            version(turbulence) {
                                foreach(j; 0 .. cqi.n_turb){
                                    U3[cqi.rhoturb+j] = fmax(U3[cqi.rhoturb+j], U_old[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                                }
                            }
                            version(MHD) {
                                if (blk.myConfig.MHD && blk.myConfig.divergence_cleaning && !blk.myConfig.MHD_static_field) {
                                    U3[cqi.psi] *= cell.divergence_damping_factor(dt, blk.myConfig.c_h, blk.myConfig.divB_damping_length);
                                }
                            }
                            cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
                        } // end foreach cell
                        break;
                    case 4:
                        double gamma_1 = 1.0/6.0; // presume classic_rk4 scheme.
                        double gamma_2 = 2.0/6.0;
                        double gamma_3 = 2.0/6.0;
                        double gamma_4 = 1.0/6.0;
                        final switch (blk.myConfig.gasdynamic_update_scheme) {
                        case GasdynamicUpdate.euler:
                        case GasdynamicUpdate.backward_euler:
                        case GasdynamicUpdate.implicit_rk1:
                        case GasdynamicUpdate.moving_grid_1_stage:
                        case GasdynamicUpdate.moving_grid_2_stage:
                        case GasdynamicUpdate.pc:
                        case GasdynamicUpdate.midpoint:
                        case GasdynamicUpdate.classic_rk3:
                        case GasdynamicUpdate.tvd_rk3:
                        case GasdynamicUpdate.denman_rk3:
                        case GasdynamicUpdate.rkl1:
                        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
                        case GasdynamicUpdate.classic_rk4:
                            gamma_1 = 1.0 / 6.0; gamma_2 = 1.0 / 3.0; gamma_3 = 1.0 / 3.0; gamma_4 = 1.0 / 6.0;
                            break;
                        }
                        foreach (cell; blk.cells) {
                            auto dt = (blklocal_with_local_time_stepping) ? cell.dt_local : blklocal_dt_global;
                            auto dUdt0 = cell.dUdt[0];
                            auto dUdt1 = cell.dUdt[1];
                            auto dUdt2 = cell.dUdt[2];
                            auto dUdt3 = cell.dUdt[3];
                            auto U_old = cell.U[0];
                            if (blk.myConfig.gasdynamic_update_scheme == GasdynamicUpdate.denman_rk3) { U_old = cell.U[2]; }
                            auto U4 = cell.U[4];
                            foreach (j; 0 .. cqi.n) {
                                U4[j] = U_old[j] + dt * (gamma_1*dUdt0[j] + gamma_2*dUdt1[j] + gamma_3*dUdt2[j] + gamma_4 * dUdt3[j]);
                            }
                            version(turbulence) {
                                foreach(j; 0 .. cqi.n_turb){
                                    U4[cqi.rhoturb+j] = fmax(U4[cqi.rhoturb+j], U_old[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                                }
                            }
                            version(MHD) {
                                if (blk.myConfig.MHD && blk.myConfig.divergence_cleaning && !blk.myConfig.MHD_static_field) {
                                    U4[cqi.psi] *= cell.divergence_damping_factor(dt, blk.myConfig.c_h, blk.myConfig.divB_damping_length);
                                }
                            }
                            cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
                        } // end foreach cell
                        break;
                    default: throw new Error("Invalid state for explicit update.");
                    } // end switch (stage)
                    //
                    local_invalid_cell_count[i] = blk.count_invalid_cells(blklocal_gtl, blklocal_ftl+1);
                } // end foreach blk
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 13 of stage %d of explicit update: %s", stage, e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                break;
            }
            // Phase 14 LOCAL
            flagTooManyBadCells = 0;
            foreach (i, blk; localFluidBlocksBySize) { // serial loop
                if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                    flagTooManyBadCells = 1;
                    writefln("Following stage %d gasdynamic update: %d bad cells in block[%d].",
                             stage, local_invalid_cell_count[i], i);
                }
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (flagTooManyBadCells > 0) {
                throw new FlowSolverException("Too many bad cells during explicit gasdynamic update.");
            }
            //
            if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
                // Phase 15 LOCAL
                // Next do solid domain update IMMEDIATELY after at same flow time level
                // exchange_ghost_cell_gas_solid_boundary_data();
                exchange_ghost_cell_solid_boundary_data();
                try {
                    if (GlobalConfig.apply_bcs_in_parallel) {
                        foreach (sblk; parallel(localSolidBlocks, 1)) {
                            if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
                        }
                        foreach (sblk; parallel(localSolidBlocks, 1)) {
                            if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl); }
                        }
                    } else {
                        foreach (sblk; localSolidBlocks) {
                            if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
                        }
                        foreach (sblk; localSolidBlocks) {
                            if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl); }
                        }
                    }
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (!sblk.active) continue;
                        sblk.averageTemperatures();
                        sblk.averageProperties();
                        sblk.clearSources();
                        sblk.computeSpatialDerivatives(ftl);
                    }
                } catch (Exception e) {
                    debug { writefln("Exception thrown in phase 15 of stage %d of explicit update: %s", stage, e.msg); }
                    step_failed = 1;
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (step_failed) {
                    // Start the step over again with a reduced time step.
                    SimState.dt_global = SimState.dt_global * 0.2;
                    break;
                }
                //
                // Phase 16 MPI
                exchange_ghost_cell_solid_boundary_data();
                //
                // Phase 17 LOCAL
                try {
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (!sblk.active) continue;
                        sblk.averageTGradients();
                        sblk.computeFluxes();
                    }
                    if (GlobalConfig.apply_bcs_in_parallel) {
                        foreach (sblk; parallel(localSolidBlocks, 1)) {
                            if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                        }
                    } else {
                        foreach (sblk; localSolidBlocks) {
                            if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                        }
                    }
                    // We need to synchronise before updating
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        foreach (scell; sblk.cells) {
                            if (GlobalConfig.udfSolidSourceTerms) {
                                addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
                            }
                            scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                            switch (stage) {
                            case 1: scell.stage1Update(SimState.dt_global); break;
                            case 2: scell.stage2Update(SimState.dt_global); break;
                            case 3: scell.stage3Update(SimState.dt_global); break;
                            case 4: scell.stage4Update(SimState.dt_global); break;
                            default: throw new Error("Invalid state for explicit update.");
                            }
                            scell.ss.e = scell.e[ftl+1];
                            sblk.stm.updateTemperature(scell.ss);
                            scell.T = scell.ss.T;
                        } // end foreach scell
                    } // end foreach sblk
                } catch (Exception e) {
                    debug { writefln("Exception thrown in phase 17 of stage %d of explicit update: %s", stage, e.msg); }
                    step_failed = 1;
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (step_failed) {
                    // Start the step over again with a reduced time step.
                    SimState.dt_global = SimState.dt_global * 0.2;
                    break;
                }
            } // end if tight solid domain coupling.
        } // end foreach stage
    } while (step_failed && (attempt_number < GlobalConfig.max_attempts_for_step));
    //
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        string msg = format("Explicit update failed after %d attempts; giving up.",
                            GlobalConfig.max_attempts_for_step);
        throw new FlowSolverException(msg);
    }
    //
    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (cell; blk.cells) { swap(cell.U[0], cell.U[end_indx]); }
        }
    } // end foreach blk
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
        foreach (sblk; localSolidBlocks) {
            if (sblk.active) {
                size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
                foreach (scell; sblk.cells) { scell.e[0] = scell.e[end_indx]; }
            }
        } // end foreach sblk
    }
    //
    // Finally, update the globally know simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_explicit_increment_with_fixed_grid()


void gasdynamic_explicit_increment_with_moving_grid()
{
    // For moving grid simulations we move the grid on the first predictor step and then
    // leave it fixed in this position for the corrector steps.
    if (GlobalConfig.reacting && GlobalConfig.chemistry_update == ChemistryUpdateMode.integral) {
        throw new Error("Not implemented: explicit gasdynamic update with integral chemistry update.");
    }
    if (GlobalConfig.with_local_time_stepping) {
        throw new Error("Not implemented: cell-local time step with moving grid.");
    }
    shared double t0 = SimState.time;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    // Set the time-step coefficients for the stages of the update scheme.
    shared double c2 = 1.0; // same for 1-stage or 2-stage update
    shared double c3 = 1.0; // ditto
    //
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level
    int step_failed = 0; // Use int because we want to reduce across MPI ranks.
    // Phase 00a LOCAL
    if (GlobalConfig.udf_source_terms) {
        try {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                int blklocal_ftl = ftl;
                int blklocal_gtl = gtl;
                double blklocal_sim_time = SimState.time;
                foreach (cell; blk.cells) {
                    size_t i_cell = cell.id; size_t j_cell = 0; size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                        i_cell = ijk_indices[0]; j_cell = ijk_indices[1]; k_cell = ijk_indices[2];
                    }
                    getUDFSourceTermsForCell(blk.myL, cell, blklocal_gtl, blklocal_sim_time,
                                             blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown while preparing user source terms"~
                             " for explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            throw new FlowSolverException("User-defined-source-term failure.");
        }
    }
    // Phase 00b (mabe) MPI
    final switch(GlobalConfig.grid_motion) {
    case GridMotion.none:
        throw new Error("Should not be setting grid velocities in with GridMotion.none");
    case GridMotion.user_defined:
        // Rely on user to set vertex velocities.
        // Note that velocities remain unchanged if the user does nothing.
        assign_vertex_velocities_via_udf(SimState.time, SimState.dt_global);
        break;
    case GridMotion.shock_fitting:
        if (SimState.time > GlobalConfig.shock_fitting_delay) {
            foreach (i, fba; fluidBlockArrays) {
                if (fba.shock_fitting) { compute_vtx_velocities_for_sf(fba); }
            }
        }
        break;
    version(FSI) {
        case GridMotion.FSI:
            foreach (FEMModel; FEMModels) {
                if (SimState.step % FEMModel.myConfig.couplingStep == 0) {
                    FEMModel.compute_vtx_velocities_for_FSI(SimState.dt_global);
                }
            }
            break;
        }
    } // end switch grid_motion
    //
    int flagTooManyBadCells;
    int attempt_number = 0;
    do {
        ++attempt_number;
        step_failed = 0;
        // Preparation for the predictor-stage of inviscid gas-dynamic flow update.
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.clear_fluxes_of_conserved_quantities();
                foreach (cell; blk.cells) {
                    cell.clear_source_vector();
                    cell.data_is_bad = false;
                }
            }
        }
        ftl = 0; gtl = 0;
        // Phase 01 LOCAL
        try {
            // Moving Grid - predict new vertex positions for moving grid
            foreach (blk; localFluidBlocksBySize) {
                if (!blk.active) continue;
                auto sblk = cast(SFluidBlock) blk;
                assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                // move vertices
                predict_vertex_positions(sblk, SimState.dt_global, gtl);
                // recalculate cell geometry with new vertex positions @ gtl = 1
                blk.compute_primary_cell_geometric_data(gtl+1);
                blk.compute_least_squares_setup(gtl+1);
                // determine interface velocities using GCL for gtl = 1
                set_gcl_interface_properties(sblk, gtl+1, SimState.dt_global);
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 01 of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        gtl = 1; // update gtl now that grid has moved
        // Phase 02 (maybe) MPI
        exchange_ghost_cell_geometry_data();
        exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
        // Phase 03 LOCAL
        try {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                }
            }
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.set_face_flowstates_to_averages_from_cells(); }
            }
            // And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
            foreach (sblk; localSolidBlocks) {
                if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
            }
            foreach (sblk; localSolidBlocks) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 03 of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 04 (maybe) MPI
        // We've put this detector step here because it needs the ghost-cell data
        // to be current, as it should be just after a call to apply_convective_bc().
        if (GlobalConfig.do_shock_detect) {
            detect_shocks(gtl, ftl);
        }
        // Phase 05a LOCAL
        try {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 05a of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 05b (maybe) MPI
        // for unstructured blocks we need to transfer the convective gradients before the flux calc
        if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
            exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
        }
        // Phase 06a LOCAL
        try {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 06a of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 06b (maybe) MPI
        // for unstructured blocks we need to transfer the convective gradients before the flux calc
        if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
            exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
        }
        // Phase 07 LOCAL
        try {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.convective_flux_phase2(allow_high_order_interpolation, gtl); }
            }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 07 of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        if (GlobalConfig.viscous) {
            // Phase 08 LOCAL
            try {
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
                            blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                            blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                        }
                    }
                } else {
                    foreach (blk; localFluidBlocksBySize) {
                        if (blk.active) {
                            blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                            blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                        }
                    }
                }
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.flow_property_spatial_derivatives(gtl);
                    }
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 08 of stage 1 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 09 (maybe) MPI
            // for unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
            // we need to transfer the viscous gradients before the flux calc
            exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, gtl, ftl);
            // Phase 10 LOCAL
            try {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        // we need to average cell-centered spatial (/viscous) gradients
                        // to get approximations of the gradients at the cell interfaces
                        // before the viscous flux calculation.
                        if (blk.myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                            foreach(f; blk.faces) { f.average_cell_deriv_values(0); }
                        }
                    }
                }
                foreach (blk; parallel(localFluidBlocks,1)) {
                    if (blk.active) { blk.estimate_turbulence_viscosity(); }
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 10 of stage 1 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 11 (maybe) MPI
            // we exchange boundary data at this point to ensure the
            // ghost cells along block-block boundaries have the most
            // recent mu_t and k_t values.
            exchange_ghost_cell_turbulent_viscosity();
            // Phase 12 LOCAL
            try {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.average_turbulent_transprops_to_faces();
                        blk.viscous_flux();
                    }
                }
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                    }
                } else {
                    foreach (blk; localFluidBlocksBySize) {
                        if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                    }
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 12 of stage 1 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if viscous
        // Phase 13 LOCAL
        try {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                int blklocal_ftl = ftl;
                int blklocal_gtl = gtl;
                double blklocal_dt_global = SimState.dt_global;
                auto cqi = blk.myConfig.cqi;
                foreach (cell; blk.cells) {
                    cell.add_inviscid_source_vector(blklocal_gtl, blk.omegaz);
                    if (blk.myConfig.viscous) { cell.add_viscous_source_vector(); }
                    if (blk.myConfig.udf_source_terms) { cell.add_udf_source_vector(); }
                    cell.time_derivatives(blklocal_gtl, blklocal_ftl);
                    //
                    auto dt = blklocal_dt_global;
                    auto dUdt0 = cell.dUdt[0];
                    auto U0 = cell.U[0];
                    auto U1 = cell.U[1];
                    number vr = cell.volume[0] / cell.volume[1];
                    foreach (j; 0 .. cqi.n) {
                        U1[j] = vr*(U0[j] + dt*dUdt0[j]);
                    }
                    version(turbulence) {
                        foreach(j; 0 .. cqi.n_turb){
                            U1[cqi.rhoturb+j] = fmax(U1[cqi.rhoturb+j], U1[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                        }
                    }
                    version(MHD) {
                        if (blk.myConfig.MHD && blk.myConfig.divergence_cleaning && !blk.myConfig.MHD_static_field) {
                            U1[cqi.psi] *= cell.divergence_damping_factor(dt, blk.myConfig.c_h, blk.myConfig.divB_damping_length);
                        }
                    }
                    cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
                } // end foreach cell
                local_invalid_cell_count[i] = blk.count_invalid_cells(blklocal_gtl, blklocal_ftl+1);
            } // end foreach blk
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 13 of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 14
        flagTooManyBadCells = 0;
        foreach (i, blk; localFluidBlocksBySize) { // serial loop
            if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                flagTooManyBadCells = 1;
                writefln("Following first-stage gasdynamic update with moving grid: %d bad cells in block[%d].",
                         local_invalid_cell_count[i], i);
            }
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (flagTooManyBadCells > 0) {
            throw new FlowSolverException("Too many bad cells following first-stage gasdynamic update with moving grid.");
        }
        //
        // Next do solid domain update IMMEDIATELY after at same flow time level
        exchange_ghost_cell_solid_boundary_data();
        // Phase 15 LOCAL
        try {
            foreach (sblk; localSolidBlocks) {
                if (!sblk.active) continue;
                sblk.averageTemperatures();
                sblk.averageProperties();
                sblk.clearSources();
                sblk.computeSpatialDerivatives(ftl);
                sblk.applyPostFluxAction(SimState.time, ftl);
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 15 of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 16 (maybe) MPI
        exchange_ghost_cell_solid_boundary_data();
        // Phase 17 LOCAL
        try {
            foreach (sblk; parallel(localSolidBlocks, 1)) {
                if (!sblk.active) continue;
                sblk.averageTGradients();
                sblk.computeFluxes();
                sblk.applyPostFluxAction(SimState.time, ftl);
                foreach (scell; sblk.cells) {
                    if (GlobalConfig.udfSolidSourceTerms) {
                        addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
                    }
                    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                    scell.stage1Update(SimState.dt_global);
                    scell.ss.e = scell.e[ftl+1];
                    sblk.stm.updateTemperature(scell.ss);
                    scell.T = scell.ss.T;
                } // end foreach scell
            } // end foreach sblk
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 17 of stage 1 of explicit update on moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        /////
        if (number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) == 2) {
            // Preparation for second-stage of gas-dynamic update.
            SimState.time = t0 + c2 * SimState.dt_global;
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.clear_fluxes_of_conserved_quantities();
                    foreach (cell; blk.cells) { cell.clear_source_vector(); }
                }
            }
            // Phase 01 LOCAL
            try {
                // Second stage of gas-dynamic update.
                // Moving Grid - update geometry to gtl 2
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        // move vertices - this is a formality since pos[2] = pos[1]
                        predict_vertex_positions(sblk, SimState.dt_global, gtl);
                        // recalculate cell geometry with new vertex positions
                        blk.compute_primary_cell_geometric_data(gtl+1);
                        blk.compute_least_squares_setup(gtl+1);
                        // grid remains at pos[gtl=1], thus let's use old interface velocities
                        // thus no need to set_gcl_interface_properties(blk, 2, dt_global);
                    }
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 01 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            ftl = 1;
            gtl = 2;
            // Phase 02 (maybe) MPI
            // We are relying on exchanging boundary data as a pre-reconstruction activity.
            exchange_ghost_cell_geometry_data();
            exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
            // Phase 03 LOCAL
            try {
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                    }
                } else {
                    foreach (blk; localFluidBlocksBySize) {
                        if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                    }
                }
                // Let's set up solid domain bc's also before changing any flow properties.
                foreach (sblk; localSolidBlocks) {
                    if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
                }
                foreach (sblk; localSolidBlocks) {
                    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                }
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, 0); } // note 0 rather then gtl
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 03 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 04 (maybe) MPI
            // for unstructured blocks we need to transfer the convective gradients before the flux calc
            if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
            }
            // Phase 04 LOCAL
            try {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, 0); } // note 0 rather then gtl
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 04 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 04 (maybe) MPI
            // for unstructured blocks we need to transfer the convective gradients before the flux calc
            if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
            }
            // Phase 05 LOCAL
            try {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.convective_flux_phase2(allow_high_order_interpolation, 0); } // note 0 rather then gtl
                }
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                    }
                } else {
                    foreach (blk; localFluidBlocksBySize) {
                        if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                    }
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 05 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            if (GlobalConfig.viscous) {
                // Phase 06 LOCAL
                try {
                    if (GlobalConfig.apply_bcs_in_parallel) {
                        foreach (blk; parallel(localFluidBlocksBySize,1)) {
                            if (blk.active) {
                                blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                                blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                            }
                        }
                    } else {
                        foreach (blk; localFluidBlocksBySize) {
                            if (blk.active) {
                                blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                                blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                            }
                        }
                    }
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
                            blk.flow_property_spatial_derivatives(gtl);
                        }
                    }
                } catch (Exception e) {
                    debug { writefln("Exception thrown in phase 06 of stage 2 of explicit update on moving grid: %s", e.msg); }
                    step_failed = 1;
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (step_failed) {
                    // Start the step over again with a reduced time step.
                    SimState.dt_global = SimState.dt_global * 0.2;
                    continue;
                }
                // Phase 07 (maybe) MPI
                // for unstructured blocks employing the cell-centered spatial (/viscous) gradient method,
                // we need to transfer the viscous gradients before the flux calc
                exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, gtl, ftl);
                // Phase 08 LOCAL
                try {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
                            // we need to average cell-centered spatial (/viscous) gradients
                            // to get approximations of the gradients at the cell interfaces
                            // before the viscous flux calculation.
                            if (blk.myConfig.spatial_deriv_locn == SpatialDerivLocn.cells) {
                                foreach(f; blk.faces) {
                                    f.average_cell_deriv_values(0);
                                }
                            }
                        }
                    }
                    foreach (blk; parallel(localFluidBlocks,1)) {
                        if (blk.active) {
                            blk.estimate_turbulence_viscosity();
                        }
                    }
                } catch (Exception e) {
                    debug { writefln("Exception thrown in phase 08 of stage 2 of explicit update on moving grid: %s", e.msg); }
                    step_failed = 1;
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (step_failed) {
                    // Start the step over again with a reduced time step.
                    SimState.dt_global = SimState.dt_global * 0.2;
                    continue;
                }
                // Phase 09 (maybe) MPI
                // we exchange boundary data at this point to ensure the
                // ghost cells along block-block boundaries have the most
                // recent mu_t and k_t values.
                exchange_ghost_cell_turbulent_viscosity();
                // Phase 10 LOCAL
                try {
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
                            blk.average_turbulent_transprops_to_faces();
                            blk.viscous_flux();
                        }
                    }
                    if (GlobalConfig.apply_bcs_in_parallel) {
                        foreach (blk; parallel(localFluidBlocksBySize,1)) {
                            if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                        }
                    } else {
                        foreach (blk; localFluidBlocksBySize) {
                            if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                        }
                    }
                } catch (Exception e) {
                    debug { writefln("Exception thrown in phase 10 of stage 2 of explicit update on moving grid: %s", e.msg); }
                    step_failed = 1;
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (step_failed) {
                    // Start the step over again with a reduced time step.
                    SimState.dt_global = SimState.dt_global * 0.2;
                    continue;
                }
            } // end if viscous
            // Phase 11 LOCAL
            try {
                foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                    if (!blk.active) continue;
                    int blklocal_ftl = ftl;
                    int blklocal_gtl = gtl;
                    double blklocal_dt_global = SimState.dt_global;
                    auto cqi = blk.myConfig.cqi;
                    foreach (cell; blk.cells) {
                        cell.add_inviscid_source_vector(blklocal_gtl, blk.omegaz);
                        if (blk.myConfig.viscous) { cell.add_viscous_source_vector(); }
                        if (blk.myConfig.udf_source_terms) { cell.add_udf_source_vector(); }
                        cell.time_derivatives(blklocal_gtl, blklocal_ftl);
                        //
                        auto dt = blklocal_dt_global;
                        auto dUdt0 = cell.dUdt[0];
                        auto dUdt1 = cell.dUdt[1];
                        auto U0 = cell.U[0];
                        auto U2 = cell.U[2];
                        //
                        number v_old = cell.volume[0];
                        number vol_inv = 1.0 / cell.volume[2];
                        number gamma_1 = 0.5 * cell.volume[0]; // Roll-in the volumes for convenience below.
                        number gamma_2 = 0.5 * cell.volume[1];
                        //
                        foreach (j; 0 .. cqi.n) {
                            U2[j] = vol_inv * (v_old * U0[j] + dt * (gamma_1 * dUdt0[j] + gamma_2 * dUdt1[j]));
                        }
                        version(turbulence) {
                            foreach(j; 0 .. cqi.n_turb){
                                U2[cqi.rhoturb+j] = fmax(U2[cqi.rhoturb+j], U2[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                            }
                        }
                        version(MHD) {
                            if (blk.myConfig.MHD && blk.myConfig.divergence_cleaning && !blk.myConfig.MHD_static_field) {
                                U2[cqi.psi] *= cell.divergence_damping_factor(dt, blk.myConfig.c_h, blk.myConfig.divB_damping_length);
                            }
                        }
                        cell.decode_conserved(blklocal_gtl, blklocal_ftl+1, blk.omegaz);
                    } // end foreach cell
                    local_invalid_cell_count[i] = blk.count_invalid_cells(blklocal_gtl, blklocal_ftl+1);
                } // end foreach blk
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 11 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 12
            flagTooManyBadCells = 0;
            foreach (i, blk; localFluidBlocksBySize) { // serial loop
                if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                    flagTooManyBadCells = 1;
                    writefln("Following second-stage gasdynamic update with moving grid: %d bad cells in block[%d].",
                             local_invalid_cell_count[i], i);
                }
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (flagTooManyBadCells > 0) {
                throw new FlowSolverException("Too many bad cells following second-stage gasdynamic update with moving grid.");
            }
            //
            // Do solid domain update IMMEDIATELY after at same flow time level
            exchange_ghost_cell_solid_boundary_data();
            // Phase 13 LOCAL
            try {
                foreach (sblk; localSolidBlocks) {
                    if (!sblk.active) continue;
                    sblk.averageTemperatures();
                    sblk.averageProperties();
                    sblk.clearSources();
                    sblk.computeSpatialDerivatives(ftl);
                    sblk.applyPostFluxAction(SimState.time, ftl);
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 13 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 14 (maybe) MPI
            exchange_ghost_cell_solid_boundary_data();
            // Phase 15 LOCAL
            try {
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
                    sblk.averageTGradients();
                    sblk.computeFluxes();
                    sblk.applyPostFluxAction(SimState.time, ftl);
                    foreach (scell; sblk.cells) {
                        if (GlobalConfig.udfSolidSourceTerms) {
                            addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
                        }
                        scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                        scell.stage2Update(SimState.dt_global);
                        scell.ss.e = scell.e[ftl+1];
                        sblk.stm.updateTemperature(scell.ss);
                        scell.T = scell.ss.T;
                    } // end foreach cell
                } // end foreach blk
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 15 of stage 2 of explicit update on moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if number_of_stages_for_update_scheme >= 2
    } while (step_failed && (attempt_number < GlobalConfig.max_attempts_for_step));
    //
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        string msg = format("Explicit update failed after %d attempts; giving up.", GlobalConfig.max_attempts_for_step);
        throw new FlowSolverException(msg);
    }
    //
    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (cell; blk.cells) { swap(cell.U[0], cell.U[end_indx]); }
        }
    }
    foreach (sblk; localSolidBlocks) {
        if (sblk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (scell; sblk.cells) { scell.e[0] = scell.e[end_indx]; }
        }
    }
    //
    // Update the latest grid level to the new step grid level 0 and recalculate geometry.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            foreach (cell; blk.cells) { cell.copy_grid_level_to_level(gtl, 0); }
            blk.compute_primary_cell_geometric_data(0);
            blk.compute_least_squares_setup(0);
        }
    }
    // Finally, update the globally known simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_explicit_increment_with_moving_grid()


void gasdynamic_implicit_increment_with_fixed_grid()
{
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    //
    immutable int ftl0 = 0;
    immutable int ftl1 = 1;
    immutable int gtl0 = 0; // grid time-level remains at zero for the non-moving grid
    int step_failed = 0; // Use int because we want to reduce across MPI ranks.
    //
    if (GlobalConfig.udf_source_terms) {
        // Phase 00 LOCAL
        try {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                double blklocal_sim_time = SimState.time;
                foreach (cell; blk.cells) {
                    size_t i_cell = cell.id; size_t j_cell = 0; size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                        i_cell = ijk_indices[0]; j_cell = ijk_indices[1]; k_cell = ijk_indices[2];
                    }
                    getUDFSourceTermsForCell(blk.myL, cell, gtl0, blklocal_sim_time,
                                             blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in Phase 00 of implicit update: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            throw new FlowSolverException("Implicit update failed when setting source terms. Giving up.");
        }
    }
    //
    double reaction_fraction = GlobalConfig.reaction_fraction_schedule.interpolate_value(SimState.time);
    int attempt_number = 0;
    int flagTooManyBadCells = 0;
    do {
        ++attempt_number;
        step_failed = 0;
        // Phase 01 LOCAL Preparation for the gas-dynamic flow update.
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.clear_fluxes_of_conserved_quantities();
                foreach (cell; blk.cells) {
                    cell.clear_source_vector();
                    cell.data_is_bad = false;
                }
            }
        }
        // Phase 02 (maybe) MPI
        exchange_ghost_cell_boundary_data(SimState.time, gtl0, ftl0);
        exchange_ghost_cell_gas_solid_boundary_data();
        if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
            exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl0, ftl0);
        }
        if (GlobalConfig.viscous) {
            exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, to!int(gtl0), to!int(ftl0));
        }

        // Phase 03 LOCAL
        try {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl0, ftl0); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl0, ftl0); }
                }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 03 of implicit update: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 04 MPI
        // We've put this detector step here because it needs the ghost-cell data
        // to be current, as it should be just after a call to apply_convective_bc().
        if ((GlobalConfig.do_shock_detect) &&
            ((!GlobalConfig.frozen_shock_detector) ||
             (GlobalConfig.shock_detector_freeze_step > SimState.step))) {
            detect_shocks(gtl0, ftl0);
        }
        // Phase 05 LOCAL
        try {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                //
                // Make sure that the workspace is of the correct size.
                auto cqi = blk.myConfig.cqi;
                if (!blk.crhs || blk.crhs.nrows != cqi.n || blk.crhs.ncols != (cqi.n+1)) {
                    // Augmented matrix for our linear system.
                    blk.crhs = new Matrix!double(cqi.n, cqi.n+1);
                }
                if (blk.U0save.length != cqi.n) { blk.U0save = new_ConservedQuantities(cqi.n); }
                if (blk.RU0.length != cqi.n) { blk.RU0 = new_ConservedQuantities(cqi.n); }
                if (blk.dRUdU.length != cqi.n) { blk.dRUdU.length = cqi.n; }
                //
                // Now, work through the cells, doing the update.
                // Note that the only difference between the backward-Euler and the implicit-RK1 schemes
                // is the factor of 2 that appears in 2 places.  We call this M.
                double M = (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.implicit_rk1) ? 2.0 : 1.0;
                bool allow_hoi_rhs = allow_high_order_interpolation;
                bool allow_hoi_matrix = allow_high_order_interpolation && GlobalConfig.allow_interpolation_for_sensitivity_matrix;
                bool blklocal_with_local_time_stepping = with_local_time_stepping;
                double blklocal_dt_global = SimState.dt_global;
                double blklocal_t0 = SimState.time;
                foreach (cell; blk.cells) {
                    double dt = (blklocal_with_local_time_stepping) ? cell.dt_local : blklocal_dt_global;
                    auto dUdt0 = cell.dUdt[0];
                    auto U0 = cell.U[0];
                    auto U1 = cell.U[1];
                    //
                    // Set up the linear system by evaluating the sensitivity matrix.
                    // Do this without high-order interpolation, just to be cheap.
                    version(complex_numbers) { U0.clear_imaginary_components(); }
                    blk.U0save.copy_values_from(U0);
                    cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                    dUdt0.clear();
                    blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_matrix, reaction_fraction);
                    bool allFinite = true;
                    foreach (k; 0 .. cqi.n) { if (!isFinite(dUdt0[k].re)) { allFinite = false; } }
                    if (!allFinite) {
                        debug { writeln("Unperturbed state U0=", U0, " dUdt0=", dUdt0); }
                        throw new Error("While evaluating initial R(U), not all dUdt elements are finite.");
                    }
                    blk.RU0.copy_values_from(dUdt0);
                    foreach (j; 0 .. cqi.n) {
                        U0.copy_values_from(blk.U0save);
                        // Perturb one quantity and get the derivative vector d(R(U))/dU[j].
                        version(complex_numbers) {
                            // Use a small enough perturbation such that the error
                            // in first derivative will be much smaller then machine precision.
                            double h = 1.0e-30;
                            number hc = Complex!double(0.0, h);
                            // Perturb one quantity.
                            U0[j] += hc;
                            cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                            // Get derivative vector.
                            dUdt0.clear();
                            blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_matrix, reaction_fraction);
                            foreach (k; 0 .. cqi.n) { if (!isFinite(dUdt0[k].re)) { allFinite = false; } }
                            if (!allFinite) {
                                debug { writeln("Perturbation j=", j, " U0", U0, " dUdt0=", dUdt0); }
                                throw new Error("While evaluating perturbed R(U), Not all dUdt elements are finite.");
                            }
                            foreach (k; 0 .. cqi.n) {
                                blk.dRUdU[k] = dUdt0[k].im / h;
                            }
                        } else {
                            // Scale the perturbation on the magnitude of the conserved quantity.
                            double h = (blk.myConfig.perturbation_for_real_differences)*(fabs(U0[j]) + 1.0);
                            // Perturb one quantity.
                            U0[j] += h;
                            cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                            // Get derivative vector.
                            dUdt0.clear();
                            blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_matrix, reaction_fraction);
                            foreach (k; 0 .. cqi.n) { if (!isFinite(dUdt0[k])) { allFinite = false; } }
                            if (!allFinite) {
                                debug { writeln("Perturbation j=", j, " U0", U0, " dUdt0=", dUdt0); }
                                throw new Error("While evaluating perturbed R(U), Not all dUdt elements are finite.");
                            }
                            foreach (k; 0 .. cqi.n) {
                                blk.dRUdU[k] = (dUdt0[k] - blk.RU0[k])/h;
                            }
                        }
                        // Assemble coefficients of the linear system in the augmented matrix.
                        foreach (k; 0 .. cqi.n) {
                            blk.crhs[k,j] = ((k==j) ? M/dt : 0.0) - blk.dRUdU[k];
                        }
                    }
                    // Evaluate the right-hand side of the linear system equations.
                    U0.copy_values_from(blk.U0save);
                    cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                    dUdt0.clear();
                    blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_rhs, reaction_fraction);
                    foreach (k; 0 .. cqi.n) { blk.crhs[k,cqi.n] = dUdt0[k].re; }
                    // Solve for dU and update U.
                    gaussJordanElimination!double(blk.crhs);
                    foreach (j; 0 .. cqi.n) { U1[j] = U0[j] + M*blk.crhs[j,cqi.n]; }
                    //
                    version(turbulence) {
                        foreach(j; 0 .. cqi.n_turb){
                            U1[cqi.rhoturb+j] = fmax(U1[cqi.rhoturb+j], U0[cqi.mass] * blk.myConfig.turb_model.turb_limits(j));
                        }
                    }
                    // [TODO] PJ 2021-05-15 MHD bits
                    cell.decode_conserved(gtl0, ftl1, blk.omegaz);
                }
                local_invalid_cell_count[i] = blk.count_invalid_cells(gtl0, ftl1);
            } // end foreach blk
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 05 of implicit update: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        //
        flagTooManyBadCells = 0;
        foreach (i, blk; localFluidBlocksBySize) { // serial loop
            if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                flagTooManyBadCells = 1;
                writefln("Following implicit gasdynamic update: %d bad cells in block[%d].",
                         local_invalid_cell_count[i], i);
            }
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (flagTooManyBadCells > 0) {
            throw new FlowSolverException("Too many bad cells following implicit gasdynamic update.");
        }
        //
        if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
            // Phase 06 LOCAL
            // Next do solid domain update IMMEDIATELY after at same flow time leve
            // exchange_ghost_cell_gas_solid_boundary_data();
            exchange_ghost_cell_solid_boundary_data();
            try {
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl0); }
                    }
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl0); }
                    }
                } else {
                    foreach (sblk; localSolidBlocks) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl0); }
                    }
                    foreach (sblk; localSolidBlocks) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl0); }
                    }
                }
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
                    sblk.averageTemperatures();
                    sblk.averageProperties();
                    sblk.clearSources();
                    sblk.computeSpatialDerivatives(ftl0);
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 06 of implicit update: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 07 MPI
            exchange_ghost_cell_solid_boundary_data();
            // Phase 08 LOCAL
            try {
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
                    sblk.averageTGradients();
                    sblk.computeFluxes();
                }
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl0); }
                    }
                } else {
                    foreach (sblk; localSolidBlocks) {
                        if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl0); }
                    }
                }
                // We need to synchronise before updating
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    foreach (scell; sblk.cells) {
                        if (GlobalConfig.udfSolidSourceTerms) {
                            addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
                        }
                        scell.timeDerivatives(ftl0, GlobalConfig.dimensions);
                        scell.stage1Update(SimState.dt_global);
                        scell.ss.e = scell.e[ftl1];
                        sblk.stm.updateTemperature(scell.ss);
                        scell.T = scell.ss.T;
                    } // end foreach scell
                } // end foreach sblk
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 08 of implicit update: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if tight solid domain coupling.
    } while (step_failed && (attempt_number < GlobalConfig.max_attempts_for_step));
    //
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        string msg = format("Implicit update failed after %d attempts; giving up.",
                            GlobalConfig.max_attempts_for_step);
        throw new FlowSolverException(msg);
    }
    //
    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            assert (end_indx == 1, "Unexpected value for end_index.");
            foreach (cell; blk.cells) { swap(cell.U[0], cell.U[end_indx]); }
        }
    } // end foreach blk
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
        foreach (sblk; localSolidBlocks) {
            if (sblk.active) {
                size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
                foreach (scell; sblk.cells) { scell.e[0] = scell.e[end_indx]; }
            }
        } // end foreach sblk
    }
    //
    // Finally, update the globally-known simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_implicit_increment_on_fixed_grid()


void gasdynamic_implicit_increment_with_moving_grid()
{
    shared double t0 = SimState.time;
    assert (!GlobalConfig.with_local_time_stepping, "Cannot use cell-local time step with moving grid.");
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    //
    immutable int ftl0 = 0;
    immutable int ftl1 = 1;
    immutable int gtl0 = 0;
    immutable int gtl1 = 1;
    int step_failed = 0; // Use int because we want to reduce across MPI ranks.
    //
    if (GlobalConfig.udf_source_terms) {
        // Phase 00 LOCAL
        try {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                double blklocal_sim_time = SimState.time;
                foreach (cell; blk.cells) {
                    size_t i_cell = cell.id; size_t j_cell = 0; size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.to_ijk_indices_for_cell(cell.id);
                        i_cell = ijk_indices[0]; j_cell = ijk_indices[1]; k_cell = ijk_indices[2];
                    }
                    getUDFSourceTermsForCell(blk.myL, cell, gtl0, blklocal_sim_time,
                                             blk.myConfig, blk.id, i_cell, j_cell, k_cell);
                }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in Phase 00 of implicit update with miving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            throw new FlowSolverException("Implicit update with moving grid failed when setting source terms. Giving up.");
        }
    }
    // Phase 01 (maybe) MPI
    // Determine grid vertex velocities.
    final switch(GlobalConfig.grid_motion) {
    case GridMotion.none:
        throw new Error("We should not be here without grid motion.");
    case GridMotion.user_defined:
        // Rely on user to set vertex velocities.
        // Note that velocities remain unchanged if the user does nothing.
        assign_vertex_velocities_via_udf(SimState.time, SimState.dt_global);
        break;
    case GridMotion.shock_fitting:
        if (SimState.time > GlobalConfig.shock_fitting_delay) {
            foreach (i, fba; fluidBlockArrays) {
                if (fba.shock_fitting) { compute_vtx_velocities_for_sf(fba); }
            }
        }
        break;
    version(FSI) {
        case GridMotion.FSI:
            foreach (FEMModel; FEMModels) {
                if (SimState.step % FEMModel.myConfig.couplingStep == 0) {
                    FEMModel.compute_vtx_velocities_for_FSI(SimState.dt_global);
                }
            }
            break;
        }
    } // end switch grid_motion
    //
    // Phase 02 LOCAL
    // Compute the consequences of the vertex motion on the grid and conserved quantities.
    try {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (!blk.active) continue;
            auto sblk = cast(SFluidBlock) blk;
            assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
            // Move vertices.
            predict_vertex_positions(sblk, SimState.dt_global, gtl0);
            foreach (vtx; blk.vertices) {
                version(complex_numbers) { vtx.pos[1].clear_imaginary_components(); }
            }
            // Recalculate cell geometry with new vertex positions.
            blk.compute_primary_cell_geometric_data(gtl1);
            blk.compute_least_squares_setup(gtl1);
            // Determine interface velocities using GCL.
            set_gcl_interface_properties(sblk, gtl1, SimState.dt_global);
        }
    } catch (Exception e) {
        debug { writefln("Exception thrown in Phase 02 of implicit update with moving grid: %s", e.msg); }
        step_failed = 1;
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        throw new FlowSolverException("Implicit update failed when computing geometric changes. Giving up.");
    }
    //
    // Begin the gasdynamic update.
    double reaction_fraction = GlobalConfig.reaction_fraction_schedule.interpolate_value(SimState.time);
    int attempt_number = 0;
    int flagTooManyBadCells = 0;
    do {
        ++attempt_number;
        step_failed = 0;
        // Phase 03 LOCAL Preparation for the gas-dynamic flow update.
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.clear_fluxes_of_conserved_quantities();
                foreach (cell; blk.cells) {
                    cell.clear_source_vector();
                    cell.data_is_bad = false;
                }
            }
        }
        // Phase 04 (maybe) MPI
        exchange_ghost_cell_geometry_data();
        exchange_ghost_cell_boundary_data(SimState.time, gtl1, ftl0);
        exchange_ghost_cell_gas_solid_boundary_data();
        if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
            exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl0, ftl0);
        }
        if (GlobalConfig.viscous) {
            exchange_ghost_cell_boundary_viscous_gradient_data(SimState.time, to!int(gtl0), to!int(ftl0));
        }

        // Phase 05 LOCAL
        try {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl0, ftl0); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl0, ftl0); }
                }
            }
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.set_face_flowstates_to_averages_from_cells(); }
            }
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 05 of implicit update with moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        // Phase 06 MPI
        // We've put this detector step here because it needs the ghost-cell data
        // to be current, as it should be just after a call to apply_convective_bc().
        if ((GlobalConfig.do_shock_detect) &&
            ((!GlobalConfig.frozen_shock_detector) ||
             (GlobalConfig.shock_detector_freeze_step > SimState.step))) {
            detect_shocks(gtl0, ftl0);
        }
        // Phase 07 LOCAL
        try {
            foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
                if (!blk.active) continue;
                //
                // Make sure that the workspace is of the correct size.
                auto cqi = blk.myConfig.cqi;
                if (!blk.crhs || blk.crhs.nrows != cqi.n || blk.crhs.ncols != (cqi.n+1)) {
                    // Augmented matrix for our linear system.
                    blk.crhs = new Matrix!double(cqi.n, cqi.n+1);
                }
                if (blk.U0save.length != cqi.n) { blk.U0save = new_ConservedQuantities(cqi.n); }
                if (blk.RU0.length != cqi.n) { blk.RU0 = new_ConservedQuantities(cqi.n); }
                if (blk.dRUdU.length != cqi.n) { blk.dRUdU.length = cqi.n; }
                //
                // Now, work through the cells, doing the update.
                // Note that the only difference between the backward-Euler and the implicit-RK1 schemes
                // is the factor of 2 that appears in 2 places.  We call this M.
                double M = (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.implicit_rk1) ? 2.0 : 1.0;
                bool allow_hoi_rhs = allow_high_order_interpolation;
                bool allow_hoi_matrix = allow_high_order_interpolation && GlobalConfig.allow_interpolation_for_sensitivity_matrix;
                double dt = SimState.dt_global;
                double blklocal_t0 = SimState.time;
                foreach (cell; blk.cells) {
                    auto dUdt0 = cell.dUdt[0];
                    auto U0 = cell.U[0];
                    auto U1 = cell.U[1];
                    //
                    // Set up the linear system by evaluating the sensitivity matrix.
                    // Do this without high-order interpolation, just to be cheap.
                    version(complex_numbers) { U0.clear_imaginary_components(); }
                    blk.U0save.copy_values_from(U0);
                    cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                    dUdt0.clear();
                    blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_matrix, reaction_fraction);
                    bool allFinite = true;
                    foreach (k; 0 .. cqi.n) { if (!isFinite(dUdt0[k].re)) { allFinite = false; } }
                    if (!allFinite) {
                        debug { writeln("Unperturbed state U0=", U0, " dUdt0=", dUdt0); }
                        throw new Error("While evaluating initial R(U), not all dUdt elements are finite.");
                    }
                    blk.RU0.copy_values_from(dUdt0);
                    foreach (j; 0 .. cqi.n) {
                        U0.copy_values_from(blk.U0save);
                        // Perturb one quantity and get the derivative vector d(R(U))/dU[j].
                        version(complex_numbers) {
                            // Use a small enough perturbation such that the error
                            // in first derivative will be much smaller then machine precision.
                            double h = 1.0e-30;
                            number hc = Complex!double(0.0, h);
                            // Perturb one quantity.
                            U0[j] += hc;
                            cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                            // Get derivative vector.
                            dUdt0.clear();
                            blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_matrix, reaction_fraction);
                            foreach (k; 0 .. cqi.n) { if (!isFinite(dUdt0[k].re)) { allFinite = false; } }
                            if (!allFinite) {
                                debug { writeln("Perturbation j=", j, " U0", U0, " dUdt0=", dUdt0); }
                                throw new Error("While evaluating perturbed R(U), Not all dUdt elements are finite.");
                            }
                            foreach (k; 0 .. cqi.n) {
                                blk.dRUdU[k] = dUdt0[k].im / h;
                            }
                        } else {
                            // Scale the perturbation on the magnitude of the conserved quantity.
                            double h = (blk.myConfig.perturbation_for_real_differences)*(fabs(U0[j]) + 1.0);
                            // Perturb one quantity.
                            U0[j] += h;
                            cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                            // Get derivative vector.
                            dUdt0.clear();
                            blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_matrix, reaction_fraction);
                            foreach (k; 0 .. cqi.n) { if (!isFinite(dUdt0[k])) { allFinite = false; } }
                            if (!allFinite) {
                                debug { writeln("Perturbation j=", j, " U0", U0, " dUdt0=", dUdt0); }
                                throw new Error("While evaluating perturbed R(U), Not all dUdt elements are finite.");
                            }
                            foreach (k; 0 .. cqi.n) {
                                blk.dRUdU[k] = (dUdt0[k] - blk.RU0[k])/h;
                            }
                        }
                        // Assemble coefficients of the linear system in the augmented matrix.
                        foreach (k; 0 .. cqi.n) {
                            blk.crhs[k,j] = ((k==j) ? M/dt : 0.0) - blk.dRUdU[k];
                        }
                    }
                    // Evaluate the right-hand side of the linear system equations.
                    U0.copy_values_from(blk.U0save);
                    cell.decode_conserved(gtl0, ftl0, blk.omegaz);
                    dUdt0.clear();
                    blk.evalRU(blklocal_t0, gtl0, ftl0, cell, allow_hoi_rhs, reaction_fraction);
                    foreach (k; 0 .. cqi.n) { blk.crhs[k,cqi.n] = dUdt0[k].re; }
                    // Solve for dU and update U.
                    gaussJordanElimination!double(blk.crhs);
                    // We apply the GCL scaling, also.
                    double vr = cell.volume[gtl0].re / cell.volume[gtl1].re;
                    foreach (k; 0 .. cqi.n) { U1[k] = to!number(vr*(U0[k].re + M*blk.crhs[k,cqi.n])); }
                    //
                    version(turbulence) {
                        foreach(k; 0 .. cqi.n_turb){
                            U1[cqi.rhoturb+k] = fmax(U1[cqi.rhoturb+k], U1[cqi.mass] * blk.myConfig.turb_model.turb_limits(k));
                        }
                    }
                    // [TODO] PJ 2021-05-15 MHD bits
                    cell.decode_conserved(gtl1, ftl1, blk.omegaz);
                }
                local_invalid_cell_count[i] = blk.count_invalid_cells(gtl1, ftl1);
            } // end foreach blk
        } catch (Exception e) {
            debug { writefln("Exception thrown in phase 07 of implicit update with moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // Start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        //
        flagTooManyBadCells = 0;
        foreach (i, blk; localFluidBlocksBySize) { // serial loop
            if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                flagTooManyBadCells = 1;
                writefln("Following implicit gasdynamic update: %d bad cells in block[%d].",
                         local_invalid_cell_count[i], i);
            }
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (flagTooManyBadCells > 0) {
            throw new FlowSolverException("Too many bad cells following implicit gasdynamic update with moving grid.");
        }
        //
        if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
            // Phase 08 LOCAL
            // Next do solid domain update IMMEDIATELY after at same flow time leve
            // exchange_ghost_cell_gas_solid_boundary_data();
            exchange_ghost_cell_solid_boundary_data();
            try {
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl0); }
                    }
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl0); }
                    }
                } else {
                    foreach (sblk; localSolidBlocks) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl0); }
                    }
                    foreach (sblk; localSolidBlocks) {
                        if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryCells(SimState.time, ftl0); }
                    }
                }
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
                    sblk.averageTemperatures();
                    sblk.averageProperties();
                    sblk.clearSources();
                    sblk.computeSpatialDerivatives(ftl0);
                }
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 08 of implicit update with moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
            // Phase 09 MPI
            exchange_ghost_cell_solid_boundary_data();
            // Phase 10 LOCAL
            try {
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
                    sblk.averageTGradients();
                    sblk.computeFluxes();
                }
                if (GlobalConfig.apply_bcs_in_parallel) {
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl0); }
                    }
                } else {
                    foreach (sblk; localSolidBlocks) {
                        if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl0); }
                    }
                }
                // We need to synchronise before updating
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    foreach (scell; sblk.cells) {
                        if (GlobalConfig.udfSolidSourceTerms) {
                            addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time, sblk);
                        }
                        scell.timeDerivatives(ftl0, GlobalConfig.dimensions);
                        scell.stage1Update(SimState.dt_global);
                        scell.ss.e = scell.e[ftl1];
                        sblk.stm.updateTemperature(scell.ss);
                        scell.T = scell.ss.T;
                    } // end foreach scell
                } // end foreach sblk
            } catch (Exception e) {
                debug { writefln("Exception thrown in phase 10 of implicit update with moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if tight solid domain coupling.
    } while (step_failed && (attempt_number < GlobalConfig.max_attempts_for_step));
    //
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        string msg = format("Implicit update with moving grid failed after %d attempts; giving up.",
                            GlobalConfig.max_attempts_for_step);
        throw new FlowSolverException(msg);
    }
    //
    // Get the end conserved data into U[0] for next step and
    // update the latest grid level to the new step grid level 0 and recalculate geometry.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            assert (end_indx == 1, "Unexpected value for end_index.");
            foreach (cell; blk.cells) {
                swap(cell.U[0], cell.U[end_indx]);
                cell.copy_grid_level_to_level(gtl1, gtl0);
            }
            blk.compute_primary_cell_geometric_data(gtl0);
            blk.compute_least_squares_setup(gtl0);
        }
    } // end foreach blk
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
        foreach (sblk; localSolidBlocks) {
            if (sblk.active) {
                size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
                foreach (scell; sblk.cells) { scell.e[0] = scell.e[end_indx]; }
            }
        } // end foreach sblk
    }
    //
    // Finally, update the globally-known simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_implicit_increment_with_moving_grid()


//---------------------------------------------------------------------------

void detect_shocks(int gtl, int ftl)
// calculate if faces/cells are considered to be influenced by a discontinuity
{
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            blk.detect_shock_points();
            blk.shock_faces_to_cells();
        }
    }
    // perform an iterative diffusion process to spread the influence of the shock detector
    if (GlobalConfig.shock_detector_smoothing) {
        foreach(i; 0 .. GlobalConfig.shock_detector_smoothing) {
            exchange_ghost_cell_shock_data(SimState.time, gtl, ftl);
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.diffuse_shock_marker(); }
            }
        }
        // mark faces as shocked based on cell values
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.shock_cells_to_faces(); }
        }
    }
    if (GlobalConfig.strict_shock_detector) {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.enforce_strict_shock_detector(); }
        }
    }
} // end detect_shocks
