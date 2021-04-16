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
import nm.complex;
import nm.number;

import geom;
import geom.misc.kdtree;
import gas;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import fluidblock;
import sfluidblock;
import ufluidblock;
import ssolidblock;
import solidprops;
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
    import mpi.util;
}
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
        foreach (i, myblk; parallel(localFluidBlocksBySize,1)) {
            // Note 'i' is not necessarily the block id but
            // that is not important here, just need a unique spot to poke into local_dt_allow.
            if (myblk.active) {
                double[3] results = myblk.determine_time_step_size(SimState.dt_global, (SimState.step > 0));
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
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);

    // compute number of super steps
    bool euler_step = false;
    double dt_global;
    double alpha;
    double s_RKL;
    int S;
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
    SimState.s_RKL = S;
    // --------------------------------------------------
    // j = 1
    // --------------------------------------------------

    // First-stage of gas-dynamic update.
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level remains at zero for the non-moving grid

    // Preparation for the inviscid gas-dynamic flow update.
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
    if (GlobalConfig.apply_bcs_in_parallel) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	}
    } else {
	foreach (blk; localFluidBlocksBySize) {
	    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	}
    }
    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
        exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
        exchange_ghost_cell_gas_solid_boundary_data();
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
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
	int local_ftl = ftl;
	int local_gtl = gtl;
	bool local_with_local_time_stepping = with_local_time_stepping;
	double local_dt_global = dt_global;
	double local_sim_time = SimState.time;
	foreach (cell; blk.cells) {
	    cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
	    if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
		cell.add_viscous_source_vector();
	    }
	    if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
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
		addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
					local_sim_time, blk.myConfig,
					blk.id, i_cell, j_cell, k_cell);
	    }
	    cell.time_derivatives(local_gtl, local_ftl);
	}
	if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
	bool force_euler = false;
	foreach (cell; blk.cells) {
	    if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 || euler_step) {
                cell.rkl1_stage_update_for_flow_on_fixed_grid1(local_dt_global, 1, SimState.s_RKL, false); // RKL1 (j=1)
            } else {
                cell.rkl2_stage_update_for_flow_on_fixed_grid1(local_dt_global, 1, SimState.s_RKL, false); // RKL2 (j=1)
            }
            cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
	} // end foreach cell
	local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
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
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
	// Next do solid domain update IMMEDIATELY after at same flow time level
	foreach (sblk; parallel(localSolidBlocks, 1)) {
	    if (!sblk.active) continue;
	    sblk.averageTemperatures();
	    sblk.clearSources();
	    sblk.computeSpatialDerivatives(ftl);
        }
        exchange_ghost_cell_solid_boundary_data();
        exchange_ghost_cell_gas_solid_boundary_data();
        foreach (sblk; parallel(localSolidBlocks, 1)) {
            if (!sblk.active) continue;
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
	    foreach (scell; sblk.activeCells) {
		if (GlobalConfig.udfSolidSourceTerms) {
		    addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
		}
		scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 || euler_step) {
                    scell.stage1RKL1Update(dt_global, 1, SimState.s_RKL); // RKL1 (j=1)
                } else {
                    scell.stage1RKL2Update(dt_global, 1, SimState.s_RKL); // RKL1 (j=1)
                }
                scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
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
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	    }
	}
        if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
            exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
            exchange_ghost_cell_gas_solid_boundary_data();
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
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
	    int local_ftl = ftl;
	    int local_gtl = gtl;
	    bool local_with_local_time_stepping = with_local_time_stepping;
	    double local_dt_global = dt_global;
	    double local_sim_time = SimState.time;
	    foreach (cell; blk.cells) {
		cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
		if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
		    cell.add_viscous_source_vector();
		}
		if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
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
		    addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
					local_sim_time, blk.myConfig,
					    blk.id, i_cell, j_cell, k_cell);
		}
		cell.time_derivatives(local_gtl, local_ftl);
	    }
	    if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
	    bool force_euler = false;
	    foreach (cell; blk.cells) {
                if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1) {
                    cell.rkl1_stage_update_for_flow_on_fixed_grid2(local_dt_global, j, SimState.s_RKL, false); // RKL1
                } else {
                    cell.rkl2_stage_update_for_flow_on_fixed_grid2(local_dt_global, j, SimState.s_RKL, false);
                }
		cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
	    } // end foreach cell
	    local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
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
	if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
	    // Next do solid domain update IMMEDIATELY after at same flow time level
	    foreach (sblk; parallel(localSolidBlocks, 1)) {
		if (!sblk.active) continue;
                sblk.averageTemperatures();
		sblk.clearSources();
		sblk.computeSpatialDerivatives(ftl);
            }
            exchange_ghost_cell_solid_boundary_data();
            foreach (sblk; parallel(localSolidBlocks, 1)) {
                if (!sblk.active) continue;
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
		foreach (scell; sblk.activeCells) {
		    if (GlobalConfig.udfSolidSourceTerms) {
			addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
		    }
		    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                    if (GlobalConfig.gasdynamic_update_scheme == GasdynamicUpdate.rkl1 || euler_step) {
                        scell.stage2RKL1Update(dt_global, j, SimState.s_RKL); // RKL1 (j=1)
                    } else {
                        scell.stage2RKL2Update(dt_global, j, SimState.s_RKL); // RKL2 (j=1)
                    }
		    scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
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
        if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
            foreach (sblk; parallel(localSolidBlocks,1)) {
                if (sblk.active) {
                    foreach (scell; sblk.activeCells) {
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

    int idx = 2;
    if (euler_step) { idx = 1; }

    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = 1; // time-level holds current solution
            foreach (cell; blk.cells) { cell.U[0].copy_values_from(cell.U[idx]); } //swap(cell.U[0], cell.U[end_indx]); }
        }
    } // end foreach blk
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
        foreach (sblk; localSolidBlocks) {
            if (sblk.active) {
                //size_t end_indx = 1; // time-level holds current solution
                foreach (scell; sblk.activeCells) { scell.e[0] = scell.e[idx]; }
            }
        } // end foreach sblk
    }
    // Finally, update the globally known simulation time for the whole step.
    SimState.time = t0 + dt_global;
} // end sts_gasdynamic_explicit_increment_with_fixed_grid()


void gasdynamic_explicit_increment_with_fixed_grid()
{
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    // Set the time-step coefficients for the stages of the update scheme.
    shared double c2 = 1.0; // default for predictor-corrector update
    shared double c3 = 1.0; // default for predictor-corrector update
    final switch ( GlobalConfig.gasdynamic_update_scheme ) {
    case GasdynamicUpdate.euler:
    case GasdynamicUpdate.pc: c2 = 1.0; c3 = 1.0; break;
    case GasdynamicUpdate.midpoint: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.classic_rk3: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.tvd_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.denman_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.rkl1:
    case GasdynamicUpdate.rkl2: assert(false, "invalid option");
    case GasdynamicUpdate.moving_grid_1_stage:
    case GasdynamicUpdate.moving_grid_2_stage: assert(false, "invalid option");
    }
    int attempt_number = 0;
    int step_failed = 0; // Use int because we want to reduce across MPI ranks.
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
        shared int ftl; // time-level within the overall convective-update
        shared int gtl; // grid time-level remains at zero for the non-moving grid
        int flagTooManyBadCells;
        try {
            // First-stage of gas-dynamic update.
            ftl = 0; gtl = 0;
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
            // We've put this detector step here because it needs the ghost-cell data
            // to be current, as it should be just after a call to apply_convective_bc().
            if ((GlobalConfig.do_shock_detect) &&
                ((!GlobalConfig.frozen_shock_detector) || (GlobalConfig.shock_detector_freeze_step > SimState.step))) {
                detect_shocks(gtl, ftl);
            }
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
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            }
            if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
                exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
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
                int local_ftl = ftl;
                int local_gtl = gtl;
                bool local_with_local_time_stepping = with_local_time_stepping;
                double local_dt_global = SimState.dt_global;
                double local_sim_time = SimState.time;
                foreach (cell; blk.cells) {
                    cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                    if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                        cell.add_viscous_source_vector();
                    }
                    if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
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
                        addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                                local_sim_time, blk.myConfig,
                                                blk.id, i_cell, j_cell, k_cell);
                    }
                    cell.time_derivatives(local_gtl, local_ftl);
                }
                if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
                bool force_euler = false;
                foreach (cell; blk.cells) {
                    cell.stage_1_update_for_flow_on_fixed_grid(local_dt_global, force_euler,
                                                               local_with_local_time_stepping);
                    cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
                } // end foreach cell
                local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
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
                throw new FlowSolverException("Too many bad cells following first-stage gasdynamic update.");
            }
            //
            if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
                // Next do solid domain update IMMEDIATELY after at same flow time leve
                //exchange_ghost_cell_gas_solid_boundary_data();
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
                    sblk.clearSources();
                    sblk.computeSpatialDerivatives(ftl);
                }
                exchange_ghost_cell_solid_boundary_data();
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
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
                    foreach (scell; sblk.activeCells) {
                        if (GlobalConfig.udfSolidSourceTerms) {
                            addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                        }
                        scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                        scell.stage1Update(SimState.dt_global);
                        scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                    } // end foreach scell
                } // end foreach sblk
            } // end if tight solid domain coupling.
        } catch (Exception e) {
            debug { writefln("Exception thrown in first-stage of explicit update: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // First-stage of update has failed for some reason,
            // so start the step over again with a reduced time step.
            SimState.dt_global = SimState.dt_global * 0.2;
            continue;
        }
        //
        if (number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) >= 2) {
            // Preparation for second-stage of gas-dynamic update.
            SimState.time = t0 + c2 * SimState.dt_global;
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.clear_fluxes_of_conserved_quantities();
                    foreach (cell; blk.cells) cell.clear_source_vector();
                }
            }
            try {
                // Second stage of gas-dynamic update.
                ftl = 1;
                // We are relying on exchanging boundary data as a pre-reconstruction activity.
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
                    if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
                }

                // for unstructured blocks we need to transfer the convective gradients before the flux calc
                if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                    exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
                }

                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
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
                if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
                    // we exchange boundary data at this point to ensure the
                    // ghost cells along block-block boundaries have the most
                    // recent mu_t and k_t values.
                    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
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
                    int local_ftl = ftl;
                    int local_gtl = gtl;
                    bool local_with_local_time_stepping = with_local_time_stepping;
                    double local_dt_global = SimState.dt_global;
                    double local_sim_time = SimState.time;
                    foreach (cell; blk.cells) {
                        cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                        if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                            cell.add_viscous_source_vector();
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
                            addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                                    local_sim_time, blk.myConfig,
                                                    blk.id, i_cell, j_cell, k_cell);
                        }
                        cell.time_derivatives(local_gtl, local_ftl);
                    }
                    if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
                    foreach (cell; blk.cells) {
                        cell.stage_2_update_for_flow_on_fixed_grid(local_dt_global, local_with_local_time_stepping);
                        cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
                    } // end foreach cell
                    local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
                } // end foreach blk
                //
                flagTooManyBadCells = 0;
                foreach (i, blk; localFluidBlocksBySize) { // serial loop
                    if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                        flagTooManyBadCells = 1;
                        writefln("Following second-stage gasdynamic update: %d bad cells in block[%d].",
                                 local_invalid_cell_count[i], i);
                    }
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (flagTooManyBadCells > 0) {
                    throw new FlowSolverException("Too many bad cells following second-stage gasdynamic update.");
                }
                //
                if ( GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ) {
                    // Do solid domain update IMMEDIATELY after at same flow time level
                    //exchange_ghost_cell_gas_solid_boundary_data();
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
                        sblk.clearSources();
                        sblk.computeSpatialDerivatives(ftl);
                    }
                    exchange_ghost_cell_solid_boundary_data();
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (!sblk.active) continue;
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
                        foreach (scell; sblk.activeCells) {
                            if (GlobalConfig.udfSolidSourceTerms) {
                                addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                            }
                            scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                            scell.stage2Update(SimState.dt_global);
                            scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                        } // end foreach cell
                    } // end foreach blk
                } // end if tight solid domain coupling.
            } catch (Exception e) {
                debug { writeln("Caught exception at end of second-stage of explicit update: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Second-stage of update has failed for some reason,
                // so start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if number_of_stages_for_update_scheme >= 2
        //
        if (number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) >= 3) {
            // Preparation for third stage of gasdynamic update.
            SimState.time = t0 + c3 * SimState.dt_global;
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.clear_fluxes_of_conserved_quantities();
                    foreach (cell; blk.cells) cell.clear_source_vector();
                }
            }
            try {
                // Third stage of gas-dynamic update.
                ftl = 2;
                // We are relying on exchanging boundary data as a pre-reconstruction activity.
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
                    if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
                }

                // for unstructured blocks we need to transfer the convective gradients before the flux calc
                if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                    exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
                }

                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
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
                if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
                    // we exchange boundary data at this point to ensure the
                    // ghost cells along block-block boundaries have the most
                    // recent mu_t and k_t values.
                    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
                    exchange_ghost_cell_gas_solid_boundary_data();
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
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
                    int local_ftl = ftl;
                    int local_gtl = gtl;
                    bool local_with_local_time_stepping = with_local_time_stepping;
                    double local_dt_global = SimState.dt_global;
                    double local_sim_time = SimState.time;
                    foreach (cell; blk.cells) {
                        cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                        if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                            cell.add_viscous_source_vector();
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
                            addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                                    local_sim_time, blk.myConfig,
                                                    blk.id, i_cell, j_cell, k_cell);
                        }
                        cell.time_derivatives(local_gtl, local_ftl);
                    }
                    if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
                    foreach (cell; blk.cells) {
                        cell.stage_3_update_for_flow_on_fixed_grid(local_dt_global, local_with_local_time_stepping);
                        cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
                    } // end foreach cell
                    local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
                } // end foreach blk
                //
                flagTooManyBadCells = 0;
                foreach (i, blk; localFluidBlocksBySize) { // serial loop
                    if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                        flagTooManyBadCells = 1;
                        writefln("Following third-stage gasdynamic update: %d bad cells in block[%d].",
                                 local_invalid_cell_count[i], i);
                    }
                }
                version(mpi_parallel) {
                    MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                }
                if (flagTooManyBadCells > 0) {
                    throw new FlowSolverException("Too many bad cells following third-stage gasdynamic update.");
                }
                //
                if ( GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ) {
                    // Do solid domain update IMMEDIATELY after at same flow time level
                    exchange_ghost_cell_gas_solid_boundary_data();
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
                        sblk.clearSources();
                        sblk.computeSpatialDerivatives(ftl);
                    }
                    exchange_ghost_cell_solid_boundary_data();
                    foreach (sblk; parallel(localSolidBlocks, 1)) {
                        if (!sblk.active) continue;
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
                        foreach (scell; sblk.activeCells) {
                            if (GlobalConfig.udfSolidSourceTerms) {
                                addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                            }
                            scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                            scell.stage3Update(SimState.dt_global);
                            scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                        } // end foreach cell
                    } // end foreach blk
                } // end if tight solid domain coupling.
            } catch (Exception e) {
                debug { writeln("Caught exception at end of third-stage of explicit update: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Third-stage of update has failed for some reason,
                // so start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if number_of_stages_for_update_scheme >= 3
    } while (step_failed && (attempt_number < 3));
    //
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        throw new FlowSolverException("Explicit update failed after 3 attempts; giving up.");
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
                foreach (scell; sblk.activeCells) { scell.e[0] = scell.e[end_indx]; }
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
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    // Set the time-step coefficients for the stages of the update scheme.
    shared double c2 = 1.0; // same for 1-stage or 2-stage update
    shared double c3 = 1.0; // ditto

    int flagTooManyBadCells;
    shared int ftl; // time-level within the overall convective-update
    shared int gtl; // grid time-level

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
    } // end switch grid_motion

    int attempt_number = 0;
    int step_failed = 0; // Use int because we want to reduce across MPI ranks.
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
        try {
            // First-stage of gas-dynamic update.
            ftl = 0; gtl = 0;
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
            gtl = 1; // update gtl now that grid has moved
            exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
                }
            }
            // And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
            foreach (sblk; localSolidBlocks) {
                if (sblk.active) { sblk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, ftl); }
            }
            foreach (sblk; localSolidBlocks) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
            // We've put this detector step here because it needs the ghost-cell data
            // to be current, as it should be just after a call to apply_convective_bc().
            if (GlobalConfig.do_shock_detect) {
                detect_shocks(gtl, ftl);
            }
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
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
                }
            }
            if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
                // we exchange boundary data at this point to ensure the
                // ghost cells along block-block boundaries have the most
                // recent mu_t and k_t values.
                exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
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
                int local_ftl = ftl;
                int local_gtl = gtl;
                bool local_with_local_time_stepping = with_local_time_stepping;
                double local_dt_global = SimState.dt_global;
                double local_sim_time = SimState.time;
                foreach (cell; blk.cells) {
                    cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                    if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                        cell.add_viscous_source_vector();
                    }
                    if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
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
                        addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                                local_sim_time, blk.myConfig,
                                                blk.id, i_cell, j_cell, k_cell);
                    }
                    cell.time_derivatives(local_gtl, local_ftl);
                    bool force_euler = false;
                    cell.stage_1_update_for_flow_on_moving_grid(local_dt_global, local_with_local_time_stepping);
                    cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
                } // end foreach cell
                local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
            } // end foreach blk
            //
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
            foreach (sblk; localSolidBlocks) {
                if (!sblk.active) continue;
                sblk.averageTemperatures();
                sblk.clearSources();
                sblk.computeSpatialDerivatives(ftl);
                sblk.applyPostFluxAction(SimState.time, ftl);
            }
            exchange_ghost_cell_solid_boundary_data();
            foreach (sblk; parallel(localSolidBlocks, 1)) {
                if (!sblk.active) continue;
                sblk.computeFluxes();
                sblk.applyPostFluxAction(SimState.time, ftl);
                foreach (scell; sblk.activeCells) {
                    if (GlobalConfig.udfSolidSourceTerms) {
                        addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                    }
                    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                    scell.stage1Update(SimState.dt_global);
                    scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                } // end foreach scell
            } // end foreach sblk
        } catch (Exception e) {
            debug { writefln("Exception thrown in first-stage of explicit update with moving grid: %s", e.msg); }
            step_failed = 1;
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (step_failed) {
            // First-stage of update has failed for some reason,
            // so start the step over again with a reduced time step.
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
                ftl = 1;
                gtl = 2;
                // We are relying on exchanging boundary data as a pre-reconstruction activity.
                exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
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
                    if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, 0); }
                    // FIX-ME PJ 2018-07-25 Should this be gtl rather than 0?
                }

                // for unstructured blocks we need to transfer the convective gradients before the flux calc
                if (allow_high_order_interpolation && (GlobalConfig.interpolation_order > 1)) {
                    exchange_ghost_cell_boundary_convective_gradient_data(SimState.time, gtl, ftl);
                }

                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, 0); }
                    // FIX-ME PJ 2018-07-25 Should this be gtl rather than 0?
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
                if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
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
                    // we exchange boundary data at this point to ensure the
                    // ghost cells along block-block boundaries have the most
                    // recent mu_t and k_t values.
                    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
                    foreach (blk; parallel(localFluidBlocksBySize,1)) {
                        if (blk.active) {
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
                    int local_ftl = ftl;
                    int local_gtl = gtl;
                    bool local_with_local_time_stepping = with_local_time_stepping;
                    double local_dt_global = SimState.dt_global;
                    double local_sim_time = SimState.time;
                    foreach (cell; blk.cells) {
                        cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                        if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                            cell.add_viscous_source_vector();
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
                            addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                                    local_sim_time, blk.myConfig,
                                                    blk.id, i_cell, j_cell, k_cell);
                        }
                        cell.time_derivatives(local_gtl, local_ftl);
                        cell.stage_2_update_for_flow_on_moving_grid(local_dt_global, local_with_local_time_stepping);
                        cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
                    } // end foreach cell
                    local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
                } // end foreach blk
                //
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
                foreach (sblk; localSolidBlocks) {
                    if (!sblk.active) continue;
                    sblk.averageTemperatures();
                    sblk.clearSources();
                    sblk.computeSpatialDerivatives(ftl);
                    sblk.applyPostFluxAction(SimState.time, ftl);
                }
                exchange_ghost_cell_solid_boundary_data();
                foreach (sblk; parallel(localSolidBlocks, 1)) {
                    if (!sblk.active) continue;
                    sblk.computeFluxes();
                    sblk.applyPostFluxAction(SimState.time, ftl);
                    foreach (scell; sblk.activeCells) {
                        if (GlobalConfig.udfSolidSourceTerms) {
                            addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                        }
                        scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                        scell.stage2Update(SimState.dt_global);
                        scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                    } // end foreach cell
                } // end foreach blk
            } catch (Exception e) {
                debug { writefln("Exception thrown in second-stage of explicit update with moving grid: %s", e.msg); }
                step_failed = 1;
            }
            version(mpi_parallel) {
                MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (step_failed) {
                // Second-stage of update has failed for some reason,
                // so start the step over again with a reduced time step.
                SimState.dt_global = SimState.dt_global * 0.2;
                continue;
            }
        } // end if number_of_stages_for_update_scheme >= 2
    } while (step_failed && (attempt_number < 3));
    //
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &step_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (step_failed) {
        throw new FlowSolverException("Explicit update with moving grid failed after 3 attempts; giving up.");
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
            foreach (scell; sblk.activeCells) { scell.e[0] = scell.e[end_indx]; }
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
