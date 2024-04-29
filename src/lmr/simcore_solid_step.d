// simcore_solid_step.d
// 2021-04-15: extracted from simcore.d (PJ)

module simcore_solid_step;

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

import geom;
import geom.misc.kdtree;
import gas;
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
import simcore_exchange;


double determine_solid_time_step_size()
{
    // Set the size of the time step to be the minimum allowed for any active block.
    double solid_dt_allow;
    double local_solid_dt_allow;
    double cfl_value = GlobalConfig.solid_domain_cfl;

    // First, check what each block thinks should be the allowable step size.
    foreach (mysblk; parallel(localSolidBlocks,1)) {
        if (mysblk.active) { local_solid_dt_allow = mysblk.determine_time_step_size(cfl_value); }
    }
    // Second, reduce this estimate across all local blocks.
    solid_dt_allow = double.max; // to be sure it is replaced.
    foreach (i, myblk; localFluidBlocks) { // serial loop
        if (myblk.active) { solid_dt_allow = min(solid_dt_allow, local_solid_dt_allow); }
    }
    version(mpi_parallel) {
        double my_solid_dt_allow = solid_dt_allow;
        MPI_Allreduce(MPI_IN_PLACE, &my_solid_dt_allow, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        solid_dt_allow = my_solid_dt_allow;
    }
    return solid_dt_allow;
} // end determine_solid_time_step_size()


void solid_step(double dt_solid)
{
    //
    // we have currently removed the implicit solid update.
    // Call Nigel's update function here.
    // solid_domains_backward_euler_update(SimState.time, SimState.dt_global);
    //
    // perform an Euler update for the solid domain
    int ftl = 0;
    // Next do solid domain update IMMEDIATELY after at same flow update
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
    //
    foreach (sblk; parallel(localSolidBlocks, 1)) {
        if (!sblk.active) continue;
        sblk.averageTemperatures();
        sblk.clearSources();
        sblk.computeSpatialDerivatives(ftl);
    }
    exchange_ghost_cell_solid_boundary_data();
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
            scell.eulerUpdate(dt_solid);
            scell.ss.e = scell.e[ftl+1];
            sblk.stm.updateTemperature(scell.ss);
            scell.T = scell.ss.T;
        } // end foreach scell
    } // end foreach sblk
    foreach (sblk; localSolidBlocks) {
        if (sblk.active) {
            //size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            size_t end_indx = 1;
            foreach (scell; sblk.cells) { scell.e[0] = scell.e[end_indx]; }
        }
    } // end foreach sblk
} // end solid_step()
