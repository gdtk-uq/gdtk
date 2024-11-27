/** simcore.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Peter J. and Rowan G.
 * First code: 2015-02-05
 * History:
 *   2024-02-12 -- moved to lmr5 area
 *                 init simulation routines are time-marching specific, so moved to that module
 *                 integrate-in-time moved to timemarching module
 *                 march-over-blocks removed; will be reinstated in own module.
 */

module simcore;

import core.memory;
import std.math;
import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.typecons;
import std.datetime;
import std.parallelism;
import std.json;
import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import lua_helper;
import fileutil;
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
import lmr.loads;
import conservedquantities;
import special_block_init;
import efield;
version (opencl_gpu_chem) {
    import opencl_gpu_chem;
}
version (cuda_gpu_chem) {
    import cuda_gpu_chem;
}
version(mpi_parallel) {
    import mpi;
}
version(FSI) { import fsi; }

import simcore_gasdynamic_step;
import simcore_solid_step;
import simcore_exchange;
import util.time_utils;
import init : readControl;

// The shared double[] flavour of GlobalConfig.userPad can give trouble,
// so we need a normal array for the MPI task to work with.
double[] userPad_copy;

//----------------------------------------------------------------------------

void check_run_time_configuration(double target_time_as_requested)
{
    // Alter configuration setting if necessary.
    if (GlobalConfig.control_count > 0 && (SimState.step % GlobalConfig.control_count) == 0) {
        readControl(); // Reparse the time-step control parameters occasionally.
        SimState.target_time = (GlobalConfig.block_marching) ? target_time_as_requested : GlobalConfig.max_time;
    }
    if (GlobalConfig.viscous && GlobalConfig.viscous_factor < 1.0 &&
        SimState.time > GlobalConfig.viscous_delay) {
        // We want to increment the viscous_factor that scales the viscous effects terms.
        double viscous_factor = GlobalConfig.viscous_factor;
        viscous_factor += GlobalConfig.viscous_factor_increment;
        viscous_factor = min(viscous_factor, 1.0);
        // Make sure that everyone is up-to-date.
        foreach (myblk; localFluidBlocksBySize) {
            myblk.myConfig.viscous_factor = viscous_factor;
        }
        GlobalConfig.viscous_factor = viscous_factor;
    }
    // We might need to activate or deactivate the IgnitionZones depending on
    // what simulation time we are up to.
    if (SimState.time >= GlobalConfig.ignition_time_start && SimState.time <= GlobalConfig.ignition_time_stop) {
        foreach (blk; localFluidBlocksBySize) { blk.myConfig.ignition_zone_active = true; }
        GlobalConfig.ignition_zone_active = true;
    } else {
        foreach (blk; localFluidBlocksBySize) { blk.myConfig.ignition_zone_active = false; }
        GlobalConfig.ignition_zone_active = false;
    }
} // end check_run_time_configuration()

void synchronize_corner_coords_for_all_blocks()
{
    // Presently, the corner coordinates are only meaningful for structured-grid blocks.
    version(mpi_parallel) {
        // In MPI context, we can only see a subset of the block data.
        foreach (blk; globalBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (!sblk) { continue; }
            if (canFind(GlobalConfig.localFluidBlockIds, sblk.id)) {
                // We can see this inside this block to get valid coordinate values.
                sblk.copy_current_corner_coords();
            } else {
                // Cannot see this block so fill in invalid coordinate values.
                sblk.set_current_corner_coords_to_infinity();
            }
        }
        // Now, propagate the valid coordinates across all tasks.
        foreach (blk; globalBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (sblk) {
                MPI_Allreduce(MPI_IN_PLACE, sblk.corner_coords.ptr, sblk.corner_coords.length,
                              MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
        }
    } else {
        // In shared-memory, we can see all blocks.
        foreach (blk; globalBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (sblk) { sblk.copy_current_corner_coords(); }
        }
    }
} // end synchronize_corner_coords_for_all_blocks()

void call_UDF_at_timestep_start()
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "atTimestepStart");
    if (lua_isnil(L, -1)) {
        // There is no suitable Lua function.
        lua_pop(L, 1); // discard the nil item
    } else {
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
        // Proceed to call the user's function.
        lua_pushnumber(L, SimState.time);
        lua_pushnumber(L, SimState.step);
        lua_pushnumber(L, SimState.dt_global);
        int number_args = 3;
        int number_results = 0;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            string errMsg = "ERROR: while running user-defined function atTimestepStart()\n";
            errMsg ~= to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        lua_getglobal(L, "dt_override");
        if (lua_isnumber(L, -1)) {
            SimState.dt_override = to!double(lua_tonumber(L, -1));
        } else {
            SimState.dt_override = 0.0;
        }
        lua_pop(L, 1); // dispose dt_override item
        //
        if (GlobalConfig.user_pad_length > 0) {
            fill_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_timestep_start()

void broadcast_master_userPad()
{
    version(mpi_parallel) {
        // The userPad data in master MPI task is broadcast to all other MPI tasks.
        int nelem = to!int(GlobalConfig.userPad.length);
        assert(nelem == GlobalConfig.user_pad_length, "Oops, wrong lengths");
        // We allocate the array once.
        if (userPad_copy.length < nelem) { userPad_copy.length = nelem; }
        foreach (i, elem; GlobalConfig.userPad) { userPad_copy[i] = elem; }
        MPI_Bcast(userPad_copy.ptr, nelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (GlobalConfig.mpi_rank_for_local_task > 0) {
            foreach (i, elem; userPad_copy) { GlobalConfig.userPad[i] = elem; }
        }
    }
} // end broadcast_master_userPad()

void copy_userPad_into_block_interpreters()
{
    // Within the one task, broadcast userPad to the Lua interpreters
    // associated with the blocks and boundary-conditions.
    foreach (blk; localFluidBlocks) {
        push_array_to_Lua(blk.myL, GlobalConfig.userPad, "userPad");
        foreach (bc; blk.bc) {
            if (bc.myL) { push_array_to_Lua(bc.myL, GlobalConfig.userPad, "userPad"); }
        }
    }
} // end copy_userPad_into_block_interpreters()

void call_UDF_at_timestep_end()
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "atTimestepEnd");
    if (lua_isnil(L, -1)) {
        // There is no suitable Lua function.
        lua_pop(L, 1); // discard the nil item
    } else {
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
        // Proceed to call the user's function.
        lua_pushnumber(L, SimState.time);
        lua_pushnumber(L, SimState.step);
        lua_pushnumber(L, SimState.dt_global);
        int number_args = 3;
        int number_results = 0;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            string errMsg = "ERROR: while running user-defined function atTimestepEnd()\n";
            errMsg ~= to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        if (GlobalConfig.user_pad_length > 0) {
            fill_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_timestep_end()

void call_UDF_at_write_to_file()
{
    // The user may also have some writing of data to do via their Lua script file.
    // This function is called just after writing the flow solution to file.
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "atWriteToFile");
    if (lua_isnil(L, -1)) {
        // There is no suitable Lua function.
        lua_pop(L, 1); // discard the nil item
    } else {
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
        //
        // Proceed to call the user's function.
        lua_pushnumber(L, SimState.time);
        lua_pushnumber(L, SimState.step);
        lua_pushnumber(L, SimState.dt_global);
        int number_args = 3;
        int number_results = 0;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            string errMsg = "ERROR: while running user-defined function atWriteToFile()\n";
            errMsg ~= to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        if (GlobalConfig.user_pad_length > 0) {
            fill_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_write_to_file()

void update_ch_for_divergence_cleaning()
{
    bool first = true;
    foreach (blk; localFluidBlocksBySize) {
        if (!blk.active) continue;
        if (first) {
            GlobalConfig.c_h = blk.update_c_h(SimState.dt_global);
            first = false;
        } else {
            GlobalConfig.c_h = fmin(blk.update_c_h(SimState.dt_global), GlobalConfig.c_h);
        }
    }
    version(mpi_parallel) {
        double my_c_h = GlobalConfig.c_h;
        MPI_Allreduce(MPI_IN_PLACE, &my_c_h, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        GlobalConfig.c_h = my_c_h;
    }
    // Now that we have a globally-reduced value, propagate that new value
    // into the block-local config structure.
    foreach (blk; localFluidBlocksBySize) {
        if (!blk.active) continue;
        blk.myConfig.c_h = GlobalConfig.c_h;
    }
} // end update_ch_for_divergence_cleaning()

void set_mu_and_k()
{
    version(turbulence){
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.flow_property_spatial_derivatives(0);
                blk.estimate_turbulence_viscosity();
            }
        }
    }
} // end set_mu_and_k()

void chemistry_step(double dt)
{
    version (gpu_chem) {
        if (GlobalConfig.with_local_time_stepping)
            assert(0, "Oops, GPU accelerated chemistry and LTS aren't currently compatible.");
        GlobalConfig.gpuChem.thermochemical_increment(dt);
    } else {
        // without GPU accelerator
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                if (GlobalConfig.with_local_time_stepping) {
                    foreach (cell; blk.cells) { cell.thermochemical_increment(cell.dt_local); }
                } else {
                    foreach (cell; blk.cells) { cell.thermochemical_increment(dt); }
                }
            }
        }
    }
} // end chemistry_half_step()


void compute_Linf_residuals(ref ConservedQuantities Linf_residuals)
{
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_Linf_residuals();
    }
    Linf_residuals.copy_values_from(localFluidBlocks[0].Linf_residuals);
    auto cqi = GlobalConfig.cqi;
    foreach (blk; localFluidBlocksBySize) {
        foreach (i; 0 .. cqi.n) {
            Linf_residuals[i] = fmax(Linf_residuals[i], fabs(blk.Linf_residuals[i]));
        }
    }
} // end compute_Linf_residuals()

void compute_L2_residual(ref number L2_residual)
{
    L2_residual = to!number(0.0);
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_L2_residual();
    }
    foreach (blk; localFluidBlocksBySize) {
        L2_residual += blk.L2_residual;
    }
} // end compute_L2_residuals()

/* This version used in timemarching.d */
void compute_mass_balance(ref number mass_balance)
{
    mass_balance = to!number(0.0);
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_mass_balance();
    }
    foreach (blk; localFluidBlocksBySize) {
        mass_balance += blk.mass_balance;
    }
} // end compute_mass_balance()

/* This version used in newtonkrylovsolver.d */
double compute_mass_balance()
{
    double mass_balance = 0.0;
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_mass_balance();
    }
    foreach (blk; localFluidBlocksBySize) {
        mass_balance += blk.mass_balance;
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(mass_balance), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    return fabs(mass_balance);
} // end compute_mass_balance()

void compute_wall_distances() {
    /*
        Compute the distance from each cell to the nearest viscous wall,
        for any turbulence models that need it.

        @author: Nick Gibbons
    */
    // First count many viscous wall faces are in our local blocks
    int nfaces = 0;
    foreach (blk; localFluidBlocksBySize) {
        foreach(bc; blk.bc) {
            if (bc.is_wall_with_viscous_effects) nfaces+=bc.faces.length;
        }
    }
    //
    // Now pack their centroid positions and normal vectors and lengths into a special buffer
    double[] facepos;
    size_t ii=0;
    facepos.length = nfaces*7;
    int this_rank = GlobalConfig.mpi_rank_for_local_task;
    //
    foreach(blk; localFluidBlocksBySize) {
        foreach(bc; blk.bc) {
            if (bc.is_wall_with_viscous_effects){
                foreach(face; bc.faces){
                    version(complex_numbers){
                    facepos[ii++] = face.pos.x.re;
                    facepos[ii++] = face.pos.y.re;
                    facepos[ii++] = face.pos.z.re;
                    facepos[ii++] = face.n.x.re;
                    facepos[ii++] = face.n.y.re;
                    facepos[ii++] = face.n.z.re;
                    facepos[ii++] = face.length.re;
                    } else {
                    facepos[ii++] = face.pos.x;
                    facepos[ii++] = face.pos.y;
                    facepos[ii++] = face.pos.z;
                    facepos[ii++] = face.n.x;
                    facepos[ii++] = face.n.y;
                    facepos[ii++] = face.n.z;
                    facepos[ii++] = face.length;
                    }
                }
            }
        }
    }
    // Now we need to accumulate faces from across the other MPI tasks
    version(mpi_parallel) {
        int my_rank = GlobalConfig.mpi_rank_for_local_task;
        int mpi_worldsize = GlobalConfig.mpi_size;
        //
        int globalsize = to!int(facepos.length);
        MPI_Allreduce(MPI_IN_PLACE,&globalsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        //
        double[] globalfacepos;
        globalfacepos.length = globalsize;
        //
        int ngathered = 0;
        double[] taskbuffer;
        int ntask;
        foreach (task; 0 .. mpi_worldsize){
            if (my_rank == task) ntask = to!int(facepos.length);
            MPI_Bcast(&ntask, 1, MPI_INT, task, MPI_COMM_WORLD);
            //
            taskbuffer.length = ntask;
            if (my_rank == task) foreach(i; 0 .. ntask) taskbuffer[i] = facepos[i];
            MPI_Bcast(taskbuffer.ptr, ntask, MPI_DOUBLE, task, MPI_COMM_WORLD);
            //
            foreach(i; 0 .. ntask) globalfacepos[ngathered+i] = taskbuffer[i];
            ngathered += ntask;
            taskbuffer.length=0;
        }
        // Now clean up by copying the globaldata back into facepos
        facepos.length = globalfacepos.length;
        foreach(i; 0 .. globalfacepos.length)  facepos[i] = globalfacepos[i];
        nfaces = to!int(facepos.length)/7;
    } // end version(mpi_parallel)
    //
    if (nfaces == 0) {
        throw new Exception("Turbulence model requires wall distance, but no walls found!");
    }
    // These centroids need to be assembled into a kdtree
    size_t totalfaces = nfaces;
    Node[] nodes;
    foreach(i; 0 .. nfaces) {
        Node node = {[facepos[7*i+0], facepos[7*i+1], facepos[7*i+2]]};
        node.cellid = i;
        nodes ~= node;
    }
    auto root = makeTree(nodes);
    //
    // Now loop over the nodes in each of our local blocks and set dwall
    foreach(blk; localFluidBlocksBySize) {
        foreach(cell; blk.cells){
            Node cellnode = {[cell.pos[0].x.re, cell.pos[0].y.re, cell.pos[0].z.re]};
            const(Node)* found = null;
            double bestDist = 0.0;
            size_t nVisited = 0;
            root.fast_nearest(cellnode, 0, found, bestDist, nVisited);
            //double dist = sqrt(bestDist); // This caused problems on distorted grid cells.

            // Unpack face info into a point and a normal, and compute the normal distance to the resulting plane
            double[3] n = [facepos[7*found.cellid+3],
                           facepos[7*found.cellid+4],
                           facepos[7*found.cellid+5]];
            double[3] px = [cell.pos[0].x.re - found.x[0],
                            cell.pos[0].y.re - found.x[1],
                            cell.pos[0].z.re - found.x[2]];
            double dnormal = fabs(px[0]*n[0] + px[1]*n[1] + px[2]*n[2]);
            // We also care about the parallel direction, since we could be far from a wall but still close to its
            // infinite plane. See derivation 28/10/21 by NNG for the new correction factor offwalldistance.
            double celllength = facepos[7*found.cellid+6];
            double dpx_squared = px[0]*px[0]  + px[1]*px[1]  + px[2]*px[2];
            double dparallel = sqrt(dpx_squared - dnormal*dnormal);
            double offwalldistance = fmax(dparallel - celllength/2, 0.0);
            double dwall = sqrt(dnormal*dnormal + offwalldistance*offwalldistance);
            cell.dwall = dwall;
        }
    }
} // end compute_wall_distances()
