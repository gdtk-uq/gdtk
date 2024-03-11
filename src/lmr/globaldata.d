/**
 * globaldata.d
 *
 * Author: Peter J. and Rowan G.
 * Versions: 2014-07-18 : first cut.
 *           2015-22-04 : added containers for solid blocks
 *           2017-01-17 : Added globalFluidBlocks for the MPI-parallel.
 *           2020-02-11 : Removed globalFluidBlocks & globalSolidBlocks;
 *                        replaced with globalBlocks (Kyle A. Damm).
 *           2021-03-30 : Added globalFluidBlockIO
 *           2024-02-12 : Moved into lmr5 area, removed globalFluidBlockIO
 */

module globaldata;

import std.datetime;

import globalconfig;
import fluidblock;
import fluidblockarray;
import ssolidblock;
import block;
import efield;
import blockio : BlockIO;
import lmr.loads : RunTimeLoads;
version(FSI) { import fsi; }

// I/O objects (thread-local) for access from any module
BlockIO fluidBlkIO;
BlockIO solidBlkIO;

// State data for simulation.
// Needs to be seen by all of the coordination functions in simcore.d.
final class SimState {
    shared static double time;  // present simulation time, tracked by code
    shared static int step;
    shared static double dt_global;     // simulation time step determined by code
    shared static double dt_allow;      // allowable global time step determined by code
    shared static bool is_restart; // Flag indicates if this run is a restart (true) or fresh start (false).
    // for STS
    shared static double dt_global_parab;
    shared static double dt_allow_parab;
    shared static int s_RKL;
    //
    shared static double cfl_max;      // current max cfl determined by code
    shared static double dt_override = 0.0;  // A positive value will override a larger computed time step.
    shared static double target_time;  // simulate_in_time will work toward this value

    // We want to write sets of output files periodically.
    // The following periods set the cadence for output.
    shared static double t_plot;
    shared static double t_history;
    shared static double t_loads;
    // Once we write some data to files, we don't want to write another set of files
    // until we have done some more stepping.  The following flags help us remember
    // the state of the solution output.
    shared static bool output_just_written = true;
    shared static bool history_just_written = true;
    shared static bool loads_just_written = true;
    // We connect the sets of files to the simulation time at which they were written
    // with an index that gets incremented each time we write a set of files.
    shared static int current_tindx;
    shared static int current_loads_tindx;
    // Keep track of the snapshot output
    shared static int nWrittenSnapshots;
    // For working out how long the simulation has been running.
    static SysTime wall_clock_start;
    static int maxWallClockSeconds;
} // end class SimState


// For the following globally-accessible data we use the __gshared storage class
// to put the variables in the classic global memory space.
// All of the threads may access the data but we take responsibility for
// not allowing the threads to step all over each other.

// Global collections of blocks for the simulation, as a whole.
// For the fluid blocks, not all may be present in the local MPI task (Linux process).
__gshared static Block[] globalBlocks;

// Collections of blocks that we can iterate over in parallel.
// The current (shared-memory) parallel code is based on having one FluidBlock object
// or SSolidBlock object per thread.
// In this context, "local" is within the local MPI task (or Linux process),
// which may have several threads running within it.
__gshared static FluidBlock[] localFluidBlocks;
__gshared static FluidBlock[] localFluidBlocksBySize; // sorted largest to smallest
__gshared static SSolidBlock[] localSolidBlocks;
__gshared static SSolidBlock[] localSolidBlocksBySize; // sorted largest to smallest

// We also need to have a dedicated set of configuration parameters for each thread
// so that there is no need to have memory barriers guarding their access.
// There will be one of these LocalConfig objects per block.
__gshared static LocalConfig[] dedicatedConfig;

// We store the computed run time loads globally
// so that we may make these available in the
// Lua environments.
__gshared static RunTimeLoads[] runTimeLoads;
__gshared static size_t[string] runTimeLoadsByName;

// Shock fitting is coordinated across arrays of FluidBlock objects.
// Here is the storage for that coordination data.
__gshared static FBArray[] fluidBlockArrays;

// An object for solving the electric field across the entire simulation
// It may cooperate with other ElectricFields in other processes
__gshared static ElectricField eField;

// Collection of FEM models for FSI
version(FSI) { __gshared static FEMModel[] FEMModels; }
