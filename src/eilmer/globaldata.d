/**
 * globaldata.d
 *
 * Author: Peter J. and Rowan G.
 * Versions: 2014-07-18 : first cut.
 *           2015-22-04 : added containers for solid blocks
 *           2017-01-17 : Added globalFluidBlocks for the MPI-parallel.
 *           2020-02-11 : Removed globalFluidBlocks & globalSolidBlocks;
 *                        replaced with globalBlocks (Kyle A. Damm).
 */

module globaldata;

import globalconfig;
import fluidblock;
import fluidblockarray;
import ssolidblock;
import block;
import loads;

// For this globally-accessible data we use the __gshared storage class
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

