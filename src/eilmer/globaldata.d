/**
 * globaldata.d
 *
 * Author: Peter J. and Rowan G.
 * Versions: 2014-07-18 : first cut.
 *           2015-22-04 : added containers for solid blocks
 *           2017-01-17 : Added globalFluidBlocks for the MPI-parallel.
 */

module globaldata;

import globalconfig;
import fluidblock;
import ssolidblock;

// Global collections of blocks for the simulation, as a whole.
// For the fluid blocks, not all may be present in the local MPI task (Linux process).
static FluidBlock[] globalFluidBlocks;
static SSolidBlock[] solidBlocks;

// Collections of blocks that we can iterate over in parallel.
// The current (shared-memory) parallel code is based on having one FluidBlock object
// or SSolidBlock object per thread.
// In this context, "local" is within the local MPI task (or Linux process),
// which may have several threads running within it.
static FluidBlock[] localFluidBlocks;
static FluidBlock[] localFluidBlocksBySize; // sorted largest to smallest  

// We also need to have a dedicated set of configuration parameters for each thread
// so that there is no need to have memory barriers guarding their access.
// There will be one of these LocalConfig objects per block.
static LocalConfig[] dedicatedConfig;
static LocalConfig[] dedicatedSolidConfig;
