/**
 * globaldata.d
 *
 * Author: Peter J. and Rowan G.
 * Versions: 2014-07-18 : first cut.
 *           2015-22-04 : added containers for solid blocks
 */

module globaldata;

import globalconfig;
import block;
import ssolidblock;

// Collections of blocks that we can iterate over in parallel.
// The current (shared-memory) parallel code is based on having one SBlock object
// or SSolidBlock object per thread.
static Block[] localFluidBlocks;
static Block[] localFluidBlocksBySize; // sorted largest to smallest  
static SSolidBlock[] solidBlocks;

// We also need to have a dedicated set of configuration parameters for each thread
// so that there is no need to have memory barriers guarding their access.
static LocalConfig[] dedicatedConfig;
static LocalConfig[] dedicatedSolidConfig;
