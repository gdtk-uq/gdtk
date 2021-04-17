-- prep-flow.lua
-- A place to put helper functions and classes for flow preparation activities.
-- This script is read by the Eilmer4 program at config and flow preparation time,
-- before reading and processing the user's input script.
--
-- Authors: PJ, RJG, Kyle D., Nick G. and Daryl B.
--
print("Loading prep-flow.lua...")

require 'lua_helper'
deepclone = lua_helper.deepclone
checkAllowedNames = lua_helper.checkAllowedNames

require 'blk_conn'
-- Let's pull the symbols out of the blk_conn module
-- and make them global in this namespace
for k,v in pairs(blk_conn) do
   _G[k] = v
end

require 'bc'
for k,v in pairs(bc) do
   _G[k] = v
end

require 'flowstate'
FlowState = flowstate.FlowState
makeFlowStateFn = flowstate.makeFlowStateFn

require 'fluidblock'
FluidBlock = fluidblock.FluidBlock
SBlock2UBlock = fluidblock.SBlock2UBlock
connectBlocks = fluidblock.connectBlocks
identifyBlockConnections = fluidblock.identifyBlockConnections
FBArray = fluidblock.FBArray
FluidBlockArray = fluidblock.FluidBlockArray

require 'solidblock'
SolidBlock = solidblock.SolidBlock
SolidBlockArray = solidblock.SolidBlockArray

require 'mpi'
mpiDistributeBlocks = mpi.mpiDistributeBlocks
mpiDistributeFBArray = mpi.mpiDistributeFBArray

require 'history'
setHistoryPoint = history.setHistoryPoint
setSolidHistoryPoint = history.setSolidHistoryPoint

require 'zones'
ReactionZone = zones.ReactionZone
IgnitionZone = zones.IgnitionZone
TurbulentZone = zones.TurbulentZone
SuppressReconstructionZone = zones.SuppressReconstructionZone

require 'output'
write_control_file = output.write_control_file
write_config_file = output.write_config_file
write_times_file = output.write_times_file
write_block_list_file = output.write_block_list_file
write_mpimap_file = output.write_mpimap_file
write_fluidBlockArrays_file = output.write_fluidBlockArrays_file
write_shock_fitting_helper_files = output.write_shock_fitting_helper_files

require 'sssoptions'
SteadyStateSolver = sssoptions.SteadyStateSolver
ShapeSensitivityCalculator = sssoptions.ShapeSensitivityCalculator
SolidDomainLooseUpdate = sssoptions.SolidDomainLooseUpdate

require 'prep_check'
initTurbulence = prep_check.initTurbulence
checkCellVolumes = prep_check.checkCellVolumes
perform_spatial_gradient_consistency_check = prep_check.perform_spatial_gradient_consistency_check
warn_if_blocks_not_connected = prep_check.warn_if_blocks_not_connected
check_DFT_settings = prep_check.check_DFT_settings

-- ---------------------------------------------------------------------------------------

-- Storage for later definitions of FluidBlock objects.
-- Note that the index for this array starts at 1, in the Lua way.
-- The block ids start at 0 to be like the indexing inside the D code.
-- Yes, this is confusing...
fluidBlocks = {}
-- Storage for later definitions of FluidBlockArray objects.
fluidBlockArrays = {}
-- We may want to look up the blocks via labels rather than numerical id
-- in user-defined procedures.
-- The following dictionaries store the connections.
fluidBlocksDict = {}
fluidBlockArraysDict = {}

-- The user may assign the MPI task id for eack block manually
-- but, if they don't, a default distribution will be made.
mpiTasks = nil

-- Storgage for later definitions of SolidBlock objects
solidBlocks = {}

-- Storage for history cells
historyCells = {}
solidHistoryCells = {}

-- Storage for special zones
ignitionZones = {}
reactionZones = {}
turbulentZones = {}
suppressReconstructionZones = {}
