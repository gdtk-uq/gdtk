-- prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- This script is read by the Eilmer4 program at preparation time,
-- before reading and processing the user's input script.
--
-- Authors: PJ, RJG, Kyle D., Nick G. and Daryl B.
--
print("Loading prep.lua...")

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

require 'configoptions'
config = configoptions.config

require 'gridpro'
-- Make these functions global so that users may call them
-- in the configuration script
applyGridproConnectivity = gridpro.applyGridproConnectivity
applyGridproBoundaryConditions = gridpro.applyGridproBoundaryConditions
to_eilmer_axis_map = gridpro.to_eilmer_axis_map

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
SuppressViscousStressesZone = zones.SuppressViscousStressesZone

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
suppressViscousStressesZones = {}

-- ---------------------------------------------------------------------------------------

function build_config_files(job)
   perform_spatial_gradient_consistency_check()
   warn_if_blocks_not_connected()
   if config.do_temporal_DFT then
       check_DFT_settings()
   end
   print("Build config files for job:", job)
   os.execute("mkdir -p config")
   write_config_file("config/" .. job .. ".config")
   write_control_file("config/" .. job .. ".control")
   write_times_file("config/" .. job .. ".times")
   write_block_list_file("config/" .. job .. ".list")
   write_mpimap_file("config/" .. job .. ".mpimap")
   write_fluidBlockArrays_file("config/" .. job .. ".fluidBlockArrays")
   if config.grid_motion == "shock_fitting" then
      write_shock_fitting_helper_files(job)
   end
   print("Done building config files.")
end

function build_grid_and_flow_files(job)
   if #fluidBlockIdsForPrep == 0 then
      -- We'll set *all* blocks for processing.
      for i=1,#fluidBlocks do
         fluidBlockIdsForPrep[i] = fluidBlocks[i].id
      end
   end

   if (config.new_flow_format) then
      os.execute("mkdir -p CellData/field/t0000") -- this violoates "Don't Repeat Yourself" (CellData/field)
   else
      os.execute("mkdir -p flow/t0000")
   end

   os.execute("mkdir -p grid/t0000")
   if #solidBlocks >= 1 then
      os.execute("mkdir -p solid-grid/t0000")
      os.execute("mkdir -p solid/t0000")
   end
   for i, id in ipairs(fluidBlockIdsForPrep) do
      if false then
         -- May activate print statement for debug.
         print("FluidBlock id=", id)
      end
      local idx = id+1
      local fileName = "grid/t0000/" .. job .. string.format(".grid.b%04d.t0000", id)
      if (config.grid_format == "gziptext") then
	 fluidBlocks[idx].grid:write_to_gzip_file(fileName .. ".gz")
      elseif (config.grid_format == "rawbinary") then
	 fluidBlocks[idx].grid:write_to_raw_binary_file(fileName .. ".bin")
      else
	 error(string.format("Oops, invalid grid_format: %s", config.grid_format))
      end

      local fileName;

      if (config.new_flow_format) then 
         fileName = "CellData/field/t0000/" .. job .. string.format(".field.b%04d.t0000", id)
         if ((config.flow_format == "eilmer4text") or (config.flow_format == "eilmer4binary")) then
            fileName = fileName .. ".zip"
         else
            error(string.format("Oops, new flow format selected, %s is not valid", config.flow_format))
         end
      else
         fileName = "flow/t0000/" .. job .. string.format(".flow.b%04d.t0000", id)
         if (config.flow_format == "gziptext") then
            fileName = fileName .. ".gz"
         elseif (config.flow_format == "rawbinary") then
            fileName = fileName .. ".bin"
         else
            error(string.format("Oops, invalid flow_format: %s", config.flow_format))
         end
      end
      local ifs = fluidBlocks[idx].initialState
      if type(ifs) == "table" and ifs.myType == "FlowState" then
	 -- We have one of the pure-Lua FlowState objects and we convert it to
	 -- a wrapped-D-language _FlowState object.
	 ifs = _FlowState:new(ifs)
      elseif type(ifs) == "function" then
	 -- leave alone
      elseif type(ifs) == "userdata" then
	 -- presume to be a wrapped-D-language _FlowState object already
      elseif type(ifs) == "string" then
         -- We are given the name of a flow file and we'll copy that in place directly
         existingFlowFile = ifs
         -- Presume file exists and let 'cp' command complain if it doesn't
         cmd = "cp " .. existingFlowFile .. " " .. fileName
         returnCode = os.execute(cmd)
         if returnCode ~= 0 then
            errMsg = "Error while trying to copy an existing flow solution as initial flow solution.\n"
            errMsg = errMsg .. "FluidBlock id= " .. id .. "\n"
            errMsg = errMsg .. "Specified existing flow file: " .. existingFlowFile .. "\n"
            errMsg = errMsg .. "Check this file exists and is readable.\n"
            errMsg = errMsg .. "Bailing out!\n"
            error(errMsg)
         end
         -- Otherwise succesful.
         str = string.format("Initialised FluidBlock id= %d with existing flow solution: \n\t%s", id, existingFlowFile)
         print(str)
      else
	 error("Unexpected type for initial flow state in block.")
      end
      if type(ifs) ~= "string" then
         local grid = fluidBlocks[idx].grid
         local omegaz = fluidBlocks[idx].omegaz
         if grid:get_type() == "structured_grid" then
            write_initial_sg_flow_file(fileName, grid, ifs, config.start_time, omegaz)
         else
            write_initial_usg_flow_file(fileName, grid, ifs, config.start_time, omegaz)
         end
      end
   end
   for i = 1, #solidBlocks do
      local id = solidBlocks[i].id
      print("SolidBlock id=", id)
      local fileName = "solid-grid/t0000/" .. job .. string.format(".solid-grid.b%04d.t0000.gz", id)
      solidBlocks[i].grid:write_to_gzip_file(fileName)
      local fileName = "solid/t0000/" .. job .. string.format(".solid.b%04d.t0000", id)
      writeInitialSolidFile(fileName, solidBlocks[i].grid,
			    solidBlocks[i].initTemperature, solidBlocks[i].properties, config.start_time)
      os.execute("gzip -f " .. fileName)
   end
   --
   if #fluidBlocks == 0 then print("Warning: number of FluidBlocks is zero.") end
   print("Done building grid and flow files.")
end

print("Done loading prep.lua")
