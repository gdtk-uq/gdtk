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

-- sleep() function copied from http://lua-users.org/wiki/SleepFunction
local clock = os.clock
function sleep(n)  -- seconds
   -- warning: clock can eventually wrap around for sufficiently large n
   -- (whose value is platform dependent).  Even for n == 1, clock() - t0
   -- might become negative on the second that clock wraps.
   local t0 = clock()
   while clock() - t0 <= n do end
end

function checkAllowedNames(myTable, allowedNames)
   local setOfNames = {}
   local namesOk = true
   for i,name in ipairs(allowedNames) do
      setOfNames[name] = true
   end
   for k,v in pairs(myTable) do
      if not setOfNames[k] then
	 print("Warning: Invalid name: ", k)
	 namesOk = false
      end
   end
   return namesOk
end

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

require 'gridpro'
-- Make these functions global so that users may call them
-- in the configuration script
applyGridproConnectivity = gridpro.applyGridproConnectivity
applyGridproBoundaryConditions = gridpro.applyGridproBoundaryConditions

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

-- The following module is not required by prep itself, but we load it here
-- to make it available in the user's script.
require 'billig'

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

-- ---------------------------------------------------------------------------------------

function initTurbulence(fs, turbulence_model_name)
    -- Setup dynamic turbulent primitive array and check user inputs
    if turbulence_model_name == "none" then
        -- Check to ensure the user hasn't tried defining any turbulence stuff
        if fs.tke~=0.0 then error(string.format("Turbulence model is none but tke set: %.18e", fs.tke)) end
        if fs.omega~=1.0 then error(string.format("Turbulence model is none but omega set: %.18e", fs.omega)) end
        if fs.nuhat~=0.0 then error(string.format("Turbulence model is none but nuhat set: %.18e", fs.nuhat)) end
        turb = {0.0, 1.0}
    elseif turbulence_model_name == "k_omega" then
        if fs.nuhat~=0.0 then error(string.format("Turbulence model is k_omega but nuhat set: %.18e", fs.nuhat)) end
        turb = {fs.tke, fs.omega}
    elseif string.find(turbulence_model_name , "spalart_allmaras") then
        if fs.tke~=0.0 then error(string.format("Turbulence model is spalart_allmaras but tke set: %.18e", fs.tke)) end
        if fs.omega~=1.0 then error(string.format("Turbulence model is spalart_allmaras but omega set: %.18e", fs.omega)) end
        turb = {fs.nuhat, 0.0}
    else
        error(string.format("Unsupported turbulence model: %s", turbulence_model_name))
    end
    return turb
end

function to_eilmer_axis_map(gridpro_ijk)
   -- Convert from GridPro axis_map string to Eilmer3 axis_map string.
   -- From GridPro manual, Section 7.3.2 Connectivity Information.
   -- Example, 123 --> '+i+j+k'
   local axis_map = {[0]='xx', [1]='+i', [2]='+j', [3]='+k',
		     [4]='-i', [5]='-j', [6]='-k'}
   if type(gridpro_ijk) == "number" then
      gridpro_ijk = string.format("%03d", gridpro_ijk)
   end
   if type(gridpro_ijk) ~= "string" then
      error("Expected a string or integer of three digits but got:"..tostring(gridpro_ijk))
   end
   local eilmer_ijk = axis_map[tonumber(string.sub(gridpro_ijk, 1, 1))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 2, 2))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 3, 3))]
   return eilmer_ijk
end

function checkCellVolumes(t)
   if not t then
      t = {}
   end
   -- Stop reporting cells after this limit
   local badCellLimit = t.badCellLimit or 20
   local badCells = {}
   local badCellCount = 0
   for ib,blk in ipairs(fluidBlocks) do
      local grid = blk.grid
      for idx=0,grid:get_ncells()-1 do
	 local vol = grid:cellVolume(idx)
	 if vol <= 0 then
	    badCellCount = badCellCount + 1
	    badCells[#badCells+1] = {ib, idx}
	    if badCellCount >= badCellLimit then
	       return false, badCells
	    end
	 end
      end
   end
   if #badCells > 0 then
      return false, badCells
   else
      return true, badCells
   end
end

function perform_spatial_gradient_consistency_check()
   -- Not all spatial gradient options are available, depending on the type of grid.
   -- First, search for any unstructured grids, since these are the most restricted.
   unstructuredGridsPresent = false
   for _,blk in ipairs(fluidBlocks) do
      if blk.grid:get_type() == "unstructured_grid" then
	 unstructuredGridsPresent = true
	 break
      end
   end
   if unstructuredGridsPresent then
      if config.spatial_deriv_calc == "divergence" then
	 print("NOTE: config.spatial_deriv_calc is being set to 'least_squares' because unstructured grids detected.")
	 config.spatial_deriv_calc = "least_squares"
      end
      if config.spatial_deriv_locn == "vertices" then
	 print("NOTE: config.spatial_deriv_locn is being set to 'cells' when using least squares.")
	 config.spatial_deriv_locn = "cells"
      end
   else
      -- Only structured grids are present.
      -- 2D structured grids have all options available.
      if config.dimensions == 3 then
         if config.spatial_deriv_calc == "divergence" then
            print("NOTE: config.spatial_deriv_calc is being set to 'least_squares' for 3D simulations.")
            config.spatial_deriv_calc = "least_squares"
         end
      end
   end
   if not config.spatial_deriv_from_many_points then
      if config.spatial_deriv_locn == "vertices" then
         print("NOTE: config.spatial_deriv_location is being set to 'faces' for 2-point derivatives.")
         config.spatial_deriv_locn = "faces"
      end
   end
end

function warn_if_blocks_not_connected()
   -- It would be very unusual to have defined multiple FluidBlocks
   -- and not have any connections between them.
   -- Such an arrangement would be unintended, almost certainly.
   if #fluidBlocks > 1 then
      local n = 0
      for _,blk in ipairs(fluidBlocks) do
         for bndry,bc in pairs(blk.bcList) do
            if string.find(bc.type, "exchange_") then
               n = n + 1
            end
         end
      end
      if n == 0 then
         print("WARNING: Did not find any inter-block connections.")
         print("         Did you forget to call identifyBlockConnections()?")
      end
   end
end

function check_DFT_settings()
    -- Check to see that the DFT has been correctly configured.

    if not ((config.DFT_n_modes * config.DFT_step_interval) == config.max_step) then
        print("WARNING: config.DFT_n_modes * config.DFT_step_interval should equal config.max_step")
    end
    if not config.fixed_time_step then
        print("WARNING: should turn config.fixed_time_step on when computing the DFT")
    end
end

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

function build_block_files(job)
   if #fluidBlocksForPrep == 0 then
      -- We'll set *all* blocks for processing.
      for i=1,#fluidBlocks do
         fluidBlocksForPrep[i] = fluidBlocks[i].id
      end
   end
   os.execute("mkdir -p grid/t0000")
   os.execute("mkdir -p flow/t0000")
   if #solidBlocks >= 1 then
      os.execute("mkdir -p solid-grid/t0000")
      os.execute("mkdir -p solid/t0000")
   end
   for i, id in ipairs(fluidBlocksForPrep) do
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
      local fileName = "flow/t0000/" .. job .. string.format(".flow.b%04d.t0000", id)
      if (config.flow_format == "gziptext") then
	 fileName = fileName .. ".gz"
      elseif (config.flow_format == "rawbinary") then
	 fileName = fileName .. ".bin"
      else
	 error(string.format("Oops, invalid flow_format: %s", config.flow_format))
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
   print("Done building block files.")
end

print("Done loading prep.lua")
