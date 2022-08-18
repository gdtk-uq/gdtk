-- prep-flow.lua
-- A place to put helper functions and classes for flow preparation activities.
-- This script is read by the Eilmer4 program at config and flow preparation time,
-- before reading and processing the user's input script.
--
-- Authors: PJ, RJG, Kyle D., Nick G. and Daryl B.
--
if false then -- debug
   print("Loading prep-flow.lua.")
end

require 'lua_helper'
require 'blk_conn'
require 'bc'

local configoptions = require 'configoptions'
config = configoptions.config

local grid = require 'grid'
Grid = grid.Grid
connectGrids = grid.connectGrids
connectionAsJSON = grid.connectionAsJSON
identifyGridConnections = grid.identifyGridConnections

local gridarray = require 'gridarray'
GridArray = gridarray.GridArray

local flowstate = require 'flowstate'
FlowState = flowstate.FlowState
makeFlowStateFn = flowstate.makeFlowStateFn

local fluidblock = require 'fluidblock'
FluidBlock = fluidblock.FluidBlock
SBlock2UBlock = fluidblock.SBlock2UBlock
connectBlocks = fluidblock.connectBlocks
identifyBlockConnections = fluidblock.identifyBlockConnections

local fbarray = require 'fbarray'
FBArray = fbarray.FBArray
FluidBlockArray = fbarray.FluidBlockArray

local solidblock = require 'solidblock'
SolidBlock = solidblock.SolidBlock
SolidBlockArray = solidblock.SolidBlockArray

local mpi = require 'mpi'
mpiDistributeBlocks = mpi.mpiDistributeBlocks
mpiDistributeFBArray = mpi.mpiDistributeFBArray

local history = require 'history'
setHistoryPoint = history.setHistoryPoint
setSolidHistoryPoint = history.setSolidHistoryPoint

local zones = require 'zones'
ReactionZone = zones.ReactionZone
IgnitionZone = zones.IgnitionZone
TurbulentZone = zones.TurbulentZone
SuppressReconstructionZone = zones.SuppressReconstructionZone
SuppressViscousStressesZone = zones.SuppressViscousStressesZone

local output = require 'output'
write_control_file = output.write_control_file
write_config_file = output.write_config_file
write_times_file = output.write_times_file
write_block_list_file = output.write_block_list_file
write_mpimap_file = output.write_mpimap_file
write_fluidBlockArrays_file = output.write_fluidBlockArrays_file
write_shock_fitting_helper_files = output.write_shock_fitting_helper_files

require 'sssoptions'

local prep_check = require 'prep_check'
initTurbulence = prep_check.initTurbulence
checkCellVolumes = prep_check.checkCellVolumes
perform_spatial_gradient_consistency_check = prep_check.perform_spatial_gradient_consistency_check
warn_if_blocks_not_connected = prep_check.warn_if_blocks_not_connected
check_DFT_settings = prep_check.check_DFT_settings

json = require 'json'

-- ---------------------------------------------------------------------------------------

-- Storage for later definitions of FluidBlock objects.
-- Note that the index for this array starts at 1, in the Lua way.
-- The block ids start at 0 to be like the indexing inside the D code.
-- Yes, this is confusing...
gridsList = {}
fluidBlocks = {}
-- Storage for later definitions of GridArray and FBArray objects.
gridArraysList = {}
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

gridConnections = {} -- Will be overwritten when the JSON data is parsed.

-- ---------------------------------------------------------------------------------------

function makeFluidBlocks(bcDict, flowDict)
   -- Using the list of grids, apply boundary-conditions and initial flow conditions
   -- to build FluidBlock objects that are ready for flow simulation.
   if false then -- debug
      print("Make FluidBlock objects.")
   end
   for idx, g in ipairs(gridsList) do
      if false then -- debug
         print("Make FluidBlock for g.id=", g.id)
      end
      local ifs
      if g.fsTag then ifs = flowDict[g.fsTag] end
      if not ifs then
         error(string.format("Grid.id=%d fsTag=%s, does not seem to have a valid initial FlowState.", g.id, g.fsTag))
      end
      if g.type == "structured_grid" then
         -- Build the bc list for this block,
         -- using the tags that the user applied when building the grid.
         bcs = {}
         for _,face in ipairs(faceList(config.dimensions)) do
            local tag = g.bcTags[face] -- get the user's spec for this face
            if tag and tag ~= "" then
               local bc = bcDict[tag]
               if bc then bcs[face] = bc end
            end
         end
      else -- unstructured_grid
         -- We just use the same user-supplied dictionary for all blocks.
         bcs = bcDict
      end
      if false then --debug
         print("  bcs=[")
         for face, bc in pairs(bcs) do print("    face=", face, " bc=", bc) end
         print("  ]")
      end
      fb = FluidBlock:new{gridMetadata=g, initialState=ifs, bcDict=bcs}
   end
   if false then -- debug
      print("Make block-to-block connections.")
   end
   for idx, c in ipairs(gridConnections) do
      -- Remember that the Lua array index will be one more than the block id.
      connectBlocks(fluidBlocks[c.idA+1], c.faceA, fluidBlocks[c.idB+1], c.faceB, c.orientation)
   end
end

-- ---------------------------------------------------------------------------------------

function readGridMetadata(jobName)
   print('Read Grid Metadata.')
   local fileName = "grid/" .. jobName .. ".grid-metadata"
   local f = assert(io.open(fileName, "r"))
   local jsonStr = f:read("*a")
   f:close()
   local jsonData = json.parse(jsonStr)
   gridConnections = jsonData["grid-connections"]
   if false then -- debug
      print('number of connections=', #gridConnections)
      for i, c in ipairs(gridConnections) do
         print("i=", i, "idA=", c.idA, "faceA=", c.faceA,
               "idB=", c.idB, "faceB=", c.faceB, "orientation=", c.orientation)
      end
   end
   print(string.format('  #connections: %d', #gridConnections))
   --
   local ngrids = jsonData["ngrids"]
   for i=1, ngrids do
      local fileName = "grid/" .. jobName .. string.format(".grid.b%04d.metadata", i-1)
      if false then --debug
         print('Set up grid object from file', fileName)
      end
      local f = assert(io.open(fileName, "r"))
      local jsonStr = f:read("*a")
      f:close()
      local gridMetadata = json.parse(jsonStr)
      gridMetadata.id = i-1
      gridsList[#gridsList+1] = gridMetadata
   end
   print(string.format('  #grids: %d', #gridsList))
end

function buildRuntimeConfigFiles(jobName)
   print(string.format('Build runtime config files for job="%s"', jobName))
   perform_spatial_gradient_consistency_check()
   warn_if_blocks_not_connected()
   if config.do_temporal_DFT then
       check_DFT_settings()
   end
   os.execute("mkdir -p config")
   write_config_file("config/" .. jobName .. ".config")
   write_control_file("config/" .. jobName .. ".control")
   write_times_file("config/" .. jobName .. ".times")
   write_block_list_file("config/" .. jobName .. ".list")
   write_mpimap_file("config/" .. jobName .. ".mpimap")
   write_fluidBlockArrays_file("config/" .. jobName .. ".fluidBlockArrays")
   if config.grid_motion == "shock_fitting" then
      write_shock_fitting_helper_files(jobName)
   end
   if false then -- debug
      print("Done buildRuntimeConfigFiles.")
   end
end

function buildFlowFiles(jobName)
   print(string.format('Build flow files for job="%s"', jobName))
   if #fluidBlockIdsForPrep == 0 then
      -- We'll set *all* blocks for processing.
      for i=1,#fluidBlocks do
         fluidBlockIdsForPrep[i] = fluidBlocks[i].id
      end
   end
   if (config.new_flow_format) then
      os.execute("mkdir -p CellData/field/t0000")
   else
      os.execute("mkdir -p flow/t0000")
   end
   for i, id in ipairs(fluidBlockIdsForPrep) do
      if false then -- May activate print statement for debug.
         print("FluidBlock id=", id)
      end
      local idx = id+1
      local blk = fluidBlocks[idx]
      local gridMetadata = gridsList[idx]
      local fileName;
      if (config.new_flow_format) then
          fileName = "CellData/field/t0000/" .. jobName .. string.format(".field.b%04d.t0000", id)
          if ((config.flow_format == "eilmer4text") or (config.flow_format == "eilmer4binary")) then
              fileName = fileName .. ".zip"
          else
              error(string.format("Oops, new flow format selected, %s is not valid", config.flow_format))
          end
      else
          fileName = "flow/t0000/" .. jobName .. string.format(".flow.b%04d.t0000", id)
          if config.flow_format == "gziptext" then
         fileName = fileName .. ".gz"
          elseif config.flow_format == "rawbinary" then
         fileName = fileName .. ".bin"
          else
         error(string.format("Oops, invalid flow_format: %s", config.flow_format))
          end
      end
      --
      if not blk.grid then
         -- We assume a direct match between FluidBlock and grid id numbers.
         local gridFileName = "grid/t0000/" .. jobName .. string.format(".grid.b%04d.t0000", id)
         if config.grid_format == "gziptext" then
            gridFileName = gridFileName .. ".gz"
         elseif config.grid_format == "rawbinary" then
            gridFileName = gridFileName .. ".bin"
         else
            error(string.format("Oops, invalid grid_format: %s", config.grid_format))
         end
         if gridMetadata.type == "structured_grid" then
            blk.grid = StructuredGrid:new{filename=gridFileName, fmt="gziptext"}
         else
            blk.grid = UnstructuredGrid:new{filename=gridFileName, fmt="gziptext"}
         end
      end
      --
      local ifs = blk.initialState
      -- Check if we need to do something special with initialState
      -- before calling the Dlang function to write the initial flow file.
      if type(ifs) == "table" and ifs.myType == "FlowState" then
	 -- We have one of the pure-Lua FlowState objects and we convert it to
	 -- a wrapped-D-language _FlowState object.
	 ifs = _FlowState.new(ifs)
      elseif type(ifs) == "function" then
	 -- leave alone
      elseif type(ifs) == "number" then
	 -- presume that we have a D-language _FlowState struct already in place
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
         str = string.format("Initialised FluidBlock id= %d with existing flow solution: \n\t%s",
                             id, existingFlowFile)
         print(str)
      else
	 error("Unexpected type for initial flow state in block.")
      end
      -- Ready to use initialState in the Dlang function.
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
   --
   if #fluidBlocks == 0 then print("Warning: number of FluidBlocks is zero.") end
   if false then -- debug
      print("Done buildingFlowFiles.")
   end
end

if false then -- debug
   print("Done loading prep-flow.lua")
end
