-- prep-steady.lua
--
-- A place to put helper functions and classes for flow preparation activities for the steady-state solver.
--
-- Authors: RJG, KAD, NNG, PJ
--

require 'lua_helper'
require 'blk_conn'
require 'bc'

-- Extract directory and file names from central config.
local lmrconfig = require 'lmrconfig'
local lmrCfg = lmrconfig.lmrCfg

local globalconfig = require 'globalconfig'
config = globalconfig.config

local nkconfig = require 'nkconfig'
NewtonKrylovGlobalConfig = nkconfig.NewtonKrylovGlobalConfig
nkPhases = nkconfig.NewtonKrylovPhases
NewtonKrylovPhase = nkconfig.NewtonKrylovPhase

local grid = require 'grid'
RegisteredGrid = grid.RegisteredGrid
connectGrids = grid.connectGrids
connectionAsJSON = grid.connectionAsJSON

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

local solidblock = require 'solidblock'
SolidBlock = solidblock.SolidBlock

local solidthermalmodel = require 'solidthermalmodel'
ConstantPropertiesModel = solidthermalmodel.ConstantPropertiesModel
registerSolidModels = solidthermalmodel.registerSolidModels

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

local prep_check = require 'prep_check'
initTurbulence = prep_check.initTurbulence
checkCellVolumes = prep_check.checkCellVolumes
perform_spatial_gradient_consistency_check = prep_check.perform_spatial_gradient_consistency_check
warn_if_blocks_not_connected = prep_check.warn_if_blocks_not_connected
check_DFT_settings = prep_check.check_DFT_settings

json = require 'json'

-- Users might want to put instructions for grid and flow preparation in one script.
-- When this is the case, we do NOT want to do anything with registering grids
-- since this should have been done earlier at the prep-grid stage.
-- So, we disable any actions related to those function calls.
function registerFluidGrid()
   if verbosity >= 1 then
      print("    registerFluidGrid(): Do NOTHING when in prep-sim mode.")
   end
end

function registerFluidGridArray()
   if verbosity >= 1 then
      print("    registerFluidGridArray(): Do NOTHING when in prep-sim mode.")
   end
end

function registerSolidGrid()
   if verbosity >= 1 then
      print("    registerSolidGrid(): Do NOTHING when in prep-sim mode.")
   end
end

function registerSolidGridArray()
   if verbosity >= 1 then
      print("    registerSolidGridArray(): Do NOTHING when in prep-sim mode.")
   end
end

function identifyGridConnections()
   if verbosity >= 1 then
      print("    identifyGridConnections(): Do NOTHING when in prep-sim mode.")
   end
end

function importGridproConnectivity()
   if verbosity >= 1 then
      print("    importGridproConnectivity(): Do NOTHING when in prep-sim mode.")
   end
end

function importGridproBCs()
   if verbosity >= 1 then
      print("    importGridproBCs(): Do NOTHING when in prep-sim mode.")
   end
end

-- ---------------------------------------------------------------------------------------

-- Storage for later definitions of FluidBlock objects.
-- Note that the index for this array starts at 1, in the Lua way.
-- The block ids start at 0 to be like the indexing inside the D code.
-- Yes, this is confusing...
gridsList = {}
fluidBlocks = {}
-- Storage for metadata on FluidBlockArrays, which correspond to gridArrays
-- as defined in prepgrids.lua.
fluidBlockArrays = {}
-- We may want to look up the blocks via labels rather than numerical id
-- in user-defined procedures.
-- The following dictionaries store the connections.
fluidBlocksDict = {}

-- The user may assign the MPI task id for eack block manually
-- but, if they don't, a default distribution will be made.
mpiTasks = nil

-- Storgage for later definitions of SolidBlock objects
solidBlocks = {}
-- Storage for solid models
_solidModels = {}

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
   -- NOTE: those fluid blocks that are connected to solid blocks won't be completely
   --       configured at the end of this function because we have not yet built solid
   --       blocks. Fluid blocks with solid block connections will be completely
   --       configured when the solid blocks are built.
   if false then -- debug
      print("Make FluidBlock objects.")
   end
   -- Build a fluid grids list to work on.
   fluidGridsList = {}
   for idx, g in ipairs(gridsList) do
      if g.fieldType == "fluid" then
         fluidGridsList[#fluidGridsList+1] = g
      end
   end
   for idx, g in ipairs(fluidGridsList) do
      if false then -- debug
         print("Make FluidBlock for g.id=", g.id)
      end
      local ifs
      if g.fsTag then ifs = flowDict[g.fsTag] end
      if not ifs then
         local msg = string.format("Grid.id=%d fsTag=%s, does not have a valid initial FlowState.\n", g.id, g.fsTag)
         msg = msg .. "Keys in flowDict are: "
         for k,v in pairs(flowDict) do
            msg = msg .. string.format(" %s", k)
         end
         error(msg)
      end
      if g.type == "structured_grid" then
         -- Build the bc list for this block,
         -- using the tags that the user applied when building the grid.
         bcs = {}
         for _,face in ipairs(faceList(config.dimensions)) do
            local tag = g.bcTags[face] -- get the user's spec for this face
            if tag and tag ~= "" then
               local bc = bcDict[tag]
               if bc then
                  bcs[face] = bc
               else
                  print(string.format("WARNING: For Grid.id=%d face=%s there is no bcDict entry for bc tag=%s", g.id, face, tag))
                  local msg = "    Keys in bcDict are:"
                  for k,v in pairs(bcDict) do
                     msg = msg .. string.format(" %s", k)
                  end
                  print(msg)
               end
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
      fb = FluidBlock:new{gridMetadata=g, initialState=ifs, bcDict=bcs, active=g.active}
   end
   if false then -- debug
      print("Make block-to-block connections.")
   end
   -- We only have fluid blocks created presently, so we can't handle grid connections
   -- that go to solid blocks. The test will be on ids. If the grid id is greater
   -- then the available blocks, we will pass on that connection.
   for idx, c in ipairs(gridConnections) do
      if (c.idA < #fluidBlocks and c.idB < #fluidBlocks) then
         -- Remember that the Lua array index will be one more than the block id.
         connectBlocks(fluidBlocks[c.idA+1], c.faceA, fluidBlocks[c.idB+1], c.faceB, c.orientation)
      end
   end
end

-- ---------------------------------------------------------------------------------------

function makeSolidBlocks(bcDict, temperatureDict)
   if false then -- debug
      print("Make SolidBlock objects.")
   end
   -- Build a list of solid grids to work on.
   solidGridsList = {}
   for idx, g in ipairs(gridsList) do
      if g.fieldType == "solid" then
         solidGridsList[#solidGridsList+1] = g
      end
   end
   for idx, g in ipairs(solidGridsList) do
      if false then -- debug
         print("Make SolidBlock for g.id=", g.id)
      end

      local initT
      if g.ssTag then initT = temperatureDict[g.ssTag] end
      if not initT then
         error(string.format("Grid.id=%d ssTag='%s': does not seem to have a valid initial temperature state.", g.id, g.ssTag))
      end

      if g.type == "structured_grid" then
         -- Build the bc list for this block,
         -- using the tags that the user applied when building the grid.
         bcs = {}
         for _,face in ipairs(faceList(config.dimensions)) do
            local tag = g.solidBCTags[face] -- get the user's spec for this face
            if tag and tag ~= "" then
               local bc = bcDict[tag]
               if bc then bcs[face] = bc end
            end
         end
      else -- unstructured_grid
         error("Unstructured grids not available for solid domains.")
      end
      if false then --debug
         print("  bcs=[")
         for face, bc in pairs(bcs) do print("    face=", face, " bc=", bc) end
         print("  ]")
      end
      sb = SolidBlock:new{gridMetadata=g, initTemperature=initT, modelTag=g.solidModelTag, bcList=bcs}
   end
   if false then -- debug
      print("Make block-to-block connections.")
   end
   -- We need to handle the various cases of solid and fluid block connections.
   nfb = #fluidBlocks
   for idx, c in ipairs(gridConnections) do
      if (c.idA >= nfb and c.idB >= nfb) then
         -- We have both solid blocks.
         sidA = c.idA - nfb
         sidB = c.idB - nfb
         -- Remember that the Lua array index will be one more than the block id.
         connectBlocks(solidBlocks[sidA+1], c.faceA, solidBlocks[sidB+1], c.faceB, c.orientation)
      elseif (c.idA < nfb and c.idB >= nfb) then
         -- A: fluid; B: solid
         sidB = c.idB - nfb
         connectBlocks(fluidBlocks[c.idA+1], c.faceA, solidBlocks[sidB+1], c.faceB, c.orientation)
      elseif (c.idA >= nfb and c.idB < nfb) then
         -- A: soilid; B: fluid
         sidA = c.idA - nfb
         connectBlocks(solidBlocks[sidA+1], c.faceA, fluidBlocks[c.idB+1], c.faceB, c.orientation)
      end
      -- Do NOT handle case when both blocks are fluid.
      -- This hase been done already. Check end of function makeFluidBlocks()
   end
end

-- ---------------------------------------------------------------------------------------

function readGridMetadata()
   print('Read Grid Metadata.')
   local fileName = lmrconfig.gridMetadataFilename()
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
   local gridArrayMetadata = jsonData["gridarrays"]
   if gridArrayMetadata then
      -- Only assign the data if it is not nil, else the initial empty table will be retained.
      fluidBlockArrays = gridArrayMetadata
      print(string.format('  #fluidBlockArrays: %d', #fluidBlockArrays))
      -- for i,fba in ipairs(fluidBlockArrays) do print(string.format("fba[%d]=%s", i, json.stringify(fba))) end
   end
   --
   local ngrids = jsonData["ngrids"]
   for i=1, ngrids do
      fileName = lmrconfig.gridMetadataFilename(i-1)
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

function buildRuntimeConfigFiles()
   print('Build runtime config files.')
   perform_spatial_gradient_consistency_check()
   warn_if_blocks_not_connected()
   if config.do_temporal_DFT then
       check_DFT_settings()
   end
   local cfgDir = lmrconfig.lmrCfg["config-directory"]
   os.execute("mkdir -p " .. cfgDir)
   write_config_file(lmrconfig.simulationConfigFilename())
   if (config.solver_mode ~= "steady") then
      write_control_file(cfgDir .. "/" .. lmrCfg["control-filename"])
   end
   write_block_list_file(lmrconfig.blockListFilename())
   write_mpimap_file(lmrconfig.mpimapFilename())
   write_fluidBlockArrays_file(cfgDir .. "/" .. lmrCfg["fluid-block-arrays-filename"])
   if (config.solver_mode == "steady") then
      nkconfig.setIgnoreFlagInPhases(nkPhases)
      nkconfig.writeNKConfigToFile(NewtonKrylovGlobalConfig, nkPhases, lmrconfig.nkConfigFilename())
   end

   if false then -- debug
      print("Done buildRuntimeConfigFiles.")
   end
end

function buildGridAndFieldFiles()
   print('Build fluid files.')
   -- 1. Write metadata file (this might be redundant in parallel, but shouldn't hurt)
   local snapshotTopLvl = lmrconfig.snapshotDirectory()
   os.execute("mkdir -p " .. snapshotTopLvl)
   writeFluidMetadata()

   -- 2. Now look for the work to do per block
   if #fluidBlockIdsForPrep == 0 then
      -- We'll set *all* blocks for processing.
      for i=1,#fluidBlocks do
         fluidBlockIdsForPrep[i] = fluidBlocks[i].id
      end
   end

   local initDir = lmrCfg["initial-field-directory"]
   local snapshotInitDir = lmrconfig.snapshotDirectory(initDir)
   local fileName
   os.execute("mkdir -p " .. snapshotInitDir)
   for i, id in ipairs(fluidBlockIdsForPrep) do
      local idx = id+1
      local blk = fluidBlocks[idx]
      local gridMetadata = gridsList[idx]
      -- We assume a direct match between FluidBlock and grid id numbers.
      local gridFilename = lmrconfig.gridFilename(id)
      local gridForSimFilename = lmrconfig.gridForSimFilename(lmrCfg["initial-field-directory"], id)
      if gridMetadata.type == "structured_grid" then
         blk.grid = StructuredGrid:new{filename=gridFilename, fmt=config.grid_format}
      else
         blk.grid = UnstructuredGrid:new{filename=gridFilename, fmt=config.grid_format}
      end
      -- Copy from staging area to initial snapshot area
      os.execute("cp " .. gridFilename .. " " .. gridForSimFilename)
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
      elseif type(ifs) == "userdata" then
         -- presume to be a wrapped-D-language _FlowState object already
      else
         error("Unexpected type for initial flow state in block.")
      end
      -- Ready to use initialState in the Dlang function.
      local grid = fluidBlocks[idx].grid
      local omegaz = fluidBlocks[idx].omegaz
      writeInitialFluidFile(id, grid, ifs, nil, omegaz)
   end
   --
   if #fluidBlocks == 0 then print("Warning: number of FluidBlocks is zero.") end
   --
   if #solidBlocks > 0 then
      print("Build solid files.")
      writeSolidMetadata()
   end
   for i=1,#solidBlocks do
      local blk = solidBlocks[i]
      local id = blk.id
      local gridMetadata = gridsList[id+1]
      local gridFilename = lmrconfig.gridFilename(id)
      local gridForSimFilename = lmrconfig.gridForSimFilename(lmrCfg["initial-field-directory"], id)
      if gridMetadata.type == "structured_grid" then
         blk.grid = StructuredGrid:new{filename=gridFilename, fmt=config.grid_format}
      else
         error("Unstructured grids not available for solid blocks.")
      end
      -- Copy from staging area to initial snapshot area
      os.execute("cp " .. gridFilename .. " " .. gridForSimFilename)
      --
      writeInitialSolidFile(id, blk.grid, blk.initTemperature, _solidModels[blk.modelTag])
   end
end

if false then -- debug
   print("Done loading prep-flow.lua")
end


