-- prep-grids.lua
-- A place to put helper functions and classes for grid-preparation activities.
-- This script is read by the Eilmer4 program at grid-preparation time,
-- before reading and processing the user's input script.
--
-- Authors: PJ, RJG, Kyle D., Nick G. and Daryl B.
--
if false then -- debug
   print("Begin loading prep-grids.lua.")
end

require 'lua_helper'
require 'blk_conn'
require 'bc'

-- Extract grid directory and names from central config.
local lmrconfig = require 'lmrconfig'
local lmrCfg = lmrconfig.lmrCfg
local gridDir = lmrCfg["grid-directory"]
local gridMD = lmrCfg["grid-metadata-filename"]
local blkIdxFmt = lmrCfg["block-index-format"]

local configoptions = require 'configoptions'
config = configoptions.config

local nkconfig = require 'nkconfig'
NewtonKrylovGlobalConfig = nkconfig.NewtonKrylovGlobalConfig
nkPhases = nkconfig.NewtonKrylovPhases
NewtonKrylovPhase = nkconfig.NewtonKrylovPhase

local gridpro = require 'gridpro'
-- Make these functions global so that users may call them
-- in the configuration script
applyGridproConnectivity = gridpro.applyGridproConnectivity
applyGridproBoundaryConditions = gridpro.applyGridproBoundaryConditions
to_eilmer_axis_map = gridpro.to_eilmer_axis_map

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

local prep_check = require 'prep_check'
initTurbulence = prep_check.initTurbulence
checkCellVolumes = prep_check.checkCellVolumes
perform_spatial_gradient_consistency_check = prep_check.perform_spatial_gradient_consistency_check
warn_if_blocks_not_connected = prep_check.warn_if_blocks_not_connected
check_DFT_settings = prep_check.check_DFT_settings


-- Users might want to put instructions for grid and flow preparation in one script.
-- When this is the case, we do NOT want to do anything with makeFluidBlocks when
-- we are in a grid-only prep mode.
-- So, we disable any action related to that function call.
function makeFluidBlocks()
   if verbosity >= 1 then
      print("    makeFluidBlocks(): Do NOTHING when in prep-grid mode.")
   end
end

-------------------------------------------------------------------------

gridsList = {} -- to hold Grid objects
gridsDict = {}
gridArraysList = {} -- to hold GridArray objects
connectionList = {}


function registerGrid(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid or UnstructuredGrid object that has been generated
   --    or imported.
   -- tag: a string to identify the grid later in the user's script
   -- fsTag: a string that will be used to select the initial flow condition from
   --    a dictionary when the FluidBlock is later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlock is later constructed.
   -- gridArrayId: needs to be supplied only if the grid is part of a larger array.
   --
   -- Returns:
   -- the grid id number so that we may assign it and use it when making connections.
   --
   local rgrid = Grid:new(o)
   return rgrid.id
end -- function registerGrid

function registerGridArray(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid object that has been generated or imported.
   -- tag: a string to identify the gridArray later in the user's script
   -- fsTag: a string that will be used to select the initial flow condition from
   --    a dictionary when the FluidBlocks are later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlocks are later constructed.
   --
   -- Returns:
   -- the id of GridArray object so that the user may use the interior pieces later in their script.
   local rga = GridArray:new(o)
   return rga.id
end -- registerGridArray


-------------------------------------------------------------------------
--
-- IO functions to write the grid and connection files.
--
function writeGridFiles()
   os.execute("mkdir -p " .. lmrconfig.gridDirectory())
   local fileName = lmrconfig.gridMetadataFilename()
   local f = assert(io.open(fileName, "w"))
   f:write('{\n')
   f:write(string.format('  "ngrids": %d,\n', #gridsList))
   f:write(string.format('  "ngridarrays": %d,\n', #gridArraysList))
   if #gridArraysList > 0 then
      f:write('  "gridarrays": [')
      for i,ga in ipairs(gridArraysList) do
         f:write('    ' .. ga:tojson())
         if i < #gridArraysList then f:write(',') end
         f:write('\n')
      end
      f:write('  ],\n')
   end
   f:write('  "grid-connections": [\n')
   for i, c in ipairs(connectionList) do
      f:write('    ' .. connectionAsJSON(c))
      if i < #connectionList then f:write(',\n') else f:write('\n') end
   end
   f:write('  ]\n') -- Note last item has no comma.
   f:write('}\n')
   f:close()
   print(string.format("  #connections: %d", #connectionList))
   --
   for i, g in ipairs(gridsList) do
      if false then -- May activate print statement for debug.
         print("grid id=", g.id)
      end
      -- Write the grid proper.
      if config.grid_format == "gziptext" then
	 fileName = lmrconfig.gridFilename(g.id, lmrCfg["gzip-extension"])
       	 g.grid:write_to_gzip_file(fileName)
      elseif config.grid_format == "rawbinary" then
	 fileName = lmrconfig.gridFilename(g.id)
	 g.grid:write_to_raw_binary_file(fileName)
      else
	 error(string.format("Oops, invalid grid_format: %s", config.grid_format))
      end
      -- Write the grid metadata.
      fileName = lmrconfig.gridMetadataFilename(g.id)
      local f = assert(io.open(fileName, "w"))
      f:write(g:tojson() .. '\n')
      f:close()
   end
   print(string.format("  #grids %d", #gridsList))
   print(string.format("  #gridArrays %d", #gridArraysList))
end
