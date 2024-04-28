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

local globalconfig = require 'globalconfig'
config = globalconfig.config

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

local gridproimport = require 'gridproimport'
importGridproConnectivity = gridproimport.importGridproConnectivity
importGridproBCs = gridproimport.importGridproBCs

local grid = require 'grid'
RegisteredGrid = grid.RegisteredGrid
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

local solidblock = require 'solidblock'
SolidBlock = solidblock.SolidBlock
SolidBlockArray = solidblock.SolidBlockArray

local solidthermalmodel = require 'solidthermalmodel'
ConstantPropertiesModel = solidthermalmodel.ConstantPropertiesModel
registerSolidModels = solidthermalmodel.registerSolidModels

local mpi = require 'mpi'
mpiDistributeBlocks = mpi.mpiDistributeBlocks
mpiDistributeFBArray = mpi.mpiDistributeFBArray

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
function makeSolidBlocks()
   if verbosity >= 1 then
      print("    makeSolidBlocks(): Do NOTHING when in prep-grid mode.")
   end
end
function identifyBlockConnections()
   if verbosity >= 1 then
      print("    identifyBlockConnections(): Do NOTHING when in prep-grid mode.")
   end
end
function setHistoryPoint()
   if verbosity >= 1 then
      print("    setHistoryPoint(): Do NOTHING when in prep-grid mode.")
   end
end
function mpiDistributeBlocks()
   if verbosity >= 1 then
      print("    mpiDistributeBlocks(): Do NOTHING when in prep-grid mode.")
   end
end
-- We also need to define the storage for special zones, otherwise prep-grid complains that they don't exist
ignitionZones = {}
reactionZones = {}
turbulentZones = {}
suppressReconstructionZones = {}
suppressViscousStressesZones = {}
_solidModels = {}

-------------------------------------------------------------------------

gridsList = {} -- to hold Grid objects
gridsDict = {}
gridArraysList = {} -- to hold GridArray objects
connectionList = {}


function registerFluidGrid(o)
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
   o.fieldType = "fluid"
   local rgrid = RegisteredGrid:new(o)
   return rgrid.id
end -- function registerFluidGrid

function registerFluidGridArray(o)
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
end -- registerFluidGridArray


function registerSolidGrid(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid or UnstructuredGrid object that has been generated
   --    or imported.
   -- tag: a string to identify the grid later in the user's script
   -- ssTag: a string that will be used to select the initial solid condition from
   --    a dictionary when the SolidBlock is later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlock is later constructed.
   -- modelTag: a string that will be used to select the solid thermal model
   --    from a dictionary when the SolidBlock is later constructed.
   --
   -- Returns:
   -- the grid id number so that we may assign it and use it when making connections.
   --

   -- This function does some translation of parameter names.
   -- The user provides what seems sensible withing a registerSolidGrid function.
   -- We re-pack the table here with the names expected by the back-end RegisteredGrid object.
   t = {}
   t.fieldType = "solid"
   t.grid = o.grid
   t.tag = o.tag
   t.ssTag = o.ssTag
   t.solidModelTag = o.modelTag
   t.solidBCTags = o.bcTags
   local rgrid = RegisteredGrid:new(t)
   return rgrid.id
end -- function registerSolidGrid

function registerSolidGridArray(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid or UnstructuredGrid object that has been generated
   --    or imported.
   -- tag: a string to identify the grid later in the user's script
   -- ssTag: a string that will be used to select the initial solid condition from
   --    a dictionary when the SolidBlock is later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlock is later constructed.
   -- modelTag: a string that will be used to select the solid thermal model
   --    from a dictionary when the SolidBlock is later constructed.
   --
   -- Returns:
   -- the id of GridArray object so that the user may use the interior pieces later in their script.
   --
   -- This function does some translation of parameter names.
   -- The user provides what seems sensible withing a registerSolidGridArray function.
   t = {}
   t.fieldType = "solid"
   t.grid = o.grid
   t.tag = o.tag
   t.ssTag = o.ssTag
   t.solidModelTag = o.modelTag
   t.solidBCTags = o.bcTags
   t.nib = o.nib or 1
   t.njb = o.njb or 1
   t.nkb = o.nkb or 1
   local rga = GridArray:new(t)
   return rga.id
end -- registerFluidGridArray
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
   --
   local any = false
   for i = 1, #gridArraysList do
      local ga = gridArraysList[i]
      if ga.shock_fitting then any = true end
   end
   if any then
      print("  For shock-fitting, write rails and weights files.")
   end
   for i = 1, #gridArraysList do
      local ga = gridArraysList[i]
      if ga.shock_fitting then
         local filename = string.format(lmrconfig.gridDirectory() .. "/gridarray-%04d.rails", ga.id)
         local f = assert(io.open(filename, "w"))
         f:write("# Rails are presently described by the initial west- and east-boundary coordinates.\n")
         for k = 0, ga.nkv-1 do
            for j = 0, ga.njv-1 do
               local pw = ga.grid:get_vtx(0,j,k)
               local pe = ga.grid:get_vtx(ga.niv-1,j,k)
               f:write(string.format("%.18e %.18e %.18e %.18e %.18e %.18e\n",
                                     pw.x, pw.y, pw.z, pe.x, pe.y, pe.z))
            end
         end
         f:close()
         local filename = string.format(lmrconfig.gridDirectory() .. "/gridarray-%04d.weights", ga.id)
         local f = assert(io.open(filename, "w"))
         f:write("# Weights represent the arc-length distance of each vertex from the east-boundary vertex.\n")
         for k = 0, ga.nkv-1 do
            for j = 0, ga.njv-1 do
               for i = 0, ga.niv-1 do
                  f:write(string.format("%.18e\n", ga.velocity_weights[i][j][k]))
               end
            end
         end
         f:close()
      end
   end
   --
   print(string.format("  #grids %d", #gridsList))
   print(string.format("  #gridArrays %d", #gridArraysList))
end
