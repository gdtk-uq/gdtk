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

local configoptions = require 'configoptions'
config = configoptions.config


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
function writeGridFiles(jobName)
   print(string.format('Write grid files for job="%s"', jobName))
   --
   os.execute("mkdir -p grid")
   local fileName = "grid/" .. jobName .. ".grid-metadata"
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
   os.execute("mkdir -p grid/t0000")
   for i, g in ipairs(gridsList) do
      if false then -- May activate print statement for debug.
         print("grid id=", g.id)
      end
      -- Write the grid proper.
      local fileName = "grid/t0000/" .. jobName .. string.format(".grid.b%04d.t0000", g.id)
      if config.grid_format == "gziptext" then
	 g.grid:write_to_gzip_file(fileName .. ".gz")
      elseif config.grid_format == "rawbinary" then
	 g.grid:write_to_raw_binary_file(fileName .. ".bin")
      else
	 error(string.format("Oops, invalid grid_format: %s", config.grid_format))
      end
      -- Write the grid metadata.
      local fileName = "grid/" .. jobName .. string.format(".grid.b%04d.metadata", g.id)
      local f = assert(io.open(fileName, "w"))
      f:write(g:tojson() .. '\n')
      f:close()
   end
   print(string.format("  #grids %d", #gridsList))
   print(string.format("  #gridArrays %d", #gridArraysList))
end
