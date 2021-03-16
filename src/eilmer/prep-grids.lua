-- prep-grids.lua
-- A place to put helper functions and classes for grid-preparation activities.
-- This script is read by the Eilmer4 program at grid-preparation time,
-- before reading and processing the user's input script.
--
-- Authors: PJ, RJG, Kyle D., Nick G. and Daryl B.
--
print("Loading prep-grids.lua...")

require 'lua_helper'
deepclone = lua_helper.deepclone
checkAllowedNames = lua_helper.checkAllowedNames

require 'blk_conn'
-- Let's pull the symbols out of the blk_conn module
-- and make them global in this namespace
for k,v in pairs(blk_conn) do
   _G[k] = v
end

require 'gridpro'
-- Make these functions global so that users may call them
-- in the configuration script
applyGridproConnectivity = gridpro.applyGridproConnectivity
applyGridproBoundaryConditions = gridpro.applyGridproBoundaryConditions

-------------------------------------------------------------------------

gridsList = {}
gridsDict = {}

--
-- We will store a table of information on each grid that is registered.
--
function registerGrid(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid or UnstructuredGrid object that has been generated
   --    or imported.
   -- tag: a string that will be used to select the initial flow condition from
   --    a dictionary when the FluidBlock is later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlock is later constructed.
   -- GridArrayId: needs to be supplied only if the grid is part of a larger array.
   --
   -- Returns:
   -- the grid id number so that we may assign it and use it when making connections.
   --
   local flag = type(o)=='table'
   if not flag then
      error("registerGrid expects a single table with named items.", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"grid", "tag", "bcTags", "GridArrayId"})
   if not flag then
      error("Invalid name for item supplied to registerGrid.", 2)
   end
   -- Make a record of the new block, for later use.
   -- Note that we want id to start at zero for the D code.
   o.id = #(gridsList)
   gridsList[#(gridsList)+1] = o
   o.tag = o.tag or string.format("Grid-%d", o.id)
   if gridsDict[o.tag] then
      error('Have previously defined a Grid with tag "' .. o.tag .. '"', 2)
   end
   gridsDict[o.tag] = o.id
   -- Set to -1 if NOT part of a grid-array, otherwise use supplied value
   o.GridArrayId = o.GridArrayId or -1
   -- Must have a grid.
   assert(o.grid, "need to supply a grid")
   -- Check the grid information.
   if config.dimensions ~= o.grid:get_dimensions() then
      local msg = string.format("Mismatch in dimensions, config %d grid %d.",
				config.dimensions, o.grid.get_dimensions())
      error(msg)
   end
   if o.grid:get_type() == "structured_grid" then
      -- Extract some information from the StructuredGrid
      -- Note 0-based indexing for vertices and cells.
      o.nic = o.grid:get_niv() - 1
      o.njc = o.grid:get_njv() - 1
      if config.dimensions == 3 then
	 o.nkc = o.grid:get_nkv() - 1
      else
	 o.nkc = 1
      end
      o.ncells = o.nic * o.njc * o.nkc
      -- The following table p for the corner locations,
      -- is to be used later for testing for grid connections.
      o.p = {}
      if config.dimensions == 3 then
	 o.p[0] = o.grid:get_vtx(0, 0, 0)
	 o.p[1] = o.grid:get_vtx(o.nic, 0, 0)
	 o.p[2] = o.grid:get_vtx(o.nic, o.njc, 0)
	 o.p[3] = o.grid:get_vtx(0, o.njc, 0)
	 o.p[4] = o.grid:get_vtx(0, 0, o.nkc)
	 o.p[5] = o.grid:get_vtx(o.nic, 0, o.nkc)
	 o.p[6] = o.grid:get_vtx(o.nic, o.njc, o.nkc)
	 o.p[7] = o.grid:get_vtx(0, o.njc, o.nkc)
      else
	 o.p[0] = o.grid:get_vtx(0, 0)
	 o.p[1] = o.grid:get_vtx(o.nic, 0)
	 o.p[2] = o.grid:get_vtx(o.nic, o.njc)
	 o.p[3] = o.grid:get_vtx(0, o.njc)
      end
      --[[print("Grid id=", o.id, "p0=", tostring(o.p[0]), "p1=", tostring(o.p[1]),
         "p2=", tostring(o.p[2]), "p3=", tostring(o.p[3])) ]]
      -- Attach default boundary conditions for those not specified.
      for _,face in ipairs(faceList(config.dimensions)) do
	 o.bcTags[face] = o.bcTags[face] or ""
      end
   end
   if o.grid:get_type() == "unstructured_grid" then
      -- Extract some information from the UnstructuredGrid
      o.ncells = o.grid:get_ncells()
      o.nvertices = o.grid:get_nvertices()
      o.nfaces = o.grid:get_nfaces()
      o.nboundaries = o.grid:get_nboundaries()
      -- Attach boundary conditions from list or from the dictionary of conditions.
      for i = 0, o.nboundaries-1 do
         o.bcTags[i] = o.grid:get_boundaryset_tag(i)
      end
   end
   return o.id
end -- function registerGrid


--
-- Structured grids may be connected full face to full face.
--
connectionList = {}

function connectGrids(idA, faceA, idB, faceB, orientation)
   print("[FIX-ME] implement connectGrids()")
   local gridA = gridList[idA]
   local gridB = gridList[idB]
   if gridA.grid:get_type() ~= "structured_grid" or gridB.grid:get_type() ~= "structured_grid" then
      error("connectGrids() Works only for structured grids.", 2)
   end
   connectionList[#connectionList+1] = {idA=idA, faceA=faceA, idB=idB, faceB=faceB, orientation=orientation}
end


function identifyGridConnections(t)
   print("[FIX-ME] implement identifyGridConnections()")
end


--
-- A GridArray will be a single, structured grid split into a number of sub-grids.
--
gridArraysList = {}

function registerGridArray(o)
   print("[FIX-ME] implement registerGridArray as per FBArray:new.")
end


--
-- IO functions to write the grid and connection files.
--
function writeGridFiles(jobName)
   print("[FIX-ME] write grid and connection files for job=", jobName)
end
