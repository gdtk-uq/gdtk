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

require 'configoptions'
config = configoptions.config

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
to_eilmer_axis_map = gridpro.to_eilmer_axis_map

require 'grid'
Grid = grid.Grid

require 'gridarray'
GridArray = gridarray.GridArray

-------------------------------------------------------------------------

gridsList = {} -- to hold RegisteredGrid objects
gridsDict = {}


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

--
-- Structured grids may be connected full face to full face.
--
connectionList = {}

function connectGrids(idA, faceA, idB, faceB, orientation)
   if false then -- debug
      print(string.format('connectGrids(idA=%d, faceA="%s", idB=%d, faceB="%s", orientation=%d)',
                          idA, faceA, idB, faceB, orientation))
   end
   local gridA = gridsList[idA+1] -- Note that the id values start at zero.
   local gridB = gridsList[idB+1]
   if gridA.grid:get_type() ~= "structured_grid" or gridB.grid:get_type() ~= "structured_grid" then
      error("connectGrids() Works only for structured grids.", 2)
   end
   connectionList[#connectionList+1] = {idA=idA, faceA=faceA, idB=idB, faceB=faceB, orientation=orientation}
end

function connectionAsJSON(c)
   str = string.format('{"idA": %d, "faceA": "%s", "idB": %d, "faceB": "%s", "orientation": %d}',
                       c.idA, c.faceA, c.idB, c.faceB, c.orientation)
   return str
end


function identifyGridConnections(includeList, excludeList, tolerance)
   -- Identify grid connections by trying to match corner points.
   -- Parameters (all optional):
   -- includeList: the list of structured grid objects to be included in the search.
   --    If nil, the whole collection is searched.
   -- excludeList: list of pairs of structured grid objects that should not be
   --    included in the search for connections.
   -- tolerance: spatial tolerance for the colocation of vertices
   --
   local myGridList = {}
   if includeList then
      -- The caller has provided a list of grids to bound the search.
      for _,v in ipairs(includeList) do myGridList[#myGridList+1] = v end
   else
      -- The caller has not provided a list; use the global grids list.
      for _,v in ipairs(gridsList) do myGridList[#myGridList+1] = v end
   end
   excludeList = excludeList or {}
   -- Put unstructured grid objects into the exclude list because they don't
   -- have a simple topology that can always be matched to a structured grid.
   for _,A in ipairs(myGridList) do
      if A.grid:get_type() == "unstructured_grid" then excludeList[#excludeList+1] = A end
   end
   if false then -- debug
      print('myGridList=[')
      for _,A in ipairs(myGridList) do print('  ', A.id, ',') end
      print(']')
   end
   tolerance = tolerance or 1.0e-6
   --
   for _,A in ipairs(myGridList) do
      for _,B in ipairs(myGridList) do
	 if (A ~= B) and (not isPairInList({A, B}, excludeList)) then
	    -- print("Proceed with test for coincident vertices.") -- DEBUG
	    local connectionCount = 0
	    if config.dimensions == 2 then
	       -- print("2D test A.id=", A.id, " B.id=", B.id) -- DEBUG
	       for vtxPairs,connection in pairs(connections2D) do
                  if false then -- debug
                     print("vtxPairs=", tostringVtxPairList(vtxPairs),
                           "connection=", tostringConnection(connection))
                  end
                  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
		     local faceA, faceB, orientation = unpack(connection)
		     connectGrids(A.id, faceA, B.id, faceB, 0)
		     connectionCount = connectionCount + 1
		  end
	       end
	    else
	       -- print("   3D test")
               for vtxPairs,connection in pairs(connections3D) do
		  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
		     local faceA, faceB, orientation = unpack(connection)
		     connectBlocks(A.id, faceA, B.id, faceB, orientation)
		     connectionCount = connectionCount + 1
		  end
	       end
	    end
	    if connectionCount > 0 then
	       -- So we don't double-up on connections.
	       excludeList[#excludeList+1] = {A,B}
	    end
	 end -- if (A ~= B...
      end -- for _,B
   end -- for _,A
end -- identifyGridConnections


--
-- A GridArray will be a single, structured grid split into a number of sub-grids.
-- GridArray objects hold a GridArray object and some metadata that is needed
-- for the staged preparation.
--
gridArraysList = {} -- to hold GridArray objects

RegisteredGridArray = {
   myType = "RegisteredGridArray"
}

function RegisteredGridArray:new(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid object that has been generated or imported.
   -- gridArray: a GridArray object or array of StructuredGrid objects
   -- tag: a string to identify the gridArray later in the user's script
   -- fsTag: a string that will be used to select the initial flow condition from
   --    a dictionary when the FluidBlocks are later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlocks are later constructed.
   --
   local flag = type(self)=='table' and self.myType=='RegisteredGridArray'
   if not flag then
      error("Make sure that you are using RegisteredGridArray:new{} and not RegisteredGridArray.new{}", 2)
   end
   local flag = type(o)=='table'
   if not flag then
      error("RegisterGridArray expected to receive a single table with named entries", 2)
   end
   local flag = checkAllowedNames(o, {"grid", "gridArray", "tag", "fsTag", "bcTags", "nib", "njb", "nkb"})
   if not flag then
      error("Invalid name for item supplied to RegisteredGridArray constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new object, for later use.
   -- Note that we want id to start at zero for the D code.
   o.id = #gridArraysList
   gridArraysList[#gridArraysList+1] = o
   --
   o.tag = o.tag or ""
   o.fsTag = o.fsTag or ""
   o.bcTags = o.bcTags or {} -- for boundary conditions to be applied to the FluidBlocks
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcTags[face] = o.bcTags[face] or ""
   end
   if (not o.grid) and (not o.gridArray) then
      error("You need to supply a grid or gridArray to registergridArray.", 2)
   end
   if o.grid then
      -- We will take a single grid and divide it into an array of subgrids.
      if (not o.grid.get_type) or o.grid:get_type() ~= "structured_grid" then
         error("You need to supply a structured_grid to registerGridArray.", 2)
      end
      -- Numbers of subblocks in each coordinate direction
      local nib = o.nib or 1
      local njb = o.njb or 1
      local nkb = o.nkb or 1
      if config.dimensions == 2 then
         nkb = 1
      end
      o.gridArray = GridArray:new{grid=o.grid, nib=nib, njb=njb, nkb=nkb}
   else
      -- We were not given a single grid,
      -- so we assume that we were given the array of subgrids.
      -- Join these into a single overall grid.
      if o.gridArray.myType and o.gridArray.myType == "GridArray" then
         -- The GridArray object is already constructed, so we leave it as is.
      else
         if type(o.gridArray) == "table" then
            o.gridArray = GridArray:new{gridArray=o.gridArray}
         else
            error("gridArray should be an array of grid objects.", 2)
         end
      end
   end
   --
   -- At this point, we have an array of grids and the overall grid.
   local gridCollection = {}
   for ib = 1, o.nib do
      for jb = 1, o.njb do
	 if config.dimensions == 2 then
	    -- 2D flow
	    local subgrid = o.gridArray.grids[ib][jb]
	    local bcTags = {north="", east="", south="", west=""}
	    if ib == 1 then bcTags['west'] = o.bcTags['west'] end
	    if ib == o.gridArray.nib then bcTags['east'] = o.bcTags['east'] end
	    if jb == 1 then bcTags['south'] = o.bcTags['south'] end
	    if jb == o.gridArray.njb then bcTags['north'] = o.bcTags['north'] end
	    local rgrid = Grid:new{grid=subgrid, fsTag=o.fsTag, bcTags=bcTags, gridArrayId=id}
	    gridCollection[#gridCollection+1] = rgrid
	 else
	    -- 3D flow, need one more level in the array
	    for kb = 1, o.nkb do
	       local subgrid = o.gridArray.grids[ib][jb][kb]
	       local bcTags = {north="", east="", south="", west="", top="", bottom=""}
	       if ib == 1 then bcTags['west'] = o.bcTags['west'] end
	       if ib == o.gridArray.nib then bcTags['east'] = o.bcTags['east'] end
	       if jb == 1 then bcTags['south'] = o.bcTags['south'] end
	       if jb == o.gridArray.njb then bcTags['north'] = o.bcTags['north'] end
	       if kb == 1 then bcTags['bottom'] = o.bcTags['bottom'] end
	       if kb == o.gridArray.nkb then bcTags['top'] = o.bctags['top'] end
	       local rgrid = Grid:new{grid=subgrid, fsTag=o.fsTag, bcTags=bcTags, gridArrayId=id}
	       gridCollection[#gridCollection+1] = rgrid
	    end -- kb loop
	 end -- dimensions
      end -- jb loop
   end -- ib loop
   -- Make the inter-subblock connections
   if #gridCollection > 1 then
      identifyGridConnections(gridCollection)
   end
   --
   return o
end -- RegisteredGridArray:new()

function RegisteredGridArray:tojson()
   -- [TODO] PJ 2021-10-04 finish this function.
   str = '{\n'
   str = str .. string.format('  "tag": "%s",\n', self.tag)
   str = str .. string.format('  "fsTag": "%s",\n', self.fsTag)
   str = str .. string.format('  "type": "%s",\n', self.myType)
   str = str .. '  ' .. self.gridArray:tojson()
   str = str .. '}\n'
   return str
end -- end RegisteredGridArray:tojson()


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
   -- the id of RegisteredGridArray object so that the user may use the interior pieces later in their script.
   local rga = RegisteredGridArray:new(o)
   return rga.id
end -- registerGridArray

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
   --
   -- [TODO] PJ 2021-10-04, write RegisteredGridArrays data.
end
