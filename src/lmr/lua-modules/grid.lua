-- grid.lua
--
-- Each Grid object, in the Lua domain, will hold a reference to a Dlang
-- StructuredGrid or UnstructuredGrid, plus some metadata needed to assign
-- boundary conditions and the like.
--
-- PJ, 2021-10-05 pulled out of prep-grids.lua
--

local Grid = {
   myType = "Grid"
}

function Grid:new(o)
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
   local flag = type(self)=='table' and self.myType=='Grid'
   if not flag then
      error("Make sure that you are using Grid:new{} and not Grid.new{}", 2)
   end
   local flag = type(o)=='table'
   if not flag then
      error("Grid constructor expects a single table with named items.", 2)
   end
   flag = checkAllowedNames(o, {"grid", "tag", "fsTag", "bcTags", "gridArrayId"})
   if not flag then
      error("Invalid name for item supplied to Grid constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new object, for later use.
   -- Note that we want id to start at zero for the D code.
   o.id = #gridsList
   gridsList[#gridsList+1] = o
   o.tag = o.tag or string.format("Grid-%d", o.id)
   if gridsDict[o.tag] then
      error('Have previously defined a Grid with tag "' .. o.tag .. '"', 2)
   end
   gridsDict[o.tag] = o.id
   -- Set to -1 if NOT part of a grid-array, otherwise use supplied value
   o.gridArrayId = o.gridArrayId or -1
   -- Initial FlowState tag
   o.fsTag = o.fsTag or ""
   -- Must have a grid.
   assert(o.grid, "need to supply a grid")
   -- Check the grid information.
   if config.dimensions ~= o.grid:get_dimensions() then
      local msg = string.format("Mismatch in dimensions, config %d grid %d.",
				config.dimensions, o.grid.get_dimensions())
      error(msg)
   end
   o.bcTags = o.bcTags or {}
   o.type = o.grid:get_type()
   if o.type == "structured_grid" then
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
	 o.bcTags[face] = o.bcTags[face] or o.grid:get_tag(face)
      end
   end
   if o.type == "unstructured_grid" then
      -- Extract some information from the UnstructuredGrid
      o.ncells = o.grid:get_ncells()
      o.nvertices = o.grid:get_nvertices()
      o.nfaces = o.grid:get_nfaces()
      o.nboundaries = o.grid:get_nboundaries()
      -- Attach boundary conditions from list or from the dictionary of conditions.
      for i = 0, o.nboundaries-1 do
         o.bcTags[i] = o.bcTags[i] or o.grid:get_boundaryset_tag(i)
      end
   end
   return o
end -- Grid:new

function Grid:tojson()
   str = '{\n'
   str = str .. string.format('  "tag": "%s",\n', self.tag)
   str = str .. string.format('  "fsTag": "%s",\n', self.fsTag)
   str = str .. string.format('  "type": "%s",\n', self.type)
   if self.type == "structured_grid" then
      str = str .. string.format('  "dimensions": %d,\n', self.grid:get_dimensions())
      str = str .. string.format('  "niv": %d,\n', self.grid:get_niv())
      str = str .. string.format('  "njv": %d,\n', self.grid:get_njv())
      str = str .. string.format('  "nkv": %d,\n', self.grid:get_nkv())
      str = str .. string.format('  "nic": %d,\n', self.nic)
      str = str .. string.format('  "njc": %d,\n', self.njc)
      str = str .. string.format('  "nkc": %d,\n', self.nkc)
      local fmt = '  "p%d": {"x":%.18e, "y":%.18e, "z":%.18e},\n'
      for i=0, 3 do
         str = str .. string.format(fmt, i, self.p[i].x, self.p[i].y, self.p[i].z)
      end
      if config.dimensions == 3 then
         for i=4, 7 do
            str = str .. string.format(fmt, i, self.p[i].x, self.p[i].y, self.p[i].z)
         end
      end
   else -- unstructured-grid
      str = str .. string.format('  "dimensions": %d,\n', self.grid:get_dimensions())
      str = str .. string.format('  "nvertices": %d,\n', self.grid:get_nvertices())
      str = str .. string.format('  "ncells": %d,\n', self.grid:get_ncells())
      str = str .. string.format('  "nfaces": %d,\n', self.grid:get_nfaces())
      str = str .. string.format('  "nboundaries": %d,\n', self.grid:get_nboundaries())
   end
   str = str .. '  "bcTags": {\n'
   if self.type == "structured_grid" then
      -- Expect named boundaries
      for k, v in pairs(self.bcTags) do
         str = str .. string.format('    "%s": "%s",\n', k, v)
      end
   else -- unstructured_grid
      -- Expect numbered boundary sets.
      -- Note that Dlang numbering starts at zero.
      for j=0, self.grid:get_nboundaries()-1 do
         str = str .. string.format('    "%d": "%s",\n', j, self.bcTags[j])
      end
   end
   str = str .. '    "dummy_entry_without_trailing_comma": "xxxx"\n'
   str = str .. '  },\n'
   str = str .. string.format('  "gridArrayId": %d\n', self.gridArrayId) -- last item, no comma
   str = str .. '}\n'
   return str
end -- end Grid:tojson()

-------------------------------------------------------------------------
--
-- Structured grids may be connected full face to full face.
--
-- needs storage: connectionList = {}

local function connectGrids(idA, faceA, idB, faceB, orientation)
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

local function connectionAsJSON(c)
   str = string.format('{"idA": %d, "faceA": "%s", "idB": %d, "faceB": "%s", "orientation": %d}',
                       c.idA, c.faceA, c.idB, c.faceB, c.orientation)
   return str
end


local function identifyGridConnections(includeList, excludeList, tolerance)
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
      for _,A in ipairs(includeList) do
         if A.grid:get_type() == "structured_grid" then
            myGridList[#myGridList+1] = A
         end
      end
   else
      -- The caller has not provided a list; use the global grids list.
      for _,A in ipairs(gridsList) do
         if A.grid:get_type() == "structured_grid" then
            myGridList[#myGridList+1] = A
         end
      end
   end
   excludeList = excludeList or {}
   -- Put unstructured grid objects into the exclude list because they don't
   -- have a simple topology that can always be matched to a structured grid.
   for _,A in ipairs(myGridList) do
      if A.grid:get_type() == "unstructured_grid" then excludeList[#excludeList+1] = A end
   end
   if false then
      print('DEBUG myGridList=[')
      for _,A in ipairs(myGridList) do print('  ', A.id, ',') end
      print(']')
   end
   if #myGridList == 0 then
      -- There are no structured grids that are to be (potentially) connected.
      return
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
		     connectGrids(A.id, faceA, B.id, faceB, orientation)
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

return {
   Grid = Grid,
   connectGrids = connectGrids,
   connectionAsJSON = connectionAsJSON,
   identifyGridConnections = identifyGridConnections
}
