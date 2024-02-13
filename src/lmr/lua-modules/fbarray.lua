-- fbarray.lua
-- Module for constructing FluidBlock objects, required by prep.lua.
--
-- Authors: PJ and RJG
--

local gridarray = require 'gridarray'
local GridArray = gridarray.GridArray

-- Class for FluidBlock-Array objects.
local FBArray = {
   myType = "FBArray"
}

function FBArray:new(o)
   local flag = type(self)=='table' and self.myType=='FBArray'
   if not flag then
      error("Make sure that you are using FBArray:new{} and not FBArray.new{}", 2)
   end
   o = o or {}
   local flag = checkAllowedNames(o, {"grid", "gridArray", "initialState", "fillCondition",
				      "active", "label", "omegaz", "may_be_turbulent",
                                      "bcList", "nib", "njb", "nkb", "xforceList"})
   if not flag then
      error("Invalid name for item supplied to FBArray constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- We will embed the FBArray identity in the individual blocks
   -- and we would like that identity to start from 0 for the D code.
   o.id = #(fluidBlockArrays)
   --
   if not o.initialState then
      -- try old name
      o.initialState = o.fillCondition
   end
   if not o.initialState then
      error("You need supply an initialState to FBArray constructor.", 2)
   end
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("FBArray-%d", o.id)
   if fluidBlockArraysDict[o.label] then
      error('Have previously defined a FBArray with label "' .. o.label .. '"', 2)
   end
   fluidBlockArraysDict[o.label] = o.id
   o.omegaz = o.omegaz or 0.0
   if o.may_be_turbulent == nil then
      o.may_be_turbulent = true
   end
   o.bcList = o.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or WallBC_WithSlip:new()
   end
   for _,face in ipairs(faceList(config.dimensions)) do
      if (face ~= "west") and (o.bcList[face].type == "inflow_shock_fitting") then
         error("Shock fitting cannot be used on " .. face .. " face.", 2)
      end
   end
   o.shock_fitting = (o.bcList["west"].type == "inflow_shock_fitting")
   o.xforceList = o.xforceList or {}
   if (not o.grid) and (not o.gridArray) then
      error("You need to supply a grid or gridArray to FBArray constructor.", 2)
   end
   if o.grid then
      -- We will take a single grid and divide it into an array of subgrids.
      if (not o.grid.get_type) or o.grid:get_type() ~= "structured_grid" then
         error("You need to supply a structured_grid to FBArray constructor.", 2)
      end
      -- Numbers of subblocks in each coordinate direction
      o.nib = o.nib or 1
      o.njb = o.njb or 1
      o.nkb = o.nkb or 1
      if config.dimensions == 2 then
         o.nkb = 1
      end
      o.gridArray = GridArray:new{grid=o.grid, nib=o.nib, njb=o.njb, nkb=o.nkb}
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
      o.nib = o.gridArray.nib
      o.njb = o.gridArray.njb
      o.nkb = o.gridArray.nkb
   end
   -- Extract some information from the StructuredGrid
   o.niv = o.gridArray.grid:get_niv()
   o.njv = o.gridArray.grid:get_njv()
   o.nkv = o.gridArray.grid:get_nkv()
   o.nics = o.gridArray.nics
   o.njcs = o.gridArray.njcs
   o.nkcs = o.gridArray.nkcs
   --
   -- At this point, we have an array of grids and the overall grid.
   -- It's time to generate the array of FluidBlocks.
   --
   o.blockArray = {} -- will be a multi-dimensional array, indexed as [ib][jb][kb],
                     -- with 1<=ib<=nib, 1<=jb<=njb, 1<=kb<=nkb
   o.blockCollection = {} -- will be a single-dimensional array, also starting at 1
   local i0 = 0
   for ib = 1, o.nib do
      o.blockArray[ib] = {}
      local j0 = 0
      for jb = 1, o.njb do
	 if config.dimensions == 2 then
	    -- 2D flow
	    local subgrid = o.gridArray.grids[ib][jb]
	    local bcList = {north=WallBC_WithSlip:new(), east=WallBC_WithSlip:new(),
			    south=WallBC_WithSlip:new(), west=WallBC_WithSlip:new()}
	    if ib == 1 then
	       bcList['west'] = o.bcList['west']
	    end
	    if ib == o.nib then
	       bcList['east'] = o.bcList['east']
	    end
	    if jb == 1 then
	       bcList['south'] = o.bcList['south']
	    end
	    if jb == o.njb then
	       bcList['north'] = o.bcList['north']
	    end
            local label = o.label .. string.format("-blk-%d-%d", ib, jb);
	    local new_block = FluidBlock:new{grid=subgrid,
                                             initialState=o.initialState,
                                             active=o.active,
                                             label=label,
                                             omegaz=o.omegaz,
                                             may_be_turbulent=o.may_be_turbulent,
                                             bcList=bcList,
                                             fluidBlockArrayId=o.id}
	    o.blockArray[ib][jb] = new_block
	    o.blockCollection[#o.blockCollection+1] = new_block
	 else
	    -- 3D flow, need one more level in the array
	    o.blockArray[ib][jb] = {}
            local k0 = 0
	    for kb = 1, o.nkb do
	       local subgrid = o.gridArray.grids[ib][jb][kb]
	       local bcList = {north=WallBC_WithSlip:new(), east=WallBC_WithSlip:new(),
			       south=WallBC_WithSlip:new(), west=WallBC_WithSlip:new(),
			       top=WallBC_WithSlip:new(), bottom=WallBC_WithSlip:new()}
	       if ib == 1 then
		  bcList['west'] = o.bcList['west']
	       end
	       if ib == o.nib then
		  bcList['east'] = o.bcList['east']
	       end
	       if jb == 1 then
		  bcList['south'] = o.bcList['south']
	       end
	       if jb == o.njb then
		  bcList['north'] = o.bcList['north']
	       end
	       if kb == 1 then
		  bcList['bottom'] = o.bcList['bottom']
	       end
	       if kb == o.nkb then
		  bcList['top'] = o.bcList['top']
	       end
               local label = o.label .. string.format("-%d-%d-%d", ib, jb, kb);
	       local new_block = FluidBlock:new{grid=subgrid,
                                                initialState=o.initialState,
                                                active=o.active,
                                                label=label,
                                                omegaz=o.omegaz,
                                                may_be_turbulent=o.may_be_turbulent,
                                                bcList=bcList,
                                                fluidBlockArrayId=o.id,}
	       o.blockArray[ib][jb][kb] = new_block
	       o.blockCollection[#o.blockCollection+1] = new_block
	    end -- kb loop
	 end -- dimensions
      end -- jb loop
   end -- ib loop
   -- Make the inter-subblock connections
   if #o.blockCollection > 1 then
      identifyBlockConnections(o.blockCollection)
   end
   --
   -- Retain meta-information about the new FluidBlockArray
   -- for use later in the user-defined functions, during simulation.
   -- Note that the index of this array starts at 1 (in the Lua way).
   fluidBlockArrays[#fluidBlockArrays+1] = o
   --
   if o.shock_fitting then
      -- Prepare the velocity-weights for later writing to file.
      -- Note that vertex indicies start at 0 in each direction.
      o.niv = o.gridArray.grid:get_niv()
      o.njv = o.gridArray.grid:get_njv()
      o.nkv = o.gridArray.grid:get_nkv()
      o.velocity_weights = {}
      for i = 0, o.niv-1 do
         o.velocity_weights[i] = {}
         for j = 0, o.njv-1 do
            o.velocity_weights[i][j] = {}
         end
      end
      for j = 0, o.njv-1 do
         for k = 0, o.nkv-1 do
            distances = {}
            i = o.niv-1
            p0 = o.gridArray.grid:get_vtx(i,j,k)
            distances[i] = 0.0
            for irev = 1, o.niv-1 do
               i = o.niv-1-irev
               p1 = o.gridArray.grid:get_vtx(i,j,k)
               ds = vabs(p1-p0)
               distances[i] = distances[i+1] + ds
               p0 = p1 -- for next step
            end
            local arc_length = distances[0]
            for i = 0, o.niv-1 do
               o.velocity_weights[i][j][k] = distances[i] / arc_length
            end
         end
      end
   else
      o.velocity_weights = nil
   end
   --
   return o
end -- FBArray:new

-- Retain the original behaviour.
local function FluidBlockArray(t)
   print("NOTE: You have called FluidBlockArray{}; prefer FBArray:new{}.")
   o = FBArray:new(t)
   return o.blockArray
end

function FBArray:tojson()
   local str = string.format('"fluid_block_array_%d": {\n', self.id)
   str = str .. string.format('    "nib": %d,\n', self.nib)
   str = str .. string.format('    "njb": %d,\n', self.njb)
   str = str .. string.format('    "nkb": %d,\n', self.nkb)
   str = str .. string.format('    "niv": %d,\n', self.niv)
   str = str .. string.format('    "njv": %d,\n', self.njv)
   str = str .. string.format('    "nkv": %d,\n', self.nkv)
   str = str .. string.format('    "shock_fitting": %s,\n', tostring(self.shock_fitting))
   --
   str = str .. string.format('    "nics": [ ')
   for i=1,#(self.nics)-1 do
      str = str .. string.format('%d, ', self.nics[i])
   end
   str = str .. string.format('%d ],\n', self.nics[#self.nics])
   --
   str = str .. string.format('    "njcs": [ ')
   for i=1,#(self.njcs)-1 do
      str = str .. string.format('%d, ', self.njcs[i])
   end
   str = str .. string.format('%d ],\n', self.njcs[#self.njcs])
   --
   str = str .. string.format('    "nkcs": [ ')
   for i=1,#(self.nkcs)-1 do
      str = str .. string.format('%d, ', self.nkcs[i])
   end
   str = str .. string.format('%d ],\n', self.nkcs[#self.nkcs])
   --
   str = str .. string.format('    "blockIds": [ ')
   for ib=1,#(self.blockCollection)-1 do
      str = str .. string.format('%d, ', self.blockCollection[ib].id)
   end
   str = str .. string.format('%d ]\n', self.blockCollection[#self.blockCollection].id)
   str = str .. '}'
   return str
end

return {
   FBArray = FBArray,
   FluidBlockArray = FluidBlockArray
}
