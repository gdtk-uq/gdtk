-- Module for constructing SolidBlock objects, required by prep.lua.
--
-- Authors: RJG and Kyle D.
--

-- Class for SolidBlock construction
local SolidBlock = {
   myType = "SolidBlock",
} -- end SSolidBlock

function SolidBlock:new(o)
   local flag = type(self)=='table' and self.myType=='SolidBlock'
   if not flag then
      error("Make sure that you are using SolidBlock:new{} and not SolidBlock.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"grid", "initTemperature", "active",
                                "label", "bcList", "properties"})
   if not flag then
      error("Invalid name for item supplied to SolidBlock constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block for later construction of the config file.
   -- Note that we want the solidblock ids to continue on from the fluidblocks.
   -- We are assuming the fluidblocks have already been created
   print("Adding new solid block.")
   o.id = #solidBlocks + #fluidBlocks
   solidBlocks[#solidBlocks+1] = o
   -- Must have a grid and initial temperature
   if not o.grid then
      error("You need to supply a grid to SolidBlock constructor.", 2)
   end
   if (not o.grid.get_type) or o.grid:get_type() ~= "structured_grid" then
      error("You need to supply a structured_grid to SolidBlock constructor.", 2)
   end
   if not o.initTemperature then
      error("You need to supply an initTemperature to SolidBlock constructor.", 2)
   end
   if not o.properties then
      error("You need to supply physical properties for the block.", 2)
   end
   if (type(o.properties) == 'table') then
      local flag2 = checkAllowedNames(o.properties, {"rho", "k", "Cp",
						     "k11", "k12", "k13",
						     "k21", "k22", "k23",
						     "k31", "k32", "k33"})
      if not flag2 then
         error("Invalid name for item supplied in SolidBlock properties table.", 2)
      end
      -- Fill in the k values as 0.0 if not set.
      local kProps = {"k", "k11", "k12", "k13", "k21", "k22", "k23", "k31", "k32", "k33"}
      for _,kName in ipairs(kProps) do
	 o.properties[kName] = o.properties[kName] or 0.0
      end
   end
   -- Fill in some defaults, if not already set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("SOLIDBLOCK-%d", o.id)
   if o.bcList then
      o.bcList = deepclone(o.bcList, false)
   else
      o.bcList =  {} -- boundary conditions
   end
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or SolidAdiabaticBC:new{}
   end
   -- Extract some information from the StructuredGrid
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if config.dimensions == 3 then
      o.nkc = o.grid:get_nkv() - 1
   else
      o.nkc = 1
   end
   o.ncells = o.nic * o.njc * o.nkc
   -- The following table p for the corner locations,
   -- is to be used later for testing for block connections.
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
   return o
end

function SolidBlock:tojson()
   local str = string.format('"solid_block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   -- Boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      if not self.bcList[face].is_solid_domain_bc then
	 local errMsg = string.format("Boundary condition problem for solid block:%d, face:%s\n", self.id, face)
	 errMsg = errMsg .. "       This boundary condition should be a solid domain b.c.\n"
	 errMsg = errMsg .. "       The preparation stage cannot complete successfully.\n"
	 error(errMsg)
      end
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '}'
   return str
end -- SolidBlock:tojson()

local function SolidBlockArray(t)
   -- Expect one table as argument, with named fields.
   -- Returns an array of blocks defined over a single region.
   local flag = checkAllowedNames(t, {"grid", "initTemperature", "active",
				      "label", "bcList", "properties",
				      "nib", "njb", "nkb"})
   if not flag then
      error("Invalid name for item supplied to SolidBlockArray.", 2)
   end
   if not t.grid then
      error("You need to supply a grid to SolidBlockArray.", 2)
   end
   if (not t.grid.get_type) or t.grid:get_type() ~= "structured_grid" then
      error("You need to supply a structured_grid to SolidBlockArray.", 2)
   end
   if not t.initTemperature then
      error("You need to supply an 'initTemperature' to SolidBlockArray.", 2)
   end
   if not t.properties then
      error("You need to supply 'properties' to SolidBlockArray.", 2)
   end
   t.bcList = t.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      t.bcList[face] = t.bcList[face] or SolidAdiabaticBC:new{}
   end
   -- Numbers of subblocks in each coordinate direction
   t.nib = t.nib or 1
   t.njb = t.njb or 1
   t.nkb = t.nkb or 1
   if config.dimensions == 2 then
      t.nkb = 1
   end
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells in the D-domain.
   local nic_total = t.grid:get_niv() - 1
   local dnic = math.floor(nic_total/t.nib)
   local njc_total = t.grid:get_njv() - 1
   local dnjc = math.floor(njc_total/t.njb)
   local nkc_total = t.grid:get_nkv() - 1
   local dnkc = math.floor(nkc_total/t.nkb)
   if config.dimensions == 2 then
      nkc_total = 1
      dnkc = 1
   end
   local blockArray = {} -- will be a multi-dimensional array indexed as [i][j][k]
   local blockCollection = {} -- will be a single-dimensional array
   local nic_remaining = nic_total
   local i0 = 0
   for ib = 1, t.nib do
      blockArray[ib] = {}
      local nic = math.floor(nic_remaining/(t.nib-ib+1))
      if (ib == t.nib) then
	 -- Last block has to pick up remaining cells.
         nic = nic_remaining
      end
      nic_remaining = nic_remaining - nic
      local njc_remaining = njc_total
      local j0 = 0
      for jb = 1, t.njb do
         local njc = math.floor(njc_remaining/(t.njb-jb+1))
	 if (jb == t.njb) then
	    njc = njc_remaining
	 end
         njc_remaining = njc_remaining - njc
	 if config.dimensions == 2 then
	    -- 2D flow
	    local subgrid = t.grid:subgrid(i0,nic+1,j0,njc+1)
	    local bcList = {north=SolidAdiabaticBC:new{}, east=SolidAdiabaticBC:new{},
			    south=SolidAdiabaticBC:new{}, west=SolidAdiabaticBC:new{}}
	    if ib == 1 then
	       bcList['west'] = t.bcList['west']
	    end
	    if ib == t.nib then
	       bcList['east'] = t.bcList['east']
	    end
	    if jb == 1 then
	       bcList['south'] = t.bcList['south']
	    end
	    if jb == t.njb then
	       bcList['north'] = t.bcList['north']
	    end
	    local new_block = SolidBlock:new{grid=subgrid, properties=t.properties,
                                             initTemperature=t.initTemperature,
                                             bcList=bcList}
	    blockArray[ib][jb] = new_block
	    blockCollection[#blockCollection+1] = new_block
	 else
	    -- error("SolidBlockArray not implemented for 3D.")
             -- 3D flow, need one more level in the array
	    blockArray[ib][jb] = {}
            local nkc_remaining = nkc_total
            local k0 = 0
	    for kb = 1, t.nkb do
               local nkc = math.floor(nkc_remaining/(t.nkb-kb+1))
               if (kb == t.nkb) then
                  nkc = nkc_remaining
               end
               nkc_remaining = nkc_remaining - nkc
	       local subgrid = t.grid:subgrid(i0,nic+1,j0,njc+1,k0,nkc+1)
	       local bcList = {north=SolidAdiabaticBC:new{}, east=SolidAdiabaticBC:new{},
			       south=SolidAdiabaticBC:new{}, west=SolidAdiabaticBC:new{},
			       top=SolidAdiabaticBC:new{}, bottom=SolidAdiabaticBC:new{}}
	       if ib == 1 then
		  bcList['west'] = t.bcList['west']
	       end
	       if ib == t.nib then
		  bcList['east'] = t.bcList['east']
	       end
	       if jb == 1 then
		  bcList['south'] = t.bcList['south']
	       end
	       if jb == t.njb then
		  bcList['north'] = t.bcList['north']
	       end
	       if kb == 1 then
		  bcList['bottom'] = t.bcList['bottom']
	       end
	       if kb == t.nkb then
		  bcList['top'] = t.bcList['top']
	       end
	       local new_block = SolidBlock:new{grid=subgrid, properties=t.properties,
                                                initTemperature=t.initTemperature,
                                                bcList=bcList}
	       blockArray[ib][jb][kb] = new_block
	       blockCollection[#blockCollection+1] = new_block
               -- Prepare k0 at end of loop, ready for next iteration
               k0 = k0 + nkc
	    end -- kb loop
	 end -- dimensions
         -- Prepare j0 at end of loop, ready for next iteration
         j0 = j0 + njc
      end -- jb loop
      -- Prepare i0 at end of loop, ready for next iteration
      i0 = i0 + nic
   end -- ib loop
   -- Make the inter-subblock connections
   if #blockCollection > 1 then
      identifyBlockConnections(blockCollection)
   end
   return blockArray
end -- SolidBlockArray

return {
   SolidBlock = SolidBlock,
   SolidBlockArray = SolidBlockArray
}
