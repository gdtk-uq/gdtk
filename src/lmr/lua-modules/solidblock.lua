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
   flag = checkAllowedNames(o, {"grid", "gridMetadata", "initTemperature", "active",
                                "label", "bcList", "modelTag"})
   if not flag then
      error("Invalid name for item supplied to SolidBlock constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block for later construction of the config file.
   -- Note that we want the solidblock ids to continue on from the fluidblocks.
   -- We are assuming the fluidblocks have already been created
   o.id = #solidBlocks + #fluidBlocks
   solidBlocks[#solidBlocks+1] = o
   o.label = o.label or string.format("SolidBlock-%d", o.id)
   -- Must have a grid, initial temperature and properties
   assert(o.grid or o.gridMetadata, "need to supply a grid or its metadata")
   assert(o.initTemperature, "need to supply an initTemperature")
   assert(o.modelTag, "need to supply a solid thermal model")
   if _solidModels[o.modelTag] == nil then
      local msg = string.format("The solid model tag '%s' has not been registered in call to registerSolidModels().", o.modelTag)
      error(msg, 2)
   end
   
   -- Fill in some defaults, if not already set
   if o.active == nil then
      o.active = true
   end
   if o.bcList then
      o.bcList = deepclone(o.bcList, false)
   else
      o.bcList =  {} -- boundary conditions
   end
   -- Check the grid information.
   if o.grid then
      if config.dimensions ~= o.grid:get_dimensions() then
         local msg = string.format("Mismatch in dimensions, config %d grid %d.",
                                   config.dimensions, o.grid:get_dimensions())
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
         for _,face in ipairs(faceList(config.dimensions)) do
            o.bcList[face] = o.bcList[face] or SolidAdiabaticBC:new{}
         end
      end
      if o.grid:get_type() == "unstructured_grid" then
         error("Unstructured grids not available for solid blocks.")
      end
   else -- We have gridMetadata
      if config.dimensions ~= o.gridMetadata.dimensions then
         local msg = string.format("Mismatch in dimensions, config %d grid %d.",
                                   config.dimensions, o.gridMetadata.dimensions)
         error(msg)
      end
      if o.gridMetadata.type == "structured_grid" then
         -- Extract some information from the StructuredGrid
         -- Note 0-based indexing for vertices and cells.
         o.nic = o.gridMetadata.nic
         o.njc = o.gridMetadata.njc
         o.nkc = o.gridMetadata.nkc
         o.ncells = o.nic * o.njc * o.nkc
         -- The following table p for the corner locations,
         -- may be used later for testing for block connections.
         o.p = {}
         if config.dimensions == 3 then
            o.p[0] = o.gridMetadata.p0
            o.p[1] = o.gridMetadata.p1
            o.p[2] = o.gridMetadata.p2
            o.p[3] = o.gridMetadata.p3
            o.p[4] = o.gridMetadata.p4
            o.p[5] = o.gridMetadata.p5
            o.p[6] = o.gridMetadata.p6
            o.p[7] = o.gridMetadata.p7
         else
            o.p[0] = o.gridMetadata.p0
            o.p[1] = o.gridMetadata.p1
            o.p[2] = o.gridMetadata.p2
            o.p[3] = o.gridMetadata.p3
         end
         for _,face in ipairs(faceList(config.dimensions)) do
            o.bcList[face] = o.bcList[face] or SolidAdiabaticBC:new{}
         end
      end
      if o.gridMetadata.type == "unstructured_grid" then
         error("Unstructured grids not available for solid blocks.")
      end
   end
   return o
end

function SolidBlock:tojson()
   local str = string.format('"solid_block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "model": "%s",\n', self.modelTag)
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

return {
   SolidBlock = SolidBlock,
}
