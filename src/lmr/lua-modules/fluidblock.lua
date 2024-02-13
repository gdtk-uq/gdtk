-- Module for constructing FluidBlock objects, required by prep.lua.
--
-- Authors: PJ and RJG
--

require 'blk_conn'

-- Class for gas dynamics FluidBlock construction.
local FluidBlock = {
   myType = "FluidBlock",
} -- end FluidBlock

function FluidBlock:new(o)
   local flag = type(self)=='table' and self.myType=='FluidBlock'
   if not flag then
      error("Make sure that you are using FluidBlock:new{} and not FluidBlock.new{}", 2)
   end
   o = o or {}
   flag = checkAllowedNames(o, {"grid", "gridMetadata", "initialState", "fillCondition", "active",
                                "label", "omegaz", "may_be_turbulent", "bcList", "bcDict",
                                "hcellList", "xforceList", "fluidBlockArrayId"})
   if not flag then
      error("Invalid name for item supplied to FluidBlock constructor.", 2)
   end
   if o.initialState == nil then o.initialState = o.fillCondition end -- try old name
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block, for later construction of the config file.
   -- Note that we want block id to start at zero for the D code.
   o.id = #(fluidBlocks)
   fluidBlocks[#(fluidBlocks)+1] = o
   o.label = o.label or string.format("FluidBlock-%d", o.id)
   if fluidBlocksDict[o.label] then
      error('Have previously defined a FluidBlock with label "' .. o.label .. '"', 2)
   end
   fluidBlocksDict[o.label] = o.id
   -- Set to -1 if NOT part of a fluid-block-array, otherwise use supplied value
   o.fluidBlockArrayId = o.fluidBlockArrayId or -1
   -- Must have a grid and initialState
   assert(o.grid or o.gridMetadata, "need to supply a grid or its metadata")
   assert(o.initialState, "need to supply an initialState")
   if getmetatable(o.initialState) == FlowSolution then
      -- Let's build a initialState function here from a FlowSolution.
      o.initialState = makeFlowStateFn(o.initialState)
   end
   -- Fill in default values, if already not set
   if o.active == nil then
      o.active = true
   end
   o.omegaz = o.omegaz or 0.0
   if o.may_be_turbulent == nil then
      o.may_be_turbulent = true
   end
   if o.bcList then
      o.bcList = deepclone(o.bcList, false)
   else
      o.bcList = {}
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
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
         -- print("FluidBlock id=", o.id, "p0=", tostring(o.p[0]), "p1=", tostring(o.p[1]),
         --       "p2=", tostring(o.p[2]), "p3=", tostring(o.p[3]))
         -- Attach default boundary conditions for those not specified.
         -- Note that, in the classic preparation, the structured-grid bcs arrive in bcList.
         for _,face in ipairs(faceList(config.dimensions)) do
            o.bcList[face] = o.bcList[face] or WallBC_WithSlip:new()
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
            local mybc = o.bcList[i]
            if (mybc == nil) and o.bcDict then
               local tag = o.grid:get_boundaryset_tag(i)
               mybc = o.bcDict[tag]
            end
            mybc = mybc or WallBC_WithSlip:new() -- default boundary condition
            o.bcList[i] = mybc
         end
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
         -- print("FluidBlock id=", o.id, "p0=", tostring(o.p[0]), "p1=", tostring(o.p[1]),
         --       "p2=", tostring(o.p[2]), "p3=", tostring(o.p[3]))
         -- Attach default boundary conditions for those not specified.
         -- Note that, in the staged-preparation, the bcs arrive as bcDict.
         for _,face in ipairs(faceList(config.dimensions)) do
            o.bcList[face] = o.bcDict[face] or WallBC_WithSlip:new()
         end
      end
      if o.gridMetadata.type == "unstructured_grid" then
         -- Extract some information from the UnstructuredGrid
         o.ncells = o.gridMetadata.ncells
         o.nvertices = o.gridMetadata.nvertices
         o.nfaces = o.gridMetadata.nfaces
         o.nboundaries = o.gridMetadata.nboundaries
         -- Attach boundary conditions from list or from the dictionary of conditions.
         for i = 0, o.nboundaries-1 do
            local mybc = o.bcList[i]
            if (mybc == nil) and o.bcDict then
               local key = tostring(i)
               local tag = o.gridMetadata.bcTags[key]
               mybc = o.bcDict[tag]
            end
            mybc = mybc or WallBC_WithSlip:new() -- default boundary condition
            o.bcList[i] = mybc
         end
      end
   end
   return o
end -- FluidBlock:new(o)

function FluidBlock:tojson()
   local str = string.format('"block_%d": {\n', self.id)
   str = str .. string.format('    "type": "%s",\n', self.myType)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "fluidBlockArrayId": %d,\n', self.fluidBlockArrayId)
   str = str .. string.format('    "omegaz": %.18e,\n', self.omegaz)
   str = str .. string.format('    "may_be_turbulent": %s,\n', tostring(self.may_be_turbulent))
   local grid_type
   if self.grid then
      grid_type = self.grid:get_type()
   else
      grid_type = self.gridMetadata.type
   end
   str = str .. string.format('    "grid_type": "%s",\n', grid_type)
   if grid_type == "structured_grid" then
      str = str .. string.format('    "nic": %d,\n', self.nic)
      str = str .. string.format('    "njc": %d,\n', self.njc)
      str = str .. string.format('    "nkc": %d,\n', self.nkc)
      -- Boundary conditions for structured grid.
      for _,face in ipairs(faceList(config.dimensions)) do
	 if not self.bcList[face].is_gas_domain_bc then
	    local msg = string.format("Boundary condition problem for block:%d, face:%s\n", self.id, face)
	    msg = msg.."   This boundary condition should be a gas domain b.c.\n"
	    msg = msg.."   The preparation stage cannot complete successfully.\n"
	    error(msg)
	 end
	 if not self.bcList[face].is_configured then
	    local msg = string.format("Boundary condition problem for block:%d, face:%s\n", self.id, face)
	    msg = msg.."   This boundary condition was not configured correctly.\n"
	    msg = msg.."   If you used one of the standard boundary conditions,\n"
	    msg = msg.."   did you remember to call the b.c constructor as bcName:new{}?\n"
	    msg = msg.."   If you have custom configured the boundary condition,\n"
	    msg = msg.."   did you remember to set the 'is_configured' flag to true?\n"
	    error(msg)
	 end
	 str = str .. string.format('    "boundary_%s": ', face) ..
	    self.bcList[face]:tojson() .. ',\n'
      end
   end
   if grid_type == "unstructured_grid" then
      str = str .. string.format('    "ncells": %d,\n', self.ncells)
      str = str .. string.format('    "nvertices": %d,\n', self.nvertices)
      str = str .. string.format('    "nfaces": %d,\n', self.nfaces)
      str = str .. string.format('    "nboundaries": %d,\n', self.nboundaries)
      -- Boundary conditions for the unstructured grid
      for i = 0, self.nboundaries-1 do
	 if not self.bcList[i].is_gas_domain_bc then
	    local msg = string.format("Boundary condition problem for block:%d, boundary:%d\n", self.id, i)
	    msg = msg.."   This boundary condition should be a gas domain b.c.\n"
	    msg = msg.."   The preparation stage cannot complete successfully.\n"
	    error(msg)
	 end
	 if not self.bcList[i].is_configured then
	    local msg = string.format("Boundary condition problem for block:%d, boundary:%d\n", self.id, i)
	    msg = msg.."   This boundary condition was not configured correctly.\n"
	    msg = msg.."   If you used one of the standard boundary conditions,\n"
	    msg = msg.."   did you remember to call the b.c constructor as bcName:new{}?\n"
	    msg = msg.."   If you have custom configured the boundary condition,\n"
	    msg = msg.."   did you remember to set the 'is_configured' flag to true?\n"
	    error(msg)
	 end
	 str = str .. string.format('    "boundary_%d": ', i) ..
	    self.bcList[i]:tojson() .. ',\n'
      end
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '}'
   return str
end -- FluidBlock:tojson()

-- ---------------------------------------------------------------------------
local function SBlock2UBlock(blk)
   local origId = blk.id
   local origLabel = blk.label
   -- Let's swap out any exchange_over_full_face BCs and replace
   -- with exchange BCs.
   local bcList = {}
   for i,face in ipairs(faceList(config.dimensions)) do
      if blk.bcList[face].type == "exchange_over_full_face" then
	 -- We'll convert any exchange_over_full_face BCs
	 bcList[i-1] = ExchangeBC_MappedCell:new{}
      else
	 -- For all other BCs, we directly copy.
	 bcList[i-1] = blk.bcList[face]
      end
   end
   ublk = FluidBlock:new{grid=UnstructuredGrid:new{sgrid=blk.grid},
			 initialState=blk.initialState,
			 active=blk.active,
			 label=nil,
			 omegaz=blk.omegaz,
			 bcList=bcList}
   local newId = ublk.id
   local newLabel = ublk.label
   -- Swap blocks in global list
   fluidBlocks[origId+1], fluidBlocks[newId+1] = fluidBlocks[newId+1], fluidBlocks[origId+1]
   -- Fix id and label of ublk
   fluidBlocks[origId+1].id = origId
   fluidBlocks[origId+1].label = origLabel
   -- Now remove original SFluidBlock, which is presently in pos ublk.id+1 and
   -- remove the latest dictionary entry, since the ublk no longer has that label.
   table.remove(fluidBlocks, newId+1)
   fluidBlocksDict[newLabel] = nil
end -- SBlock2UBlock()

local function connectBlocks(blkA, faceA, blkB, faceB, orientation)
   -- Make a "full-face" connection between a pair of Block objects.
   -- The connection is made by attaching boundary conditions
   -- to each block that reference the other block.
   if false then
      -- To reduce visual clutter at prep time, don't use the following print statement.
      -- It may still be useful for debugging.
      print("connectBlocks: blkA.id=", blkA.id, "faceA=", faceA,
            "blkB.id=", blkB.id, "faceB=", faceB, "orientation=", orientation)
   end
   -- Note that, when doing staged preparation, we may not have the grid object
   -- actually present at the point of making a FluidBlock connection.
   if ((blkA.grid and (blkA.grid:get_type() ~= "structured_grid")) or
      (blkA.grid and (blkB.grid:get_type() ~= "structured_grid"))) then
      error("connectBlocks() Works only for structured-grid blocks.", 2)
   end
   assert(type(faceA) == 'string' and type(faceB) == 'string',
          "connectBlocks requires faceA and faceB to be strings.")
   --
   if blkA.myType == "FluidBlock" and blkB.myType == "FluidBlock" then
      blkA.bcList[faceA] = ExchangeBC_FullFace:new{otherBlock=blkB.id, otherFace=faceB,
						   orientation=orientation}
      blkB.bcList[faceB] = ExchangeBC_FullFace:new{otherBlock=blkA.id, otherFace=faceA,
						   orientation=orientation}
      -- [TODO] need to test for matching corner locations and consistent numbers of cells
   elseif blkA.myType == "FluidBlock" and blkB.myType == "SolidBlock" then
      -- The preprocessor will set the time-accurate fluid/solid coupling BC by default unless the user explicitly requests the alternate
      flag = config.coupling_with_solid_domains == "steady_fluid_transient_solid"
      if not flag then
         blkA.bcList[faceA] = WallBC_AdjacentToSolid:new{otherBlock=blkB.id,
                                                         otherFace=faceB,
                                                         orientation=orientation}
         blkB.bcList[faceB] = SolidAdjacentToGasBC:new{otherBlock=blkA.id,
                                                       otherFace=faceA,
                                                       orientation=orientation}
      else
         blkA.bcList[faceA] = WallBC_AdjacentToSolid2:new{otherBlock=blkB.id,
                                                          otherFace=faceB,
                                                          orientation=orientation,
                                                          catalytic_type=blkA.bcList[faceA].catalytic_type,
                                                          wall_massf_composition=blkA.bcList[faceA].wall_massf_composition}
         blkB.bcList[faceB] = SolidAdjacentToGasBC2:new{otherBlock=blkA.id,
                                                        otherFace=faceA,
                                                        orientation=orientation}
      end
   elseif blkA.myType == "SolidBlock" and blkB.myType == "FluidBlock" then
      -- The preprocessor will set the time-accurate fluid/solid coupling BC by default unless the user explicitly requests the alternate
      flag = config.coupling_with_solid_domains == "steady_fluid_transient_solid"
      if not flag then
         blkA.bcList[faceA] = SolidAdjacentToGasBC:new{otherBlock=blkB.id,
                                                       otherFace=faceB,
                                                       orientation=orientation}
         blkB.bcList[faceB] = WallBC_AdjacentToSolid:new{otherBlock=blkA.id,
                                                         otherFace=faceA,
                                                         orientation=orientation}
      else
         blkA.bcList[faceA] = SolidAdjacentToGasBC2:new{otherBlock=blkB.id,
                                                        otherFace=faceB,
                                                        orientation=orientation}
         blkB.bcList[faceB] = WallBC_AdjacentToSolid2:new{otherBlock=blkA.id,
                                                          otherFace=faceA,
                                                          orientation=orientation,
                                                          catalytic_type=blkB.bcList[faceB].catalytic_type,
                                                          wall_massf_composition=blkB.bcList[faceB].wall_massf_composition}
      end
   elseif blkA.myType == "SolidBlock" and blkB.myType == "SolidBlock" then
      blkA.bcList[faceA] = SolidFullFaceCopyBoundaryBC:new{otherBlock=blkB.id,
                                                           otherFace=faceB,
                                                           orientation=orientation}
      blkB.bcList[faceB] = SolidFullFaceCopyBoundaryBC:new{otherBlock=blkA.id,
                                                           otherFace=faceA,
                                                           orientation=orientation}
   end
end -- connectBlocks()

local function identifyBlockConnections(blockList, excludeList, tolerance, imported_cht_grid)
   -- Identify block connections by trying to match corner points.
   -- Parameters (all optional):
   -- blockList: the list of SFluidBlock objects to be included in the search.
   --    If nil, the whole collection is searched.
   -- excludeList: list of pairs of SFluidBlock objects that should not be
   --    included in the search for connections.
   -- tolerance: spatial tolerance for the colocation of vertices
   --
   local myBlockList = {}

   if blockList then
      -- The caller has provided a list of Blocks to bound the search.
      for _,v in ipairs(blockList) do myBlockList[#myBlockList+1] = v end
   else
      if imported_cht_grid then
         -- The caller has flagged that this is an imported fluid/solid grid,
         -- so we only need to search for connections between fluid/solid blocks,
         -- and potentially any solid/solid blocks that don't have a connectivity map
         -- from GridPro (e.g. if the solid domain consists of several different materials).
         for _,v in ipairs(fluidBlocks) do
            for _,bc in pairs(v.bcList) do
               if bc.type == "wall_adjacent_to_solid" or bc.type == "wall_adjacent_to_solid2" then
                  myBlockList[#myBlockList+1] = v
               end
            end
         end
         for _,v in ipairs(solidBlocks) do
            for _,bc in pairs(v.bcList) do
               if bc.type == "SolidAdjacentToGas" or bc.type == "SolidAdjacentToGas2" or bc.type == "SolidAdiabatic" then
                  myBlockList[#myBlockList+1] = v
               end
            end
         end
      else
         -- The caller has not provided a list; use the global blocks lists.
         for _,v in ipairs(fluidBlocks) do myBlockList[#myBlockList+1] = v end
         for _,v in ipairs(solidBlocks) do myBlockList[#myBlockList+1] = v end
      end
   end
   excludeList = excludeList or {}
   -- Put UFluidBlock objects into the exclude list because they don't
   -- have a simple topology that can always be matched to an SFluidBlock.
   for _,A in ipairs(myBlockList) do
      if A.grid:get_type() == "unstructured_grid" then excludeList[#excludeList+1] = A end
   end
   tolerance = tolerance or 1.0e-6
   for _,A in ipairs(myBlockList) do
      for _,B in ipairs(myBlockList) do
	 if (A ~= B) and (not isPairInList({A, B}, excludeList)) then
	    -- print("Proceed with test for coincident vertices.") -- DEBUG
	    local connectionCount = 0
	    if config.dimensions == 2 then
	       -- print("2D test A.id=", A.id, " B.id=", B.id) -- DEBUG
	       for vtxPairs,connection in pairs(connections2D) do
		  -- print("vtxPairs=", tostringVtxPairList(vtxPairs),
		  --       "connection=", tostringConnection(connection)) -- DEBUG
                  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
		     local faceA, faceB = unpack(connection)
		     connectBlocks(A, faceA, B, faceB, 0)
		     connectionCount = connectionCount + 1
		  end
	       end
	    else
	       -- print("   3D test")
               for vtxPairs,connection in pairs(connections3D) do
		  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
		     local faceA, faceB, orientation = unpack(connection)
		     connectBlocks(A, faceA, B, faceB, orientation)
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
end -- identifyBlockConnections()

return {
   FluidBlock = FluidBlock,
   SBlock2UBlock = SBlock2UBlock,
   connectBlocks = connectBlocks,
   identifyBlockConnections = identifyBlockConnections
}
