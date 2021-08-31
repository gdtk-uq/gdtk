-- Module for constructing FluidBlock objects, required by prep.lua.
--
-- Authors: PJ and RJG
--

module(..., package.seeall)

-- Class for gas dynamics FluidBlock construction.
FluidBlock = {
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
               error("FIX-ME set the unstructured-gird bcs")
               local tag = o.gridMetadata.get_boundaryset_tag(i) -- [FIX-ME]
               mybc = o.bcDict[tag]
            end
            mybc = mybc or WallBC_WithSlip:new() -- default boundary condition
            o.bcList[i] = mybc
         end
      end
   end
   return o
end

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
   str = str .. '},\n'
   return str
end

-- ---------------------------------------------------------------------------
function SBlock2UBlock(blk)
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
end

function connectBlocks(blkA, faceA, blkB, faceB, orientation)
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
   if blkA.myType == "FluidBlock" and blkB.myType == "FluidBlock" then
      blkA.bcList[faceA] = ExchangeBC_FullFace:new{otherBlock=blkB.id, otherFace=faceB,
						   orientation=orientation}
      blkB.bcList[faceB] = ExchangeBC_FullFace:new{otherBlock=blkA.id, otherFace=faceA,
						   orientation=orientation}
      -- [TODO] need to test for matching corner locations and consistent numbers of cells
   elseif blkA.myType == "FluidBlock" and blkB.myType == "SolidBlock" then
      blkA.bcList[faceA] = WallBC_AdjacentToSolid:new{otherBlock=blkB.id,
                                                      otherFace=faceB,
                                                      orientation=orientation}
      blkB.bcList[faceB] = SolidAdjacentToGasBC:new{otherBlock=blkA.id,
                                                    otherFace=faceA,
                                                    orientation=orientation}
   elseif blkA.myType == "SolidBlock" and blkB.myType == "FluidBlock" then
      blkA.bcList[faceA] = SolidAdjacentToGasBC:new{otherBlock=blkB.id,
                                                    otherFace=faceB,
                                                    orientation=orientation}
      blkB.bcList[faceB] = WallBC_AdjacentToSolid:new{otherBlock=blkA.id,
                                                      otherFace=faceA,
                                                      orientation=orientation}
   elseif blkA.myType == "SolidBlock" and blkB.myType == "SolidBlock" then
      blkA.bcList[faceA] = SolidFullFaceCopyBoundaryBC:new{otherBlock=blkB.id,
                                                           otherFace=faceB,
                                                           orientation=orientation}
      blkB.bcList[faceB] = SolidFullFaceCopyBoundaryBC:new{otherBlock=blkA.id,
                                                           otherFace=faceA,
                                                           orientation=orientation}
   end
end

function identifyBlockConnections(blockList, excludeList, tolerance)
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
      -- The caller has not provided a list; use the global blocks lists.
      for _,v in ipairs(fluidBlocks) do myBlockList[#myBlockList+1] = v end
      for _,v in ipairs(solidBlocks) do myBlockList[#myBlockList+1] = v end
   end
   excludeList = excludeList or {}
   -- Put UFluidBlock objects into the exclude list because they don't
   -- have a simple topology that can always be matched to an SFluidBlock.
   for _,A in ipairs(myBlockList) do
      if A.grid:get_type() == "unstructured_grid" then excludeList[#excludeList+1] = A end
   end
   tolerance = tolerance or 1.0e-6
   --
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
end -- identifyBlockConnections

-- Class for FluidBlock-Array objects.
FBArray = {
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
                                      "bcList", "nib", "njb", "nkb"})
   if not flag then
      error("Invalid name for item supplied to FBArray constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- We will embed the FBArray identity in the individual blocks
   -- and we would like that identity to start from 0 for the D code.
   o.id = #(fluidBlockArrays)
   --
   o.label = o.label or string.format("FluidBlockArray-%d", o.id)
   if fluidBlockArraysDict[o.label] then
      error('Have previously defined a FBArray with label "' .. o.label .. '"', 2)
   end
   fluidBlockArraysDict[o.label] = o.id
   if not o.initialState then
      -- try old name
      o.initialState = o.fillCondition
   end
   if not o.initialState then
      error("You need supply an initialState to FBArray constructor.", 2)
   end
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
      o.gridArray = {} -- will be a multi-dimensional array, indexed as [ib][jb][kb]
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
      -- Extract some information from the StructuredGrid
      o.niv = o.grid:get_niv()
      o.njv = o.grid:get_njv()
      o.nkv = o.grid:get_nkv()
      -- Subdivide the single grid based on numbers of cells.
      -- Note 0-based indexing for vertices and cells in the D-domain.
      local nic_total = o.niv - 1
      local dnic = math.floor(nic_total/o.nib)
      local njc_total = o.njv - 1
      local dnjc = math.floor(njc_total/o.njb)
      local nkc_total = o.nkv - 1
      local dnkc = math.floor(nkc_total/o.nkb)
      if config.dimensions == 2 then
         nkc_total = 1
         dnkc = 1
      end
      -- Work along each index direction and work out numbers of cells in subgrid.
      o.nics = {} -- numbers of cells in each subgrid
      local nic_remaining = nic_total
      for ib = 1, o.nib do
         local nic = math.floor(nic_remaining/(o.nib-ib+1))
         if (ib == o.nib) then
            -- On last subgrid, just use what's left
            nic = nic_remaining
         end
         o.nics[#o.nics+1] = nic
         nic_remaining = nic_remaining - nic
      end
      o.njcs = {}
      local njc_remaining = njc_total
      for jb = 1, o.njb do
         local njc = math.floor(njc_remaining/(o.njb-jb+1))
         if (jb == o.njb) then
            njc = njc_remaining
         end
         o.njcs[#o.njcs+1] = njc
         njc_remaining = njc_remaining - njc
      end
      o.nkcs = {}
      if config.dimensions == 2 then
         o.nkcs[1] = 1
      else
         local nkc_remaining = nkc_total
         for kb = 1, o.nkb do
            local nkc = math.floor(nkc_remaining/(o.nkb-kb+1))
            if (kb == o.nkb) then
               nkc = nkc_remaining
            end
            o.nkcs[#o.nkcs+1] = nkc
            nkc_remaining = nkc_remaining - nkc
         end
      end
      -- Now, generate the actual subgrids.
      local i0 = 0
      for ib = 1, o.nib do
         o.gridArray[ib] = {}
         local nic = o.nics[ib]
         local j0 = 0
         for jb = 1, o.njb do
            local njc = o.njcs[jb]
            if config.dimensions == 2 then
               -- 2D flow
               if false then
                  -- May activate print statements for debug.
                  print("ib=", ib, "jb= ", jb)
                  print("i0= ", i0, " nic= ", nic, " j0= ", j0, " njc= ", njc)
               end
               if nic < 1 then
                  error(string.format("Invalid nic=%d while making subgrid ib=%d, jb=%d", nic, ib, jb), 2)
               end
               if njc < 1 then
                  error(string.format("Invalid njc=%d while making subgrid ib=%d, jb=%d", njc, ib, jb), 2)
               end
               o.gridArray[ib][jb] = o.grid:subgrid(i0,nic+1,j0,njc+1)
            else
               -- 3D flow, need one more level in the array
               o.gridArray[ib][jb] = {}
               local k0 = 0
               for kb = 1, o.nkb do
                  local nkc = o.nkcs[kb]
                  if nic < 1 then
                     error(string.format("Invalid nic=%d while making subgrid ib=%d, jb=%d, kb=%d", nic, ib, jb, kb), 2)
                  end
                  if njc < 1 then
                     error(string.format("Invalid njc=%d while making subgrid ib=%d, jb=%d, kb=%d", njc, ib, jb, kb), 2)
                  end
                  if nkc < 1 then
                     error(string.format("Invalid nkc=%d while making subgrid ib=%d, jb=%d, kb=%d", nkc, ib, jb, kb), 2)
                  end
                  o.gridArray[ib][jb][kb] = o.grid:subgrid(i0,nic+1,j0,njc+1,k0,nkc+1)
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
      -- Finished generating subgrids
   else
      -- We were not given a single grid,
      -- so we assume that we were given the array of subgrids.
      --
      -- [TODO] 2021-08-31 PJ
      -- We need to check that the gridArray tables are regular and the individual grids align.
      -- Maybe we should do that while assembling a single overall grid.
      --
      o.nib = #(o.gridArray)
      o.njb = #(o.gridArray[1])
      if config.dimensions == 2 then
         -- Make stacks of the original subgrids.
         local stacksOfGrids = {}
         for ib = 1, o.nib do
            stacksOfGrids[ib] = StructuredGrid:clone(o.gridArray[ib][1]) -- [TODO] need this method
            for jb = 2, o.njb do
               stacksOfGrids[ib].joinGrid(o.gridArray[ib][jb], "north")
            end
         end
         o.grid = stacksOfGrids[1]
         for ib = 2, o.nib do
            o.grid.joinGrid(stacksOfGrids[ib], "east")
         end
      else
         -- For 3D
         o.nkb = #(o.gridArray[1][1])
         error("[TODO] In 3D, we need to assemble a single grid from the gridArray")
         -- Make stacks of the original subgrids, starting in the k-index direction,
         -- then building slabs spanning the jk directions from those stacks and,
         -- finally, joining the slabs in the i-direction.
      end
      -- Extract some information from the assembled StructuredGrid
      o.niv = o.grid:get_niv()
      o.njv = o.grid:get_njv()
      o.nkv = o.grid:get_nkv()
   end
   --
   -- At this point, we have an array of grids so generate the array of FluidBlocks.
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
	    local subgrid = o.gridArray[ib][jb]
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
	    local new_block = FluidBlock:new{grid=subgrid, omegaz=o.omegaz,
                                             initialState=o.initialState,
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
	       local subgrid = o.gridArray[ib][jb][kb]
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
	       local new_block = FluidBlock:new{grid=subgrid, omegaz=o.omegaz,
                                                initialState=o.initialState,
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
      o.niv = o.grid:get_niv()
      o.njv = o.grid:get_njv()
      o.nkv = o.grid:get_nkv()
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
            p0 = o.grid:get_vtx(i,j,k)
            distances[i] = 0.0
            for irev = 1, o.niv-1 do
               i = o.niv-1-irev
               p1 = o.grid:get_vtx(i,j,k)
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
function FluidBlockArray(t)
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
   str = str .. '},\n'
   return str
end
