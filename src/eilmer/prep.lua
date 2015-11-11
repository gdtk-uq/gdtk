-- prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- 
print("Loading prep.lua...")

require 'blk_conn'
-- Let's pull the symbols out of the blk_conn module
-- and make them global in this namespace
for k,v in pairs(blk_conn) do
   _G[k] = v
end

require 'bc'
for k,v in pairs(bc) do
   _G[k] = v
end

require 'gridpro'
-- Make these functions global so that users may call them
-- in the configuration script
applyGridproConnectivity = gridpro.applyGridproConnectivity

-- Storage for later definitions of Block objects
blocks = {}

-- Storgage for later definitions of SolidBlock objects
solidBlocks = {}


function to_eilmer_axis_map(gridpro_ijk)
   -- Convert from GridPro axis_map string to Eilmer3 axis_map string.
   -- From GridPro manual, Section 7.3.2 Connectivity Information.
   -- Example, 123 --> '+i+j+k'
   local axis_map = {[0]='xx', [1]='+i', [2]='+j', [3]='+k',
		     [4]='-i', [5]='-j', [6]='-k'}
   if type(gridpro_ijk) == "number" then
      gridpro_ijk = string.format("%03d", gridpro_ijk)
   end
   if type(gridpro_ijk) ~= "string" then
      print("Expected a string or integer of three digits but got:", gridpro_ijk)
      os.exit(-1)
   end
   eilmer_ijk = axis_map[tonumber(string.sub(gridpro_ijk, 1, 1))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 2, 2))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 3, 3))]
   return eilmer_ijk
end

-- -----------------------------------------------------------------------

-- Class for Block construction (based on a StructuredGrid).
SBlock = {
   myType = "SBlock" -- structured_grid block for gas dynamics
} -- end Block

function SBlock:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block, for later construction of the config file.
   -- Note that we want block id to start at zero for the D code.
   o.id = #(blocks)
   blocks[#(blocks)+1] = o
   -- Must have a grid and fillCondition
   assert(o.grid, "need to supply a grid")
   assert(o.fillCondition, "need to supply a fillCondition")
   -- Fill in default values, if already not set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("BLOCK-%d", o.id)
   o.omegaz = o.omegaz or 0.0
   o.bcList = o.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or WallBC_WithSlip:new()
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells.
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if config.dimensions == 3 then
      o.nkc = o.grid:get_nkv() - 1
   else
      o.nkc = 1
   end
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
   -- print("Block id=", o.id, "p0=", tostring(o.p[0]), "p1=", tostring(o.p[1]),
   --       "p2=", tostring(o.p[2]), "p3=", tostring(o.p[3]))
   return o
end

function SBlock:tojson()
   local str = string.format('"block_%d": {\n', self.id)
   str = str .. string.format('    "type": "%s",\n', self.myType)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   str = str .. string.format('    "omegaz": %f,\n', self.omegaz)
   str = str .. string.format('    "nhcell": %d,\n', #(self.hcellList))
   for i = 1, #(self.hcellList) do
      local hcell = self.hcellList[i]
      if config.dimensions == 3 then
	 str = str .. string.format('    "history-cell-%d": [%d, %d, %d],\n', 
				    i-1, hcell[1], hcell[2], hcell[3])
      else
	 str = str .. string.format('    "history-cell-%d": [%d, %d],\n',
				    i-1, hcell[1], hcell[2])
      end
   end
   -- Boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

-- -----------------------------------------------------------------------

-- Class for Block construction (based on an UnstructuredGrid).
UBlock = {
   myType = "UBlock" -- unstructured_grid block for gas dynamics
} -- end Block

function UBlock:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block, for later construction of the config file.
   -- Note that we want block id to start at zero for the D code.
   o.id = #(blocks)
   blocks[#(blocks)+1] = o
   -- Must have a grid and fillCondition
   assert(o.grid, "need to supply a grid")
   assert(o.fillCondition, "need to supply a fillCondition")
   -- Extract some information from the UnstructuredGrid
   o.ncells = o.grid:get_ncells()
   o.nvertices = o.grid:get_nvertices()
   o.nfaces = o.grid:get_nfaces()
   o.nboundaries = o.grid:get_nboundaries()
   -- Fill in default values, if already not set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("BLOCK-%d", o.id)
   o.omegaz = o.omegaz or 0.0
   o.bcList = o.bcList or {} -- boundary conditions
   -- [TODO] think about attaching boundary conditions via boundary labels or tags.
   for i = 0, o.nboundaries-1 do
      o.bcList[i] = o.bcList[i] or WallBC_WithSlip:new()
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
   return o
end

function UBlock:tojson()
   local str = string.format('"block_%d": {\n', self.id)
   str = str .. string.format('    "type": "%s",\n', self.myType)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "ncells": %d,\n', self.ncells)
   str = str .. string.format('    "nvertices": %d,\n', self.nvertices)
   str = str .. string.format('    "nfaces": %d,\n', self.nfaces)
   str = str .. string.format('    "nboundaries": %d,\n', self.nboundaries)
   str = str .. string.format('    "omegaz": %f,\n', self.omegaz)
   str = str .. string.format('    "nhcell": %d,\n', #(self.hcellList))
   for i = 1, #(self.hcellList) do
      local hcell = self.hcellList[i]
      str = str .. string.format('    "history-cell-%d": %d,\n', i-1, hcell)
   end
   -- Boundary conditions
   for i = 0, self.nboundaries-1 do
      str = str .. string.format('    "boundary_%d": ', i) ..
	 self.bcList[i]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

-- ---------------------------------------------------------------------------

function closeEnough(vA, vB, tolerance)
   -- Decide if two Vector quantities are close enough to being equal.
   -- This will be used to test that the block corners coincide.
   tolerance = tolerance or 1.0e-4
   return (vabs(vA - vB)/(vabs(vA + vB)+1.0)) <= tolerance
end

function connectBlocks(blkA, faceA, blkB, faceB, orientation)
   print("connectBlocks: blkA=", blkA.id, "faceA=", faceA, 
	 "blkB=", blkB.id, "faceB=", faceB, "orientation=", orientation)
   blkA.bcList[faceA] = ExchangeBC_FullFace:new{otherBlock=blkB.id, otherFace=faceB,
					       orientation=orientation}
   blkB.bcList[faceB] = ExchangeBC_FullFace:new{otherBlock=blkA.id, otherFace=faceA,
					       orientation=orientation}
   -- [TODO] need to test for matching corner locations and 
   -- consistent numbers of cells
end

function isPairInList(targetPair, pairList)
   local count = 0
   for _,v in ipairs(pairList) do
      if (v[1] == targetPair[1] and v[2] == targetPair[2]) or
	 (v[2] == targetPair[1] and v[1] == targetPair[2]) 
      then
	 count = count + 1
      end
   end
   return count > 0
end

function verticesAreCoincident(A, B, vtxPairs, tolerance)
   tolerance = tolerance or 1.0e-6
   local allVerticesAreClose = true
   for _,v in ipairs(vtxPairs) do
      -- print("A.id=", A.id, "B.id=", B.id, "vtxPair=", tostringVtxPair(v))
      -- print("  A.p=", tostring(A.p[v[1]]), "B.p=", tostring(B.p[v[2]]))
      if vabs(A.p[v[1]] - B.p[v[2]]) > tolerance then
	 allVerticesAreClose = false
      end
   end
   return allVerticesAreClose
end

function identifyBlockConnections(blockList, excludeList, tolerance)
   -- Identify block connections by trying to match corner points.
   -- Parameters:
   -- blockList: the list of SBlock objects to be included in the search.
   --    If nil, the whole collection is searched.
   -- excludeList: list of pairs of SBlock objects that should not be
   --    included in the search for connections.
   -- tolerance: spatial tolerance for the colocation of vertices
   blockList = blockList or blocks
   excludeList = excludeList or {}
   -- Put UBlock objects into the exclude list because they don't
   -- have a simple topology that can always be matched to an SBlock.
   for _,A in ipairs(blockList) do
      if A.myType == "UBlock" then excludeList[#excludeList+1] = A end
   end
   tolerance = tolerance or 1.0e-6
   for _,A in ipairs(blockList) do
      for _,B in ipairs(blockList) do
	 if (A ~= B) and (not isPairInList({A, B}, excludeList)) then
	    -- print("Proceed with test for coincident vertices.")
	    local connectionCount = 0
	    if config.dimensions == 2 then
	       -- print("2D test")
	       for vtxPairs,connection in pairs(connections2D) do
		  -- print("vtxPairs=", tostringVtxPairList(vtxPairs),
		  --       "connection=", tostringConnection(connection))
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
end

-- ---------------------------------------------------------------------------

function SBlockArray(t)
   -- Expect one table as argument, with named fields.
   -- Returns an array of blocks defined over a single region.
   assert(t.grid, "need to supply a grid")
   assert(t.fillCondition, "need to supply a fillCondition")
   t.omegaz = t.omegaz or 0.0
   t.bcList = t.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      t.bcList[face] = t.bcList[face] or WallBC_WithSlip:new()
   end
   t.xforceList = t.xforceList or {}
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
   for ib = 1, t.nib do
      blockArray[ib] = {}
      local i0 = (ib-1) * dnic
      if (ib == t.nib) then
	 -- Last block has to pick up remaining cells.
	 dnic = nic_total - i0
      end
      for jb = 1, t.njb do
	 local j0 = (jb-1) * dnjc
	 if (jb == t.njb) then
	    dnjc = njc_total - j0
	 end
	 if config.dimensions == 2 then
	    -- 2D flow
	    local subgrid = t.grid:subgrid(i0,dnic+1,j0,dnjc+1)
	    local bcList = {north=WallBC_WithSlip:new(), east=WallBC_WithSlip:new(),
			    south=WallBC_WithSlip:new(), west=WallBC_WithSlip:new()}
	    if ib == 1 then
	       bcList[west] = t.bcList[west]
	    end
	    if ib == t.nib then
	       bcList[east] = t.bcList[east]
	    end
	    if jb == 1 then
	       bcList[south] = t.bcList[south]
	    end
	    if jb == t.njb then
	       bcList[north] = t.bcList[north]
	    end
	    new_block = SBlock:new{grid=subgrid, omegaz=t.omegaz,
				   fillCondition=t.fillCondition, bcList=bcList}
	    blockArray[ib][jb] = new_block
	    blockCollection[#blockCollection+1] = new_block
	 else
	    -- 3D flow, need one more level in the array
	    blockArray[ib][jb] = {}
	    for kb = 1, t.nkb do
	       local k0 = (kb-1) * dnkc
	       if (kb == t.nkb) then
		  dnkc = nkc_total - k0
	       end
	       local subgrid = t.grid:subgrid(i0,dnic+1,j0,dnjc+1,k0,dnkc+1)
	       local bcList = {north=WallBC_WithSlip:new(), east=WallBC_WithSlip:new(),
			       south=WallBC_WithSlip:new(), west=WallBC_WithSlip:new(),
			       top=WallBC_WithSlip:new(), bottom=WallBC_WithSlip:new()}
	       if ib == 1 then
		  bcList[west] = t.bcList[west]
	       end
	       if ib == t.nib then
		  bcList[east] = t.bcList[east]
	       end
	       if jb == 1 then
		  bcList[south] = t.bcList[south]
	       end
	       if jb == t.njb then
		  bcList[north] = t.bcList[north]
	       end
	       if kb == 1 then
		  bcList[bottom] = t.bcList[bottom]
	       end
	       if kb == t.nkb then
		  bcList[top] = t.bcList[top]
	       end
	       new_block = SBlock:new{grid=subgrid, omegaz=t.omegaz,
				      fillCondition=t.fillCondition,
				      bcList=bcList}
	       blockArray[ib][jb][kb] = new_block
	       blockCollection[#blockCollection+1] = new_block
	    end -- kb loop
	 end -- dimensions
      end -- jb loop
   end -- ib loop
   -- Make the inter-subblock connections
   if #blockCollection > 1 then
      identifyBlockConnections(blockCollection)
   end
   return blockArray
end


-- Class for SolidBlock construction (based on a StructuredGrid)
SSolidBlock = {
   myType = "SSolidBlock"
} -- end SSolidBlock

function SSolidBlock:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new block for later construction of the config file.
   -- Note that we want the block id to start at zero for the D code.
   o.id = #solidBlocks
   solidBlocks[#solidBlocks+1] = o
   -- Must have a grid and initial temperature
   assert(o.grid, "need to supply a grid")
   assert(o.initTemperature, "need to supply an initTemperature")
   assert(o.properties, "need to supply physical properties for the block")
   -- Fill in some defaults, if not already set
   if o.active == nil then
      o.active = true
   end
   o.label = o.label or string.format("SOLIDBLOCK-%d", o.id)
   o.bcList = o.bcList or {} -- boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or SolidFixedTBC:new{Twall=300}
   end
   -- Extract some information from the StructuredGrid
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if config.dimensions == 3 then
      o.nkc = o.grid.get_nkv() - 1
   else
      o.nkc = 1
   end
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

function SSolidBlock:tojson()
   local str = string.format('"solid_block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   str = str .. '    "properties": {\n'
   str = str .. string.format('       "rho": %.6e,\n', self.properties.rho)
   str = str .. string.format('       "k": %.6e,\n', self.properties.k)
   str = str .. string.format('       "Cp": %.6e\n', self.properties.Cp)
   str = str .. '    },\n'
   -- Boundary conditions
      for _,face in ipairs(faceList(config.dimensions)) do
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
end

-- --------------------------------------------------------------------

function write_control_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"dt_init": %e,\n', config.dt_init))
   f:write(string.format('"dt_max": %e,\n', config.dt_max))
   f:write(string.format('"cfl_value": %e,\n', config.cfl_value))
   f:write(string.format('"stringent_cfl": %s,\n', tostring(config.stringent_cfl)))
   f:write(string.format('"fixed_time_step": %s,\n', tostring(config.fixed_time_step)))
   f:write(string.format('"dt_reduction_factor": %e,\n', config.dt_reduction_factor))
   f:write(string.format('"print_count": %d,\n', config.print_count))
   f:write(string.format('"cfl_count": %d,\n', config.cfl_count))
   f:write(string.format('"max_time": %e,\n', config.max_time))
   f:write(string.format('"max_step": %d,\n', config.max_step))
   f:write(string.format('"dt_plot": %e,\n', config.dt_plot))
   f:write(string.format('"dt_history": %e,\n', config.dt_history))
   f:write(string.format('"write_at_step": %d,\n', config.write_at_step))
   f:write(string.format('"halt_now": %d\n', config.halt_now))
   -- Note, also, no comma on last entry in JSON object. (^^^: Look up one line and check!)
   f:write("}\n")
   f:close()
end

function write_config_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"title": "%s",\n', config.title))
   f:write(string.format('"gas_model_file": "%s",\n', config.gas_model_file))
   f:write(string.format('"include_quality": %s,\n',
			 tostring(config.include_quality)))
   f:write(string.format('"dimensions": %d,\n', config.dimensions))
   f:write(string.format('"axisymmetric": %s,\n',
			 tostring(config.axisymmetric)))
   f:write(string.format('"interpolation_order": %d,\n', config.interpolation_order))
   f:write(string.format('"gasdynamic_update_scheme": "%s",\n',
			 config.gasdynamic_update_scheme))
   f:write(string.format('"MHD": %s,\n', tostring(config.MHD)))
   f:write(string.format('"separate_update_for_viscous_terms": %s,\n',
			 tostring(config.separate_update_for_viscous_terms)))
   f:write(string.format('"separate_update_for_k_omega_source": %s,\n', 
			 tostring(config.separate_update_for_k_omega_source)))
   f:write(string.format('"apply_bcs_in_parallel": %s,\n',
			 tostring(config.apply_bcs_in_parallel)))
   f:write(string.format('"max_invalid_cells": %d,\n', config.max_invalid_cells))
   f:write(string.format('"thermo_interpolator": "%s",\n', 
			 string.lower(config.thermo_interpolator)))
   f:write(string.format('"interpolate_in_local_frame": %s,\n', 
			 tostring(config.interpolate_in_local_frame)))
   f:write(string.format('"apply_limiter": %s,\n', tostring(config.apply_limiter)))
   f:write(string.format('"extrema_clipping": %s,\n', tostring(config.extrema_clipping)))

   f:write(string.format('"flux_calculator": "%s",\n', config.flux_calculator))
   f:write(string.format('"compression_tolerance": %e,\n', config.compression_tolerance))
   f:write(string.format('"shear_tolerance": %e,\n', config.shear_tolerance))
   f:write(string.format('"M_inf": %e,\n', config.M_inf))

   f:write(string.format('"grid_motion": "%s",\n', tostring(config.grid_motion)));
   f:write(string.format('"write_vertex_velocities": %s,\n', tostring(config.write_vertex_velocities)));
   f:write(string.format('"udf_grid_motion_file": "%s",\n', tostring(config.udf_grid_motion_file)));

   f:write(string.format('"viscous": %s,\n', tostring(config.viscous)))
   f:write(string.format('"spatial_deriv_calc": "%s",\n', config.spatial_deriv_calc))
   f:write(string.format('"viscous_signal_factor": %e,\n', config.viscous_signal_factor))

   f:write(string.format('"turbulence_model": "%s",\n',
			 string.lower(config.turbulence_model)))
   f:write(string.format('"turbulence_prandtl_number": %g,\n',
			 config.turbulence_prandtl_number))
   f:write(string.format('"turbulence_schmidt_number": %g,\n',
			 config.turbulence_schmidt_number))
   f:write(string.format('"max_mu_t_factor": %e,\n', config.max_mu_t_factor))
   f:write(string.format('"transient_mu_t_factor": %e,\n', config.transient_mu_t_factor))

   f:write(string.format('"udf_source_terms_file": "%s",\n', config.udf_source_terms_file))
   f:write(string.format('"udf_source_terms": %s,\n', tostring(config.udf_source_terms)))

   f:write(string.format('"reacting": %s,\n', tostring(config.reacting)))
   f:write(string.format('"reactions_file": "%s",\n', config.reactions_file))

   f:write(string.format('"control_count": %d,\n', config.control_count))
   f:write(string.format('"nblock": %d,\n', #(blocks)))

   f:write(string.format('"block_marching": %s,\n',
			 tostring(config.block_marching)))
   f:write(string.format('"nib": %d,\n', config.nib))
   f:write(string.format('"njb": %d,\n', config.njb))
   f:write(string.format('"nkb": %d,\n', config.nkb))
   f:write(string.format('"propagate_inflow_data": %s,\n',
			 tostring(config.propagate_inflow_data)))
 
   f:write(string.format('"udf_solid_source_terms_file": "%s",\n', config.udf_solid_source_terms_file))
   f:write(string.format('"udf_solid_source_terms": %s,\n', tostring(config.udf_solid_source_terms)))
   f:write(string.format('"nsolidblock": %d,\n', #solidBlocks))

   for i = 1, #blocks do
      f:write(blocks[i]:tojson())
   end
   for i = 1, #solidBlocks do
      f:write(solidBlocks[i]:tojson())
   end
  
   f:write('"dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
   f:write("}\n")

   f:close()
end

function write_times_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# tindx sim_time dt_global\n");
   f:write(string.format("%04d %12.6e %12.6e\n", 0, 0.0, config.dt_init))
   f:close()
end

function write_block_list_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# indx label\n")
   for i = 1, #(blocks) do
      f:write(string.format("%4d %s\n", blocks[i].id, blocks[i].label))
   end
   f:close()
end

function build_job_files(job)
   print("Build job files for ", job)
   write_config_file(job .. ".config")
   write_control_file(job .. ".control")
   write_times_file(job .. ".times")
   write_block_list_file(job .. ".list")
   os.execute("mkdir -p grid/t0000")
   os.execute("mkdir -p flow/t0000")
   if #solidBlocks >= 1 then
      os.execute("mkdir -p solid-grid/t0000")
      os.execute("mkdir -p solid/t0000")
   end
   for i = 1, #blocks do
      local id = blocks[i].id
      print("Block id=", id)
      local fileName = "grid/t0000/" .. job .. string.format(".grid.b%04d.t0000.gz", id)
      blocks[i].grid:write_to_gzip_file(fileName)
      local fileName = "flow/t0000/" .. job .. string.format(".flow.b%04d.t0000.gz", id)
      if blocks[i].myType == "SBlock" then
	 write_initial_sg_flow_file(fileName, blocks[i].grid, blocks[i].fillCondition, 0.0)
      else
	 write_initial_usg_flow_file(fileName, blocks[i].grid, blocks[i].fillCondition, 0.0)
      end
   end
   for i = 1, #solidBlocks do
      local id = solidBlocks[i].id
      print("SolidBlock id=", id)
      local fileName = "solid-grid/t0000/" .. job .. string.format(".solid-grid.b%04d.t0000.gz", id)
      solidBlocks[i].grid:write_to_gzip_file(fileName)
      local fileName = "solid/t0000/" .. job .. string.format(".solid.b%04d.t0000", id)
      writeInitialSolidFile(fileName, solidBlocks[i].grid,
			    solidBlocks[i].initTemperature, solidBlocks[i].properties, 0.0)
      os.execute("gzip -f " .. fileName)
   end

   print("Done building files.")
end

print("Done loading e4prep.lua")
