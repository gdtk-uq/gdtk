-- Module to collect functions related to importing gridpro
-- boundary conditions and connectivity information.
--
-- Authors: RJG and PJ
-- Date: 2015-10-01
--

local gridpro = {}

-- Keep the following consistent with ws_ptymap.eilmer
-- The master copy of this is in examples/eilmer/3D/gridpro-import
local gproBCMap = {
  [4] = "WALL_SLIP",
  [5] = "WALL_ADIABATIC",
  [6] = "WALL_FIXED_T",
  [7] = "INFLOW_SUPERSONIC",
  [8] = "INFLOW_SUBSONIC",
  [9] = "INFLOW_SHOCKFITTING",
 [10] = "OUTFLOW_SIMPLE",
 [11] = "OUTFLOW_SUBSONIC",
 [12] = "OUTFLOW_FIXED_P",
 [13] = "USER_DEFINED1",
 [14] = "USER_DEFINED2",
 [15] = "USER_DEFINED3",
 [16] = "USER_DEFINED4",
 [17] = "USER_DEFINED5"
}

local function split_string(str)
   tokens = {}
   for tk in string.gmatch(str, "%S+") do
      tokens[#tokens+1] = tk
   end
   return tokens
end

local function to_eilmer_axis_map(gridpro_ijk)
   -- Convert from GridPro axis_map string to Eilmer3 axis_map string.
   -- From GridPro manual, Section 7.3.2 Connectivity Information.
   -- Example, 123 --> '+i+j+k'
   local axis_map = {[0]='xx', [1]='+i', [2]='+j', [3]='+k',
		     [4]='-i', [5]='-j', [6]='-k'}
   if type(gridpro_ijk) == "number" then
      gridpro_ijk = string.format("%03d", gridpro_ijk)
   end
   if type(gridpro_ijk) ~= "string" then
      error("Expected a string or integer of three digits but got:"..tostring(gridpro_ijk))
   end
   eilmer_ijk = axis_map[tonumber(string.sub(gridpro_ijk, 1, 1))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 2, 2))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 3, 3))]
   return eilmer_ijk
end

function gridpro.applyGridproConnectivity(fname, blks, gridpro2eilmer_id_map)
   print("Applying block connections from GridPro connectivity file: ", fname)
   -- the user may provide an id_map which defines the mapping of GridPro block ids to Eilmer block ids
   -- if an id_map is not provided, we will fill in an id_map with a 1:1 mapping
   if gridpro2eilmer_id_map == nil then
      gridpro2eilmer_id_map = {}
      for i,blk in ipairs(blks) do
         gridpro2eilmer_id_map[i] = i
      end
   end
   local f = assert(io.open(fname, 'r'))
   local line
   while true do
      line = f:read("*line")
      tks = split_string(line)
      if string.sub(tks[1], 1, 1) ~= '#' then
	 break
      end
   end
   local nb = tonumber(split_string(line)[1])
   local conns = {}
   for ib=1,nb do
      conns[gridpro2eilmer_id_map[ib]] = {}
      while true do
	 line = f:read("*line")
	 tks = split_string(line)
	 if string.sub(tks[1], 1, 1) ~= '#' then
	    break
	 end
      end
      -- Work on faces in order.
      -- Gridpro imin ==> Eilmer WEST face
      local otherBlk = gridpro2eilmer_id_map[tonumber(tks[5])]
      if otherBlk ~= nil and otherBlk > 0 then -- there is a connection
	 conns[gridpro2eilmer_id_map[ib]]["west"] = {otherBlk, tks[6]}
      end
      -- Gridpro imax ==> Eilmer EAST face
      otherBlk = gridpro2eilmer_id_map[tonumber(tks[9])]
      if otherBlk ~= nil and otherBlk > 0 then
	 conns[gridpro2eilmer_id_map[ib]]["east"] = {otherBlk, tks[10]}
      end
      -- Gridpro jmin ==> Eilmer SOUTH face
      otherBlk = gridpro2eilmer_id_map[tonumber(tks[13])]
      if otherBlk ~= nil and otherBlk > 0 then
	 conns[gridpro2eilmer_id_map[ib]]["south"] = {otherBlk, tks[14]}
      end
      -- Gridpro jmax ==> Eilmer NORTH face
      otherBlk = gridpro2eilmer_id_map[tonumber(tks[17])]
      if otherBlk ~= nil and otherBlk > 0 then
	 conns[gridpro2eilmer_id_map[ib]]["north"] = {otherBlk, tks[18]}
      end
      -- Gridpro kmin ==> Eilmer BOTTOM face
      otherBlk = gridpro2eilmer_id_map[tonumber(tks[21])]
      if otherBlk ~= nil and otherBlk > 0 then
	 conns[gridpro2eilmer_id_map[ib]]["bottom"] = {otherBlk, tks[22]}
      end
      -- Gridpro kmax ==> Eilmer TOP face
      otherBlk = gridpro2eilmer_id_map[tonumber(tks[25])]
      if otherBlk ~= nil and otherBlk > 0 then
	 conns[gridpro2eilmer_id_map[ib]]["top"] = {otherBlk, tks[26]}
      end
   end
   f:close()

   for ib=1,nb do
      for faceA, conn in pairs(conns[ib]) do
	 oblk = conn[1]
	 axisMap = conn[2]
	 A = blks[ib]
	 B = blks[oblk]
	 local faceB
	 for face, t in pairs(conns[oblk]) do
	    if t[1] == ib then
	       faceB = face
	       break
	    end
	 end
	 orientation = eilmer_orientation[faceA..faceB..to_eilmer_axis_map(axisMap)]
	 connectBlocks(A, faceA, B, faceB, orientation)
      end
   end

end

function gridpro.applyGridproBoundaryConditions(fname, blks, bcMap, dim, gridpro2eilmer_id_map)
   -- the user may provide an id_map which defines the mapping of GridPro block ids to Eilmer block ids
   -- if an id_map is not provided, we will fill in an id_map with a 1:1 mapping
   if gridpro2eilmer_id_map == nil then
      gridpro2eilmer_id_map = {}
      for i,blk in ipairs(blks) do
         gridpro2eilmer_id_map[i] = i
      end
   end
   local dim = dim or 3
   f = assert(io.open(fname, "r"))
   local line = f:read("*line")
   local tks = split_string(line)
   nBlocks = tonumber(tks[1])
   if nBlocks ~= #blks then
      print("Error in applyGridproBoundaryConditions(): mismatch in number of blocks.")
      print("The number of blocks given in the Gridpro property file (.pty) is ", nBlocks)
      print("But the number of blocks supplied in the block list is ", #blks)
      print("Bailing out.")
      os.exit(1)
   end
   bcs = {}
   for i=1,nBlocks do
      -- Loop past comment lines
      while true do
	 line = f:read("*line")
	 tks = split_string(line)
	 if string.sub(tks[1], 1, 1) ~= '#' then
	    -- We have a valid line.
	    break
	 end
      end
      bcs[gridpro2eilmer_id_map[i]] = {west=tonumber(tks[5]),
                                       east=tonumber(tks[7]),
                                       south=tonumber(tks[9]),
                                       north=tonumber(tks[11])}
      if dim == 3 then
	 bcs[gridpro2eilmer_id_map[i]].bottom = tonumber(tks[13])
	 bcs[gridpro2eilmer_id_map[i]].top = tonumber(tks[15])
      end
   end
   -- Read labels and discard
   line = f:read("*line")
   tks = split_string(line)
   nLabels = tonumber(tks[1])
   for i=1,nLabels do
      f:read("*line")
   end
   -- Read bcTypes and do something with them.
   line = f:read("*line")
   tks = split_string(line)
   nBCTypes = tonumber(tks[1])
   BCTypeMap = {}
   for ibc=1,nBCTypes do
      line = f:read("*line")
      tks = split_string(line)
      bcIdx = tonumber(tks[1])
      -- Gridpro seems to give the index as either an integer <= 32
      -- or a 4-digit integer. In this 4-digit integers, the last two
      -- digits encode the BC information
      bcInt = 0
      if #(tks[1]) == 4 then
          bcInt = tonumber(string.sub(tks[1], 3, 4))
      else
          bcInt = tonumber(tks[1])
      end
      BCTypeMap[tonumber(tks[1])] = bcInt
   end
   f:close()
   -- At this point all of the information has been gathered.
   -- Now loop over the blocks, and apply the BCs as appropriate.
   for ib, blk in ipairs(blks) do
      for face, bcID in pairs(bcs[ib]) do
	 bcInt = BCTypeMap[bcID]
	 bcLabel = gproBCMap[bcInt]
	 if bcLabel == nil then
	    bcLabel = string.format("user%d", bcInt)
	 end
	 if bcMap[bcLabel] then
	    blk.bcList[face] = bcMap[bcLabel]
         end
      end
   end
end

function gridpro.returnGridproGridBlockType(fname, fluid_domain_id)
   f = assert(io.open(fname, "r"))
   local line = f:read("*line")
   local tks = split_string(line)
   nBlocks = tonumber(tks[1])
   blk_type = {}
   for i=1,nBlocks do
      -- Loop past comment lines
      while true do
	 line = f:read("*line")
	 tks = split_string(line)
	 if string.sub(tks[1], 1, 1) ~= '#' then
	    -- We have a valid line.
	    break
	 end
      end
      if (tonumber(tks[3]) == fluid_domain_id) then
         blk_type[#blk_type+1] = "fluid"
      else
         blk_type[#blk_type+1] = "solid"
      end
   end
   return blk_type
end

function gridpro.to_eilmer_axis_map(gridpro_ijk)
   -- Convert from GridPro axis_map string to Eilmer3 axis_map string.
   -- From GridPro manual, Section 7.3.2 Connectivity Information.
   -- Example, 123 --> '+i+j+k'
   local axis_map = {[0]='xx', [1]='+i', [2]='+j', [3]='+k',
		     [4]='-i', [5]='-j', [6]='-k'}
   if type(gridpro_ijk) == "number" then
      gridpro_ijk = string.format("%03d", gridpro_ijk)
   end
   if type(gridpro_ijk) ~= "string" then
      error("Expected a string or integer of three digits but got:"..tostring(gridpro_ijk))
   end
   local eilmer_ijk = axis_map[tonumber(string.sub(gridpro_ijk, 1, 1))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 2, 2))] ..
      axis_map[tonumber(string.sub(gridpro_ijk, 3, 3))]
   return eilmer_ijk
end

return gridpro
