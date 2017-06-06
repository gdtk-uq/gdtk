-- Module to collect functions related to importing gridpro
-- boundary conditions and connectivity information.
--
-- Authors: RJG and PJ
-- Date: 2015-10-01
--
-- TODO: Port apply BCs function

module(..., package.seeall)

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

function applyGridproConnectivity(fname, blks)
   print("Applying block connections from GridPro connectivity file: ", fname)
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
      conns[#conns+1] = {}
      while true do
	 line = f:read("*line")
	 tks = split_string(line)
	 if string.sub(tks[1], 1, 1) ~= '#' then
	    break
	 end
      end
      -- Work on faces in order.
      -- Gridpro imin ==> Eilmer WEST face
      local otherBlk = tonumber(tks[5])
      if otherBlk > 0 then -- there is a connection
	 conns[ib]["west"] = {otherBlk, tks[6]}
      end
      -- Gridpro imax ==> Eilmer EAST face
      otherBlk = tonumber(tks[9])
      if otherBlk > 0 then 
	 conns[ib]["east"] = {otherBlk, tks[10]}
      end
      -- Gridpro jmin ==> Eilmer SOUTH face
      otherBlk = tonumber(tks[13])
      if otherBlk > 0 then
	 conns[ib]["south"] = {otherBlk, tks[14]}
      end
      -- Gridpro jmax ==> Eilmer NORTH face
      otherBlk = tonumber(tks[17])
      if otherBlk > 0 then
	 conns[ib]["north"] = {otherBlk, tks[18]}
      end
      -- Gridpro kmin ==> Eilmer BOTTOM face
      otherBlk = tonumber(tks[21])
      if otherBlk > 0 then
	 conns[ib]["bottom"] = {otherBlk, tks[22]}
      end
      -- Gridpro kmax ==> Eilmer TOP face
      otherBlk = tonumber(tks[25])
      if otherBlk > 0 then
	 conns[ib]["top"] = {otherBlk, tks[26]}
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

	 
      

      
