-- Module to collect functions related to importing gridpro
-- boundary conditions and connectivity information.
--
-- Authors: RJG and PJ
-- Date: 2015-10-01
-- Update: 2023-11-28
--         This update is related to the staged-prep style
--         favoured in lmr5. We need to work at the grid level,
--         not the block level.
--

local gridproimport = {}

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

--[
-- We assume that the user has already register grids because we will look in the gridsList.
--]
function gridproimport.importGridproConnectivity(fname)
   print("Import and apply block connections from GridPro connectivity file: ", fname)
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
      conns[ib] = {}
      while true do
	 line = f:read("*line")
	 tks = split_string(line)
	 if string.sub(tks[1], 1, 1) ~= '#' then
	    break
	 end
      end
      -- Work on faces in order.
      -- Gridpro imin ==> Eilmer WEST face
      local otherGrid = tonumber(tks[5])
      if otherGrid ~= nil and otherGrid > 0 then -- there is a connection
	 conns[ib]["west"] = {otherGrid, tks[6]}
      end
      -- Gridpro imax ==> Eilmer EAST face
      otherGrid = tonumber(tks[9])
      if otherGrid ~= nil and otherGrid > 0 then
	 conns[ib]["east"] = {otherGrid, tks[10]}
      end
      -- Gridpro jmin ==> Eilmer SOUTH face
      otherGrid = tonumber(tks[13])
      if otherGrid ~= nil and otherGrid > 0 then
	 conns[ib]["south"] = {otherGrid, tks[14]}
      end
      -- Gridpro jmax ==> Eilmer NORTH face
      otherGrid = tonumber(tks[17])
      if otherGrid ~= nil and otherGrid > 0 then
	 conns[ib]["north"] = {otherGrid, tks[18]}
      end
      -- Gridpro kmin ==> Eilmer BOTTOM face
      otherGrid = tonumber(tks[21])
      if otherGrid ~= nil and otherGrid > 0 then
	 conns[ib]["bottom"] = {otherGrid, tks[22]}
      end
      -- Gridpro kmax ==> Eilmer TOP face
      otherGrid = tonumber(tks[25])
      if otherGrid ~= nil and otherGrid > 0 then
	 conns[ib]["top"] = {otherGrid, tks[26]}
      end
   end
   f:close()

   for ib=1,nb do
      for faceA, conn in pairs(conns[ib]) do
	 oGrid = conn[1]
	 axisMap = conn[2]
	 idA = ib-1 -- count from zero
	 idB = oGrid-1
	 local faceB
	 for face, t in pairs(conns[oGrid]) do
	    if t[1] == ib then
	       faceB = face
	       break
	    end
	 end
	 orientation = eilmer_orientation[faceA..faceB..to_eilmer_axis_map(axisMap)]
	 connectGrids(idA, faceA, idB, faceB, orientation)
      end
   end

end

--[
-- We assume that the user has already register grids because we will look in the gridsList.
-- As per importing connections, we assume the user wants all BCs applied.
--]
function gridproimport.importGridproBCs(fname)
   f = assert(io.open(fname, "r"))
   local line = f:read("*line")
   local tks = split_string(line)
   nGrids = tonumber(tks[1])
   if nGrids ~= #gridsList then
      print("Error in importGridproBoundaryConditions(): mismatch in number of grids.")
      print("The number of grids given in the Gridpro property file (.pty) is ", nGrids)
      print("But the number of grids available in gridsList is ", #gridsList)
      print("Bailing out.")
      os.exit(1)
   end
   bcs = {}
   for i=1,nGrids do
      -- Loop past comment lines
      while true do
	 line = f:read("*line")
	 tks = split_string(line)
	 if string.sub(tks[1], 1, 1) ~= '#' then
	    -- We have a valid line.
	    break
	 end
      end
      bcs[i] = {west=tonumber(tks[5]),
	        east=tonumber(tks[7]),
                south=tonumber(tks[9]),
                north=tonumber(tks[11]),
                bottom = tonumber(tks[13]),
	        top = tonumber(tks[15])}
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
   -- Now loop over the grids, and apply the BCs as appropriate.
   for ig, grid in ipairs(gridsList) do
      for face, bcID in pairs(bcs[ig]) do
	 bcInt = BCTypeMap[bcID]
	 bcLabel = gproBCMap[bcInt]
	 if bcLabel == nil then
	    bcLabel = string.format("user%d", bcInt)
	 end
	 grid.bcTags[face] = bcLabel
      end
   end
end

return gridproimport
