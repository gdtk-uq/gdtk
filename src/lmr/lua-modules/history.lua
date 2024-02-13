-- Module for setting history points, required by prep.lua.
--
-- Authors: PJ, RJG and Kyle D.
--

local history = {}

function history.setHistoryPoint(args)
   -- Accepts a variety of arguments:
   --  1. x, y, z coordinates
   --  setHistoryPoint{x=7.9, y=8.2, z=0.0}
   --  2. block and single-index for cell
   --  setHistoryPoint{ib=2, i=102}
   --  3. block and structured grid indices
   --  setHistoryPoint{ib=0, i=20, j=10, k=0}
   --
   local flag = checkAllowedNames(args, {"x", "y", "z", "ib", "i", "j", "k"})
   if not flag then
      error("Invalid name for item supplied to setHistoryPoint.", 2)
   end
   -- First look for x,y,z
   local found = false
   if args.x then
      local x = args.x
      local y = args.y
      local z = args.z or 0.0
      local minDist = 1.0e9 -- something very large
      local blkId = 0
      local cellId = 0
      for ib,blk in ipairs(fluidBlocks) do
	 local indx, dist = blk.grid:find_nearest_cell_centre{x=x, y=y, z=z}
	 if (dist < minDist) then
	    minDist = dist
	    blkId = ib
	    cellId = indx
	 end
      end
      -- Convert blkId to 0-offset
      blkId = blkId - 1
      historyCells[#historyCells+1] = {ib=blkId, i=cellId}
      found = true
   end
   -- Still trying; look for integer indices.
   if (not found) and args.j and args.ib then
      local ib = args.ib
      local i = args.i
      local j = args.j
      local k = args.k or 0
      local nic = fluidBlocks[ib+1].nic
      local njc = fluidBlocks[ib+1].njc
      local nkc = fluidBlocks[ib+1].nkc
      -- Check indices are valid.
      if i >= nic then
         errMsg = string.format("setHistoryPoint: i value invalid; valid i --> 0 <= i < nic.\n i= %d, nic= %d", i, nic)
         error(errMsg, 2)
      end
      if j >= njc then
         errMsg = string.format("setHistoryPoint: j value invalid. valid j --> 0 <= j < njc.\n j= %d, njc= %d", j, njc)
         error(errMsg, 2)
      end
      if k >= nkc then
         errMsg = string.format("setHistoryPoint: k value invalid. valid k --> 0 <= k < nkc. \n k= %d, nkc= %d", k, nkc)
         error(errMsg, 2)
      end
      -- Convert back to single_index
      local cellId = k * (njc * nic) + j * nic + i
      historyCells[#historyCells+1] = {ib=args.ib, i=cellId}
      found = true
   end
   if not found then
      -- Final option is to directly set the identity block and cell.
      if args.ib and args.i then
         historyCells[#historyCells+1] = {ib=args.ib, i=args.i}
      else
         error("Could not identify cell for setHistoryPoint.", 2)
      end
   end
   -- If we arrive here, we have successfully set the cell identity.
   -- Print a summary of its identity and location.
   local n = #historyCells
   local ib = historyCells[n].ib -- zero-based block index
   local i = historyCells[n].i
   local pos = fluidBlocks[ib+1].grid:cellCentroid(i) -- one-based block array
   print(string.format("Fluid History point [%d] ib=%d i=%d x=%g y=%g z=%g",
                       n, ib, i, pos.x, pos.y, pos.z))
   return
end

function history.setSolidHistoryPoint(args)
   -- Accepts a variety of arguments:
   --  1. x, y, z coordinates
   --  setSolidHistoryPoint{x=7.9, y=8.2, z=0.0}
   --  2. block and single-index for cell
   --  setSolidHistoryPoint{ib=2, i=102}
   --  3. block and structured grid indices
   --  setSolidHistoryPoint{ib=0, i=20, j=10, k=0}
   --
   local flag = checkAllowedNames(args, {"x", "y", "z", "ib", "i", "j", "k"})
   if not flag then
      error("Invalid name for item supplied to setSolidHistoryPoint.", 2)
   end
   -- First look for x,y,z
   local found = false
   if args.x then
      local x = args.x
      local y = args.y
      local z = args.z or 0.0
      local minDist = 1.0e9 -- something very large
      local blkId = 0
      local cellId = 0
      for ib,blk in ipairs(solidBlocks) do
	 local indx, dist = blk.grid:find_nearest_cell_centre{x=x, y=y, z=z}
	 if (dist < minDist) then
	    minDist = dist
	    blkId = ib
	    cellId = indx
	 end
      end
      -- Convert blkId to 0-offset
      blkId = blkId - 1
      solidHistoryCells[#solidHistoryCells+1] = {ib=blkId, i=cellId}
      found = true
   end
   -- Still trying; look for integer indices.
   if (not found) and args.j and args.ib then
      local ib = args.ib
      local i = args.i
      local j = args.j
      local k = args.k or 0
      -- Convert back to single_index
      local nic = solidBlocks[ib+1].nic
      local njc = solidBlocks[ib+1].njc
      local cellId = k * (njc * nic) + j * nic + i
      solidHistoryCells[#solidHistoryCells+1] = {ib=args.ib, i=cellId}
      found = true
   end
   if not found then
      -- Final option is to directly set the identity block and cell.
      if args.ib and args.i then
         solidHistoryCells[#solidHistoryCells+1] = {ib=args.ib, i=args.i}
      else
         error("Could not identify cell for setSolidHistoryPoint.", 2)
      end
   end
   -- If we arrive here, we have successfully set the cell identity.
   -- Print a summary of its identity and location.
   local n = #solidHistoryCells
   local ib = solidHistoryCells[n].ib -- zero-based block index
   local i = solidHistoryCells[n].i
   local pos = solidBlocks[ib+1].grid:cellCentroid(i) -- one-based block array
   print(string.format("Solid History point [%d] ib=%d i=%d x=%g y=%g z=%g",
                       n, ib, i, pos.x, pos.y, pos.z))
   return
end

return history
