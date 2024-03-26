-- Module for setting history points, required by prep.lua.
--
-- Authors: PJ, RJG and Kyle D.
--

local history = {}

local lmrconfig = require 'lmrconfig'
local globalconfig = require 'globalconfig'
config = globalconfig.config

function file_exists(name)
   -- Utility function copied from
   -- https://stackoverflow.com/questions/4990990/check-if-a-file-exists-with-lua
   local f = io.open(name, "r")
   return f ~= nil and io.close(f)
end

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
         if (blk.grid == nil) then
            -- Likely we've initialised only based on metadata
            if (blk.gridMetadata['type'] == 'structured_grid') then
               blk.grid = StructuredGrid:new{filename=lmrconfig.gridFilename(blk.id), fmt=config.grid_format}
            end
            if (blk.gridMetadata['type'] == 'unstructured_grid') then
               blk.grid = UnstructuredGrid:new{filename=lmrconfig.gridFilename(blk.id), fmt=config.grid_format}
            end
         end
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
      local blk = fluidBlocks[ib+1] -- one-based block array
      local nic = blk.nic
      local njc = blk.njc
      local nkc = blk.nkc
      -- Allow the user to specify the last index value in each direction
      -- by simply clipping the index values to actual range in the grid
      if i >= nic then
         i = nic-1
      end
      if j >= njc then
         j = njc-1
      end
      if k >= nkc then
         k = nkc-1
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
   local blk = fluidBlocks[ib+1] -- one-based block array
   if (blk.grid == nil) then
      -- Likely we've initialised only based on metadata
      if (blk.gridMetadata['type'] == 'structured_grid') then
         blk.grid = StructuredGrid:new{filename=lmrconfig.gridFilename(blk.id), fmt=config.grid_format}
      end
      if (blk.gridMetadata['type'] == 'unstructured_grid') then
         blk.grid = UnstructuredGrid:new{filename=lmrconfig.gridFilename(blk.id), fmt=config.grid_format}
      end
   end
   local pos = blk.grid:cellCentroid(i)
   print(string.format("Fluid History point [%d] ib=%d i=%d x=%g y=%g z=%g",
                       n, ib, i, pos.x, pos.y, pos.z))
   local fhpfile = "./lmrsim/fluid-history-points.list"
   if not file_exists(fhpfile) then
      local file = io.open(fhpfile, "w")
      file:write("blockid cellid pos.x pos.y pos.z\n")
      file:close()
   end
   if file_exists(fhpfile) then
      local file = io.open(fhpfile, "a")
      file:write(string.format("%d %d %g %g %g\n", ib, i, pos.x, pos.y, pos.z))
      file:close()
   end
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
