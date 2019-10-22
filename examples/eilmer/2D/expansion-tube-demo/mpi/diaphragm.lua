-- diaphragm.lua
--
-- These user-defined functions supplement the Exchange boundary condition
-- that would have been applied to the diaphragm boundary that sits between
-- upstream and downstream FluidBlocks, as identified in the preparation stage.
--
-- Initially it provides a reflecting-wall BC for each of the blocks, holding 
-- This is done by overwriting the previously-exchanged ghost-cell data.
-- When a high pressure signal arrives at the right end of the shock tube,
-- the diaphragm "ruptures" as indicated by the supervisor function.
-- These boundary-condition functions then stop overwriting the exchanged data,
-- such that it appears that the diaphragm has been instantly removed.
--
-- PJ, 2019-06-14, build for David Gildfind's simulations.
-- modified HH, 2019-10-14, PJ 2019-10-17

print("Initialize UDF for diaphragm for blkId=", blkId)
dofile('config/exptube.fluidBlockArrays')
local shock_tube = fluidBlockArrays[1]
blk_up_id = shock_tube.blockArray[shock_tube.nib][1][1]
local acc_tube = fluidBlockArrays[2]
blk_down_id = acc_tube.blockArray[1][1][1]
print("Block upstream of diaphragm=", blk_up_id, "downstream=", blk_down_id)
blk_up = infoFluidBlock(blk_up_id)
blk_down = infoFluidBlock(blk_down_id)

function ghostCells(args)
   local ruptureFlag = userPad[1]
   local cell0, cell1
   if blkId == blk_down_id and args.boundaryId == west then
      if ruptureFlag >= 1 then
	 -- Leave previous exchange data from upstream-block, east boundary.
	 cell0 = {}
	 cell1 = {}
      else
 	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(blk_down_id, blk_down.imin, args.j, blk_down.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(blk_down_id, blk_down.imin+1, args.j, blk_down.kmin)
	 cell1.velx = -cell1.velx
     end
   elseif blkId == blk_up_id and args.boundaryId == east then
      if ruptureFlag >= 1 then
	 -- Leave previous exchange data from downstream-block, west boundary.
	 cell0 = {}
	 cell1 = {}
      else
	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(blk_up_id, blk_up.imax, args.j, blk_up.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(blk_up_id, blk_up.imax-1, args.j, blk_up.kmin)
	 cell1.velx = -cell1.velx
      end
   else
      error(string.format("No ghost cell b.c. function available"..
			     " for block %d, on boundary: %d",
			     blkId, args.boundaryId))
   end 
   return cell0, cell1
end

function interface(args)
   return sampleFluidCell(blkId, args.i, args.j, args.k)
end
