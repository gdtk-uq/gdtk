-- diaphragm.lua
--
-- These user-defined functions supplement the Exchange boundary condition
-- that would have been applied to the diaphragm boundary that sits between
-- upstream and downstream FluidBlocks, as identified in the grid-preparation stage.
--
-- Initially it provides a reflecting-wall BC for each of the blocks,
-- holding the test-gas and acceleration-gas slugs in place.
-- This is done by overwriting the previously-exchanged ghost-cell data.
-- When a high pressure signal arrives at the right end of the shock tube,
-- the diaphragm "ruptures" as indicated by the supervisor function.
-- These boundary-condition functions then stop overwriting the exchanged data,
-- such that it appears that the diaphragm has been instantly removed.
--
-- PJ, 2019-06-14, build for David Gildfind's simulations.
--     2024-07-29, Langley expansion tube.

print("Initialize UDF for diaphragm for blkId=", blkId)
dofile('./blkIds.lua')
print("Block upstream of diaphragm=", upstreamBlk, "downstream=", downstreamBlk)
blk_up = infoFluidBlock(upstreamBlk)
blk_down = infoFluidBlock(downstreamBlk)

function ghostCells(args)
   local ruptureFlag = userPad[1]
   local cell0, cell1, cell2
   if blkId == downstreamBlk and args.boundaryId == west then
      if ruptureFlag >= 1 then
	 -- Leave previous exchange data from upstream-block, east boundary.
	 cell0 = {}
	 cell1 = {}
	 cell2 = {}
      else
 	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(downstreamBlk, blk_down.imin, args.j, blk_down.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(downstreamBlk, blk_down.imin+1, args.j, blk_down.kmin)
	 cell1.velx = -cell1.velx
	 cell2 = sampleFluidCell(downstreamBlk, blk_down.imin+2, args.j, blk_down.kmin)
	 cell2.velx = -cell2.velx
     end
   elseif blkId == upstreamBlk and args.boundaryId == east then
      if ruptureFlag >= 1 then
	 -- Leave previous exchange data from downstream-block, west boundary.
	 cell0 = {}
	 cell1 = {}
         cell2 = {}
      else
	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(upstreamBlk, blk_up.imax, args.j, blk_up.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(upstreamBlk, blk_up.imax-1, args.j, blk_up.kmin)
	 cell1.velx = -cell1.velx
	 cell2 = sampleFluidCell(upstreamBlk, blk_up.imax-2, args.j, blk_up.kmin)
	 cell2.velx = -cell2.velx
      end
   else
      error(string.format("No ghost cell b.c. function available"..
			  " for block %d, on boundary: %d", blkId, args.boundaryId))
   end
   return cell0, cell1, cell2
end

function interface(args)
   return sampleFluidCell(blkId, args.i, args.j, args.k)
end
