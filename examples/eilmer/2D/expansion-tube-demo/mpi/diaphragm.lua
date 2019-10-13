-- diaphragm.lua
--
-- These user-defined functions supplement the Exchange boundary condition
-- that would have been applied to the diaphragm boundary that sits between
-- FluidBlocks [1] and [2].
-- Initially it provides a reflecting-wall BC for each of the blocks, holding 
-- This is done by overwriting the previously-exchanged ghost-cell data.
-- When a high pressure signal arrives at the right end of FluidBlock [1],
-- the diaphragm "ruptures" as indicated by the supervisor function.
-- These boundary-condition functions then stop overwriting the exchanged data,
-- such that it appears that the diaphragm has been instantly removed.
--
-- PJ, 2019-06-14, build for David Gildfind's simulations.

print("Initialize UDF for diaphragm for blkId=", blkId)
blk1 = infoFluidBlock(1)
blk2 = infoFluidBlock(2)

function ghostCells(args)
   -- Notes
   -- (1) Because we are using sampleFluidCell to read the cell state,
   -- we should not apply these boundary conditions in parallel.
   -- (2) This ghostCells function is defined at BC initialization and
   -- later called by the UserDefinedBC at the start of each phase
   -- of the time step.
   --
   --
   local ruptureFlag = userPad[1]
   local cell0, cell1
   if blkId == 2 and args.boundaryId == west then
      if ruptureFlag >= 1 then
	 -- Leave previous exchange data from block 1, east boundary.
	 cell0 = {}
	 cell1 = {}
      else
 	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(2, blk2.imin, args.j, blk2.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(2, blk2.imin+1, args.j, blk2.kmin)
	 cell1.velx = -cell1.velx
     end
   elseif blkId == 1 and args.boundaryId == east then
      if ruptureFlag >= 1 then
	 -- Leave previous exchange data from block 2, west boundary.
	 cell0 = {}
	 cell1 = {}
      else
	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(1, blk1.imax, args.j, blk1.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(1, blk1.imax-1, args.j, blk1.kmin)
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

