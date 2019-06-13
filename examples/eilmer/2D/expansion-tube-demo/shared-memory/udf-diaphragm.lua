-- udf-diaphragm.lua
--
-- This diaphragm-like boundary condition connects FluidBlocks [1] and [2].
-- Initially it provides a reflecting-wall BC for each of the blocks, holding 
-- the gases in place.  When a high pressure signal arrives at the right end 
-- of FluidBlock [1], the diaphragm "ruptures".  The boundary condition
-- then changes to exchanging the flow data between the blocks, such that it 
-- appears that the diaphragm has been instantly removed.
--
-- PJ, 2017-05-07: Built for Kristin Stewart's thesis.

print("Initialize UDF for diaphragm.")
ruptureFlag = false
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
   -- Check for arrival of shock at end of shock tube.
   local testCell = sampleFluidCell(1, blk1.imax, blk1.jmin, blk1.kmin)
   if not ruptureFlag and testCell.p > 1.0e5 then
      ruptureFlag = true
   end
   --
   local cell0, cell1
   if blkId == 2 and args.boundaryId == west then
      if ruptureFlag then
	 -- Copy flow data from block 1, east boundary.
	 cell0 = sampleFluidCell(1, blk1.imax, args.j, blk1.kmin)
	 cell1 = sampleFluidCell(1, blk1.imax-1, args.j, blk1.kmin)
      else
 	 -- A simple reflecting boundary condition.
	 cell0 = sampleFluidCell(2, blk2.imin, args.j, blk2.kmin)
	 cell0.velx = -cell0.velx
	 cell1 = sampleFluidCell(2, blk2.imin+1, args.j, blk2.kmin)
	 cell1.velx = -cell1.velx
     end
   elseif blkId == 1 and args.boundaryId == east then
      if ruptureFlag then
	 -- Copy flow data from block 2, west boundary.
	 cell0 = sampleFluidCell(2, blk2.imin, args.j, blk2.kmin)
	 cell1 = sampleFluidCell(2, blk2.imin+1, args.j, blk2.kmin)
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


