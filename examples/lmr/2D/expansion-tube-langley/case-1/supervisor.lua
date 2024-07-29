-- supervisor.lua
-- Authors: PeterJ
-- Date: 2019-06-14, 2019-10-17, 2024-07-29
--
-- Determine the FluidBlock against the upstream side of the diaphragm.
dofile('./blkIds.lua')
print("Block against diaphragm:", upstreamBlk)
blk = infoFluidBlock(upstreamBlk)
pRupture = 2*3.45e3

function atTimestepStart(sim_time, steps, dt)
   if is_master_task then
      -- We are in the shock tube and can see the pressure wave arrive
      -- at the upstream-side of the diaphragm.
      local testCell = sampleFluidCell(upstreamBlk, blk.imax, blk.jmin, blk.kmin)
      -- We store ruptureFlag as a floating-point number, 0 or 1.
      local ruptureFlag = userPad[1]
      if ruptureFlag < 1 and testCell.p > pRupture then
         ruptureFlag = 1
         print("-- Diaphragm rupture at t=", sim_time)
         print("-- pressure=", testCell.p)
      end
      userPad[1] = ruptureFlag
   end
   return
end
