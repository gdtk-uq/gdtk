-- supervisor.lua
-- Authors: PeterJ
-- Date: 2019-06-14, 2019-10-17
--
-- Determine the FluidBlock against the upstream side of the diaphragm.
dofile('config/exptube.fluidBlockArrays')
local shock_tube = fluidBlockArrays[1]
blk_n = shock_tube.blockArray[shock_tube.nib][1][1]
print("Block against diaphragm:", blk_n)
blk = infoFluidBlock(blk_n)

function atTimestepStart(sim_time, steps, dt)
   if is_master_task then
      -- We are in the shock tube and can see the pressure wave arrive
      -- at the upstream-side of the diaphragm.
      local testCell = sampleFluidCell(blk_n, blk.imax, blk.jmin, blk.kmin)
      -- We store ruptureFlag as a floating-point number, 0 or 1.
      local ruptureFlag = userPad[1]
      if ruptureFlag < 1 and testCell.p > 1.0e4 then
         ruptureFlag = 1
         print("-- Diaphragm rupture at t=", sim_time)
         print("-- pressure=", testCell.p)
      end
      userPad[1] = ruptureFlag
   end
   return
end
