-- supervisor.lua
-- Authors: PeterJ
-- Date: 2019-06-14
--

blk1 = infoFluidBlock(1)

function atTimestepStart(sim_time, steps, dt)
   if is_master_task then
      -- We are in the shock tube and can see the pressure wave arrive
      -- at the upstream-side of the diaphragm.
      local testCell = sampleFluidCell(1, blk1.imax, blk1.jmin, blk1.kmin)
      -- We store ruptureFlag as a floating-point number, 0 or 1.
      local ruptureFlag = userPad[1]
      if ruptureFlag < 1 and testCell.p > 1.0e5 then
         ruptureFlag = 1
         print("-- Diaphragm rupture at t=", sim_time)
      end
      userPad[1] = ruptureFlag
   end
   return
end
