-- Authors: Rowan J. Gollan, Fabian Zander, Peter J. Ingo J.
-- Date: 2017-01-04 -- 2019-05-17
--
dofile('sim-config.lua')

function assignVtxVelocities(t, dt)
   xdot = userPad[2]
   zeroVel = Vector3:new{x=0, y=0}
   pSpeedVec = Vector3:new{x=xdot, y=0}
   for _,blkId in ipairs(localBlockIds) do
      if blkId == 0 then
         setVtxVelocitiesByCorners(blkId, zeroVel, pSpeedVec, pSpeedVec, zeroVel)
      end
      if blkId == 1 then
         setVtxVelocitiesByCorners(blkId, pSpeedVec, zeroVel, zeroVel, pSpeedVec)
      end
   end
end
