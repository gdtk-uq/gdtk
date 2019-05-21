-- Authors: Rowan J. Gollan, Fabian Zander, Peter J. Ingo J.
-- Date: 2017-01-04 -- 2019-05-21
--
function assignVtxVelocities(t, dt)
   local xdot = userPad[2]
   local zeroVel = Vector3:new{x=0, y=0}
   local pSpeedVec = Vector3:new{x=xdot, y=0}
   local blkId = 0
   setVtxVelocitiesByCorners(blkId, zeroVel, pSpeedVec, pSpeedVec, zeroVel)
end
