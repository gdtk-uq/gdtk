-- Authors: Rowan J. Gollan & Peter J.
-- Dates: 2017-01-04 -- 2019-05-21
--
zero = Vector3:new{x=0, y=0}
pSpeed = Vector3:new{x=293.5, y=0}

function assignVtxVelocities(t, dt)
   local blkId = 0
   -- Corner labels -               p00     p10   p11   p01
   setVtxVelocitiesByCorners(blkId, pSpeed, zero, zero, pSpeed)
end
