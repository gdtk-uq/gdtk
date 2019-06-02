-- Authors: Rowan J. Gollan & Peter J.
-- Dates: 2017-01-04 -- 2019-06-02
-- Test the setting of vertex velocities by the more general call.
--
zero = Vector3:new{x=0, y=0}
pSpeed = Vector3:new{x=293.5, y=0}

function assignVtxVelocities(t, dt)
   local blkId = 0
   blk = infoFluidBlock(blkId)
   p00 = getVtxPositionVector3(blkId, blk.imin, blk.jmin, blk.kmin)
   p10 = getVtxPositionVector3(blkId, blk.imax+1, blk.jmin, blk.kmin)
   p11 = getVtxPositionVector3(blkId, blk.imax+1, blk.jmax+1, blk.kmin)
   p01 = getVtxPositionVector3(blkId, blk.imin, blk.jmax+1, blk.kmin)
   setVtxVelocitiesByQuad{blkId=blkId, p00=p00, p10=p10, p11=p11, p01=p01,
                          vel00=pSpeed, vel10=zero, vel11=zero, vel01=pSpeed}
end
