-- Author: Rowan J. Gollan
-- Date: 2017-01-04
--
-- Modified for the new projectile test case
-- Author: Fabian Zander
-- Data: 2019-05-03
--
-- NOTES:
-- We make use of some global data and functions that Eilmer
-- provides use for doing the work. In particular, we use:
-- 1. userPad:
--    User pad gives us a way of storing data between different interpreters
--    within our simulation, i.e. the velocity data I have calculated in
--    udf-process.lua is available here to move the mesh
-- 2. localFluidBlockIds:
--    This allows us to cycle through the blocks within this particular
--    process. This is required for the MPI implementation where we do not
--    have access to the blocks that are not being computed within this
--    process.
-- 3. setVtxVelocity: 
--    This is the most important function. We need to use this
--    to actually have some effect on the grid motion.
--    This function takes 4 arguments (or 5 in 3D):
--    velVector, blkId, i, j
--    and sets vertex[blkId,i,j] to have the velocity
--    given by velVector.
--

-- Read in the config parameters as global variables.
dofile('sim-config.lua')

function assignVtxVelocities(sim_time, dt)
   --for n in pairs(_G) do print(n) end
   -- Read in the data table with the required info
   xdot = userPad[2]
   zeroVel = Vector3:new{x=0.0, y=0.0}
   pSpeedVec = Vector3:new{x=xdot, y=0.0}
   for _,blkId in ipairs(localFluidBlockIds) do
      if blkId == 0 then
         setVtxVelocitiesByCorners(blkId, zeroVel, pSpeedVec, pSpeedVec, zeroVel)
      end
      if blkId == 1 then
         setVtxVelocitiesByCorners(blkId, pSpeedVec, zeroVel, zeroVel, pSpeedVec)
      end
   end
end
