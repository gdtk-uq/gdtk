-- Author: Ingo JH Jahn
-- Date: 2017-01-04
--
-- This script is used to assign velocities to all of
-- of the vertices in the grid as a function of time.
-- At each time step, we'll do the following:
-- 1. Load the file that was created in "udf-process.lua".
-- 2. Assign the velocity to the corners of the grid using
--    setVtxVelocitiesByCorners(blkId, p0vel, p1vel, p2vel, p3vel)
--    grid internal velocities are then calculated 
--    automatically

--

-- Read in the config parameters as global variables.
dofile('sim-config.lua')


function assignVtxVelocities(sim_time, dt)
   -- Calculate current wall velcoity
   ux = 0
   --ux = u_max * math.cos( omega * sim_time)

   -- Assign velocities to entire block
   setVtxVelocitiesForBlockXYZ(ux, 0., 0., 0)

end
