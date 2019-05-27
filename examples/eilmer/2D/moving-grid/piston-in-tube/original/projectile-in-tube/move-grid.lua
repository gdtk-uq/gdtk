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
   -- Load instantaneous velocity
   TABLE,err = table.load('data.lua')
   x    = TABLE[0]
   xdot = TABLE[1]


   -- Assign velocities at corners of block
   setVtxVelocitiesByCorners(0, Vector3:new{x=0, y=0}, Vector3:new{x=xdot, y=0}, 
                                Vector3:new{x=xdot, y=0}, Vector3:new{x=0, y=0})


end
