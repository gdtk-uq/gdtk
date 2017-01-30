-- Author: Rowan J. Gollan
-- Date: 2017-01-04
--
-- This script is used to assign velocities to all of
-- of the vertices in the grid as a function of time.
-- At each time step, we'll do the following:
-- 1. Compute where the piston face is.
-- 2. Compute the length of the gas domain, given that the
--    piston has compressed some of the domain.
-- 3. Assign velocities to vertices in a linearly scaled manner.
--    The vertices at the piston face will have the velocity
--    of the piston. The vertices at the fixed wall will have
--    zero velocity. All the vertices in between will have a 
--    a velocity that is scaled based on how close they are
--    to the end wall. Vertices close to the end wall have
--    low velocities, while vertices close to the piston face
--    receive high velocities. The main 'trick' in the function
--    below is to loop over the vertices asking each what its
--    position is. Based on that information, we can compute
--    a scaling of the velocity to assign to the vertex.
--
-- NOTES:
-- We make use of some global data and functions that Eilmer
-- provides use for doing the work. In particular, we use:
-- 1. blockData :
--    By digging into blockData[0] to get at the data in
--    block[0], we can find out the indices for the min and
--    and max in the logical coordinates of 'i' and 'j'.
-- 2. getVtxPosition:
--    This function takes as arguments four integers:
--    blkId, i, j, k
--    and uses these to locate vertex in the grid.
--    The function then returns the position of the
--    vertex in physical space (Cartesian coordinates)
--    as a Vector3 object.
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
endDomain = L

function assignVtxVelocities(sim_time)
   -- Compute present position of piston
   pPos = pSpeed * sim_time
   -- Compute current length of domain
   L = endDomain - pPos
   -- Loop over all cells, assigning vertex velocities
   imin = blockData[0].vtxImin
   imax = blockData[0].vtxImax
   jmin = blockData[0].vtxJmin
   jmax = blockData[0].vtxJmax
   for j=jmin,jmax do
      for i=imin,imax do
	 -- Find position in duct
	 pos = getVtxPosition(0, i, j, 0)
	 -- Linearly scale vertex speed based on how far the vertex is from the
	 -- fixed end of the duct.
	 vtxSpeed = ((endDomain - pos.x)/L)*pSpeed
	 -- Set vertex as Vector3 object
	 setVtxVelocity(Vector3:new{x=vtxSpeed}, 0, i, j)
      end
   end
end
