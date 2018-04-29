-- Adoptation of 2-D piston-w-const-vel test case to evaluate 
-- Eilmer 4 moving mesh capabilties in 3-D and in multiple directions.
--
-- Author: Ingo Jahn
-- Date: 2018-04-29
--
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

function assignVtxVelocities(sim_time, dt)
   -- Compute present position of piston
   pPos = pSpeed * sim_time

   -- Compute current length of domain
   if Direction =='X' then
      L = Lx - pPos
   elseif Direction =='Y' then 
      L = Ly - pPos
   elseif Direction =='Z' then 
      L = Lz - pPos
   else
      print("Only Direction = 'X' and 'Y' and 'Z' is supported.")
   end

   -- Loop over all cells, assigning vertex velocities
   blk = 0
   imin = blockData[0].vtxImin
   imax = blockData[0].vtxImax
   jmin = blockData[0].vtxJmin
   jmax = blockData[0].vtxJmax
   kmin = blockData[0].vtxKmin
   kmax = blockData[0].vtxKmax
   for j=jmin,jmax do
      for i=imin,imax do
         for k=kmin,kmax do
	    -- Find position in duct
	    pos = getVtxPosition(blk, i, j, k)
	    -- Linearly scale vertex speed based on how far the vertex is from the
	    -- fixed end of the duct.
            if Direction =='X' then
	       vtxSpeedx = ((Lx - pos.x)/L)*pSpeed
               vtxSpeedy = 0.
               vtxSpeedz = 0.
            elseif Direction =='Y' then 
	       vtxSpeedx = 0.
               vtxSpeedy = ((Ly - pos.y)/L)*pSpeed
               vtxSpeedz = 0.
            elseif Direction =='Z' then 
	       vtxSpeedx = 0.
               vtxSpeedy = 0.
               vtxSpeedz = ((Lz - pos.z)/L)*pSpeed
            end
	    -- Set vertex as Vector3 object
	    setVtxVelocity(Vector3:new{x=vtxSpeedx, y=vtxSpeedy, z=vtxSpeedz}, 0, i, j, k)
         end
      end
   end
end
