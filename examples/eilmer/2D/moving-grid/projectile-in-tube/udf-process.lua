
require 'lua_helper' 

-- Initialise variables
local Fx      = 0.
local x       = 0.
local xdot    = 0.
local xdotdot = 0.

local delta_t = 0.
local T_last  = 0.
local P_mean0 = 0.
local P_mean1 = 0.

function atTimestepStart (sim_time,steps,gasBlocks)

   -- Load configuration data
   dofile('sim-config.lua') 

   -- calculate time passed since last velocity update
   delta_t =    sim_time-T_last
   T_last  =    sim_time

   Fx = 0.
   blkId = 0
   k = 0
   P_mean0  = 0.
   P_mean1  = 0.   

   imin = blockData[blkId].vtxImin
   imax = blockData[blkId].vtxImax
   jmin = blockData[blkId].vtxJmin
   jmax = blockData[blkId].vtxJmax

   for j=jmin,(jmax-1) do
       pos0 = getVtxPosition(blkId, imax, j, k)
       pos1 = getVtxPosition(blkId, imax, j+1, k)

       -- get conditions at end of tube (nearest to projectile)
       Dat  = sampleFlow(blkId,imax-1,jmin,k)
       Fx   = Fx + Dat.p * (pos1.y - pos0.y)
       P_mean1 = P_mean1 + Dat.p * (pos1.y - pos0.y)

       -- get conditions at start of tube
       Dat  = sampleFlow(blkId,imin,jmin,k)
       P_mean0 = P_mean0 + Dat.p * (pos1.y - pos0.y)
   end
   pos0 = getVtxPosition(blkId, imax, jmin, k)
   pos1 = getVtxPosition(blkId, imax, jmax, k) 
   P_mean0 = P_mean0 / ( pos1.y - pos0.y )
   P_mean1 = P_mean1 / ( pos1.y - pos0.y )

   Dat0  = sampleFlow(blkId,imax-1,jmin,k)
   Dat1  = sampleFlow(blkId,imax-1,jmin,k)
   Dat2  = sampleFlow(blkId,imax-2,jmin,k)
   --print("P, imax-1:", Dat0.p, Dat1.p, Dat2.p)


   -- solve momentum equation to get acceleration of projectile 
   xdotdot = Fx / massPiston
   xdot    = xdot + xdotdot * delta_t
   x       = x + xdot * delta_t + 0.5 * xdotdot * delta_t * delta_t 

   -- save data to file for vtxSpeed Assignment in grid-motion
   TABLE = {}
   TABLE[0] = x
   TABLE[1] = xdot
   assert( table.save(TABLE, "data.lua") == nil )

   -- print progress to screen
   if (steps % 500) ==0 then
      print('++++++++++++++++')
      print('Steps: ', steps, '   sim-time (sec): ',sim_time)
      print(string.format('X-position %.4f (m) Speed %.4f (m/s) Acceleration %.4f (m/s^2) Pressure %.4f (Pa)', x, xdot, xdotdot, P_mean1) )    
   end

   -- write data to file
   if steps==0 then
      -- check if file exists
      local infile=io.open(outfile_name,"r")
      if infile~=nil then -- yes file already exists
         instr = infile:read("*a")
         infile:close() 
         -- write contents to a new file
         outfile = io.open(outfile_name .. ".old", "w")
         outfile:write(instr)
         outfile:close()
         print('Output file=' .. outfile_name .. ' already exists. File copied to ' .. outfile_name .. '.old')
         os.remove(outfile_name) -- now dellete file
      end

      file=io.open(outfile_name,"a")
      file:write('# pos1:sim_time(s) pos2:Pressure0(Pa) pos3:Pressure1(Pa) pos4:Position(m) pos5:Velocity(m/s) pos6:TotalForce(N/m) \n')
      file:close()
   end
   if (steps % 500)==0 then
      file=io.open(outfile_name,"a")
      --file:write(steps, ' ', sim_time, ' ', P_mean, ' ', x, ' ', xdot, ' ', Fx, '\n')
      file:write(string.format('%.18e %.18e %.18e %.18e %.18e %.18e \n', sim_time, P_mean0, P_mean1, x, xdot, Fx) )   
      file:close()	        
   end

  return nil 
end
