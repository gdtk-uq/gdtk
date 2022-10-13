-- udf supervisor file for moving square (3DOF)
-- idea will be to prescribe a rotation while keeping the body location the same....

function omega(sim_time)
   -- We are imposing a trig function on the motion of the square
    
   max_time = 2e-3
   -- period = 2 pi / B
	
	A = 2000 -- rad/s (initial rotational vel)

   return 2000 * math.cos(2 * math.pi / max_time * sim_time)
end


function atTimestepStart(sim_time, steps, delta_t)

   --[[
        
        This function can be used to compute calcs such as integrating forces
        at each time step start etc. Or really, whatever we want.
        
        For this example, i'll mainly be using it to measure things and write to
        a file, since our dynamics is described.
        
        General Notation:
       
        b - body
        g - global
        
        xb_g --> x position of body, wrt global coordinates.
        
        note: using global coords for pos, vel etc as we dont need body coords
        
   --]]

   -- ==============================================================
   -- Use our userpad values for plotting
   -- ==============================================================
    
   thetab_g = userPad[1] -- rad   
    
   -- ==============================================================
   -- Forces 
   -- ==============================================================

   Fb_g,Mc_g = getRunTimeLoads('walls') -- body forces, moments (global frame)
    
   -- Newtonian Mechanics
    
   Fb_gx = Fb_g.x
   Fb_gy = Fb_g.y  -- No gravity
     
   -- here's where we could integrate forces etc.
    
   -- Eularian Mechanics
    
   Mc_gz = Mc_g.z
   Mb_gz = Mc_gz -- No deviation from {0,0,0} for the body centre...
      
   -- ==============================================================
   -- Integrate relevant states and update to userPad
   -- ==============================================================
    
   omegab_g = omega(sim_time)    
       
   thetab_g_new = thetab_g + omegab_g * delta_t
    
   userPad[1] = thetab_g_new
   userPad[2] = omegab_g
   userPad[3] = steps
    
   if (steps % 2) == 0 then

      if in_mpi_context and not is_master_task then return end
	
	   if (steps % 500) == 0 then
         userPad_file = io.open("square_userpad_data.dat", "a")
        	userPad_file:write(string.format('%.18e %.18e %.18e\n',
                  userPad[1],userPad[2],userPad[3]))
        	userPad_file:close()
	   end
   -- Save data for postprocessing.

      sFname = 'square_data.csv'
      local f = io.open(sFname, "r")
      if not f then
          -- File does not already exist, so we can write header line.
          f = io.open(sFname, "w")
          f:write("time(s) thetab omegab Fxb Fyb Mzbg\n")
          f:close()
      end
      f = io.open(sFname, "a")
      f:write(string.format('%.18e %.18e %.18e %.18e %.18e %.18e\n',
              sim_time, userPad[1], userPad[2], Fb_gx, Fb_gy, Mb_gz))
      f:close()
   end
   return
end

