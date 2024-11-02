-- udf-supervisor.lua
-- Authors: Rowan Gollan / Ingo Jahn / Peter Jacobs / Fabian Zander
-- Date: 2019-05-03 -- 2019-05-17
--
dofile('sim-config.lua')

function atTimestepStart(sim_time, steps, dt)
   -- Unpack current piston state.
   local x = userPad[1]
   local xdot = userPad[2]
   -- Get the surface loads on the piston using getRunTimeLoads.
   local uForce, uMoment = getRunTimeLoads("pistonUpstream")
   local dForce, dMoment = getRunTimeLoads("pistonDownstream")
   -- Acceleration of the piston.
   local xdotdot = ((uForce.x + dForce.x)*2*math.pi) / pMass
   -- Update piston state using semi-implicit Euler update.
   xdot = xdot + xdotdot * dt
   x    = x + xdot * dt
   -- Save data to userPad for vtxSpeed assignment in grid-motion.
   userPad[1] = x
   userPad[2] = xdot
   return
end

pFileName = "piston.data"

function atWriteToFile(sim_time, steps)
   -- We want to save our userPad data every time we write out a flow solution
   -- but only let the master task do this.
   if in_mpi_context and not is_master_task then return end
   -- Save piston data for postprocessing.
   local f = io.open(pFileName, "r")
   if not f then
      -- File does not already exist, so we can write header line.
      f = io.open(pFileName, "w")
      f:write("#   time(s)      x(m)      v(m/s)\n")
      f:close()
   end
   f = io.open(pFileName, "a")
   f:write(string.format('%.6e %.6e %.6e \n', sim_time, userPad[1], userPad[2]) )
   f:close()
   return
end
