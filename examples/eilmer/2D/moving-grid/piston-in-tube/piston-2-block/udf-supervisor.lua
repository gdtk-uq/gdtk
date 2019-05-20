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
   local upstreamForce, upstreamMoment = getRunTimeLoads("pistonUpstream")
   local downstreamForce, downstreamMoment = getRunTimeLoads("pistonDownstream")
   -- Acceleration of the piston.
   local xdotdot = ((upstreamForce.x + downstreamForce.x)*2*math.pi) / pMass
   -- Update piston state using simple Euler update.
   x    = x + xdot * dt
   xdot = xdot + xdotdot * dt
   -- Save data to userPad for vtxSpeed Assignment in grid-motion.
   userPad[1] = x
   userPad[2] = xdot
   return
end

pFileName = "piston.data"

function atWriteToFile(sim_time, steps)
   -- We want to save our userPad data every time we write out a flow solution
   -- but only let the master task do this.
   if in_mpi_context and not is_master_task then return end
   --
   -- Save userPad for possible restart of simulation.
   assert( table.save(userPad, "userPadSave.lua") == nil )
   -- Save piston data for postprocessing.
   local f = io.open(pFileName, "r")
   if not f then
      -- File does not already exist, so we can write header line.
      f = io.open(pFileName, "w")
      f:write("#         time(s)                  x(m)                  v(m/s)\n")   
      f:close()
   end
   f = io.open(pFileName, "a")
   f:write(string.format('%.18e %.18e %.18e \n', sim_time, userPad[1], userPad[2]) )   
   f:close()
   return
end
