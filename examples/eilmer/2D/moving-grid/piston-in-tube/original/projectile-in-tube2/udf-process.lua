-- Author: Rowan Gollan / Ingo Jahn / Peter Jacobs / Fabian Zander
-- Data: 2019-05-03
--
-- This particular UDF is a culmination of on-going work by all the
-- authors named above.
--
-- NOTES:
-- Within this script we are demonstrating the use of two significant
-- functions for UDF processing
--
-- atTimestepStart(sim_time, steps, dt):
--    This function is where we do the calculations of the piston
--    piston movement and write out some data to the screen and
--    to data files.
--    Here we demonstrate the use of
--    1. userPad: We utilise the userPad to store our data and then
--       transmit it to our other lua scripts which require this data.
--    2. getRunTimeLoads: The run time loads return the forces and
--       moments imparted on the surfaces by the flow. The surfaces
--       computed are defined in our prep script.
--
-- atWriteToFile(sim_time, steps):
--    This function runs whenever the code is writing a solution. This
--    allows us to write additional data, or do whatever is required
--    at the same time. This can be very useful, in particular for
--    restarting moving mesh simulations. In this case we are writing
--    out the same data twice. Once we are just saving the userPad data
--    as it is, however, this gets overwritten each time the function
--    is called. Secondly, we are saving the same data accumulatively
--    in a different file so that we have the piston position and velocity
--    corresponding to each flow solution from the main solver. These
--    data sets allow us to restart the simulation from an arbitrary time
--    step if required.
-- 
--  NOTE: All the file writing steps are only done by the master task.
--        This prevents confusion and undesired overwriting of data. This
--        is controlled by the code
--        if in_mpi_context and is_master_task then
--
-- Initialise variables
local x       = 0.
local xdot    = 0.
local xdotdot = 0.

local upstreamForce = 0.
local downstreamForce = 0.

-- Load configuration data
dofile('sim-config.lua')

function atTimestepStart(sim_time, steps, delta_t)
  -- Read in table data
  x = userPad[1]
  xdot = userPad[2]
  -- Grab the forces using getRunTimeLoads
  upstreamForce, upstreamMoment = getRunTimeLoads("pistonUpstream")
  downstreamForce, downstreamMoment = getRunTimeLoads("pistonDownstream")
  -- solve momentum equation to get acceleration of the piston 
  xdotdot = ((upstreamForce.x + downstreamForce.x)*2*math.pi) / pMass -- 3.63N friction force
  x       = x + xdot * delta_t --+ 0.5 * xdotdot * delta_t * delta_t
  xdot    = xdot + xdotdot * delta_t
  -- save data to userPad for vtxSpeed Assignment in grid-motion
  userPad[1] = x
  userPad[2] = xdot
  -- The next bit should only happen once, so when the master task is running
  if in_mpi_context and is_master_task then
  -- print progress to screen
    if (steps % 500) ==0 then
      print('++++++++++++++++')
      print('Steps: ', steps, '   sim-time (sec): ',sim_time)
      print(string.format('Piston pushing force %.4f (N) Piston braking force %.4f (N)', upstreamForce.x*2*math.pi, downstreamForce.x*2*math.pi) )
      print(string.format('X-position %.4f (m) Speed %.4f (m/s) Acceleration %.4f (m/s^2)', x, xdot, xdotdot) )
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
        os.remove(outfile_name) -- now delete file
      end

      file=io.open(outfile_name,"a")
      file:write('# pos1:sim_time(s) pos2:Position(m) pos3:Velocity(m/s) \n')
      file:close()
    end
    if (steps % 500)==0 then
      file=io.open(outfile_name,"a")
      --file:write(steps, ' ', sim_time, ' ', P_mean, ' ', x, ' ', xdot, ' ', Fx, '\n')
      file:write(string.format('%.18e %.18e %.18e \n', sim_time, x, xdot) )   
      file:close()          
    end
  end

  return nil 
end

function atWriteToFile(sim_time, steps)
  -- We want to save our userPad data every time we write out a flow solution
  -- Only let the master task do this
  if in_mpi_context and not is_master_task then return end
  -- I want to save my userPad for restarts
  assert( table.save(userPad, "userPadSave.lua") == nil )
  -- I also want to save my piston data at these points to allow a proper energy balance
  file=io.open(piston_data_flow_solution_times,"a")
  file:write(string.format('%.18e %.18e %.18e \n', sim_time, userPad[1], userPad[2]) )   
  file:close()
end