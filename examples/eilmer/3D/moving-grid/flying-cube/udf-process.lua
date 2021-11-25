-- Initialise variables
local x       = 0.
local xdot    = 0.
local xdotdot = 0.
local y       = 0.
local ydot    = 0.
local ydotdot = 0.
local z       = 0.
local zdot    = 0.
local zdotdot = 0.

local upstreamForce = 0.
local downstreamForce = 0.

-- Set a file name for saving data
local outfile_name = 'cube_output.dat'

-- Configure the cube parameters
cubeLength = 25.4e-3
cubeMass = 0.1278
cubeI = cubeMass*cubeLength*cubeLength/6

-- At each time step we need to...
function atTimestepStart(sim_time, steps, delta_t)
  
  -- Get the variables we need from userPad
  alpha = userPad[1]
  alphadot = userPad[2]
  x = userPad[3]
  xdot = userPad[4]
  y = userPad[5]
  ydot = userPad[6]
  z = userPad[7]
  zdot = userPad[8]

  -- Grab the forces using getRunTimeLoads
  cubeForce, cubeMoment = getRunTimeLoads("walls")

  -- Solve momentum equation to get acceleration of the cube in both axes
  -- This is just for our data recording, in the simulation, we are only
  -- going to rotate the cube
  alphadotdot = cubeMoment.z / cubeI
  alpha   = alpha + alphadot * delta_t
  alphadot = alphadot + alphadotdot * delta_t
  xdotdot = cubeForce.x / cubeMass
  x       = x + xdot * delta_t
  xdot    = xdot + xdotdot * delta_t
  ydotdot = (cubeForce.y) / cubeMass
  y       = y + ydot * delta_t
  ydot    = ydot + ydotdot * delta_t
  zdotdot = (cubeForce.z) / cubeMass - 9.81
  z       = z + zdot * delta_t
  zdot    = zdot + zdotdot * delta_t

  -- save data to userPad for vtxSpeed Assignment in grid-motion
  userPad[1] = alpha
  userPad[2] = alphadot
  userPad[3] = x
  userPad[4] = xdot
  userPad[5] = y
  userPad[6] = ydot
  userPad[7] = z
  userPad[8] = zdot

  -- The rest is 'just' for outputting and saving data
  if in_mpi_context and not is_master_task then
  else
    -- print progress to screen
    if (steps % 500) == 0 then
      print('++++++++++++++++')
      print('Steps: ', steps, '   sim-time (sec): ',sim_time)
      print(string.format('alpha %.4f (rad) alphadot %.4f (rad/s)', alpha, alphadot) )
      print(string.format('x %.4f (m) xdot %.4f (m/s)', x, xdot) )
      print(string.format('y %.4f (m) ydot %.4f (m/s)', y, ydot) )
      print(string.format('z %.4f (m) zdot %.4f (m/s)', z, zdot) )
    end

    -- write data to file
    if steps == 0 then
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
      file:write('# pos1:sim_time(s) pos2:Alpha(rad) pos3:AlphaDot(rad/s) pos4:PositionX(m) pos5:VelocityX(m/s) pos6:PositionY(m) pos7:VelocityY(m/s) pos8:PositionZ(m) pos3:VelocityZ(m/s) \n')
      file:close()
    end
    if (steps % 500) == 0 then
      file=io.open(outfile_name,"a")
      file:write(string.format('%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e\n', sim_time, alpha, alphadot, x, xdot, y, ydot, z, zdot) )   
      file:close()	        
    end
  end

  return nil 
end

-- The is for saving our userpad data when we write flow data so that we can do a restart
function atWriteToFile(sim_time, steps)
  -- We want to save our userPad data every time we write out a flow solution
  -- Only let the master task do this
  if in_mpi_context and not is_master_task then return end
  -- I want to save my userPad for restarts
  assert( table.save(userPad, "userPadSave.lua") == nil )
  -- I also want to save my cube data at these points to allow an
  -- arbitrary restart
  file=io.open("cube_userpad_data.dat", "a")
  file:write(string.format('%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e\n',
    sim_time, userPad[1], userPad[2], userPad[3], userPad[4], userPad[5], userPad[6], userPad[7], userPad[8]) )   
  file:close()
end
