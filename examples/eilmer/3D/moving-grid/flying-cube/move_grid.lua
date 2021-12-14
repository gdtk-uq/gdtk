-- Cube flying mesh movement for TUSQ simulation
--
-- Author: Fabian Zander
-- Last Modified: 26/08/2021

function assignVtxVelocities(sim_time, dt)
  -- I'm doing a dodgy on the movement here - one us is translation,
  -- the next us is rotation. Hence double the values for each us.
  
  -- I need a time which can use the modulo function
  tmpTime = math.floor(sim_time*1e6)

  if (tmpTime % 2 == 0) then
    -- First we translate
    xdot = userPad[4]*2 -- Multiply by 2 for every other us
    ydot = userPad[6]*2 -- Multiply by 2 for every other us
    zdot = userPad[8]*2
    -- Noting the use of the Domain movement here
    setVtxVelocitiesForDomain(Vector3:new{x=xdot, y=ydot, z=zdot})

  else
    -- Then we rotate
    alphadot = userPad[2]*2 -- Multiply by 2 for every other us
    x = userPad[3]
    y = userPad[5]
    z = userPad[7]
    for _,blkId in ipairs(localFluidBlockIds) do
      -- This time we rotate each block around the current 'centre' location
      setVtxVelocitiesForRotatingBlock(blkId, alphadot, Vector3:new{x=x, y=y, z=z})
    end
  end
end
