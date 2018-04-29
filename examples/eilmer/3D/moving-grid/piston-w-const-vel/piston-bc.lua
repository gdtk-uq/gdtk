-- Author: Rowan J. Gollan
-- Date: 2017-01-04
--
-- A user-defined b.c. to specify the convective flux
-- at a moving piston boundary.

-- Read in some config parameters and set as global variables
dofile('sim-config.lua')

function convective_flux(args)
   -- We need to get the pressure at the cell adjacent to the boundary.
   -- We aren't going to do any fancy reconstruction of the piston face pressure,
   -- we are simply going to take the adjacent cell pressure.
   cell = sampleFlow(blkId, args.i, args.j, args.k)
   p = cell.p
   flux = {}
   flux.mass = 0
   if Direction =='X' then
      flux.momentum_x = p
      flux.momentum_y = 0.0
      flux.momentum_z = 0.0
   elseif Direction =='Y' then
      flux.momentum_x = 0.0
      flux.momentum_y = p
      flux.momentum_z = 0.0
   elseif Direction =='Z' then
      flux.momentum_x = 0.0
      flux.momentum_y = 0.0
      flux.momentum_z = p
   end
   flux.total_energy = p*pSpeed
   return flux   
end

