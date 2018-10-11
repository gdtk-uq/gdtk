-- Author: Ingo HJ Jahn
-- Date: 2018-10-11
--
-- A user-defined b.c. to specify constant pressure at top but and no shear
require 'lua_helper' 

-- Read in some config parameters and set as global variables
dofile('sim-config.lua')

--[[
function convective_flux(args)
   -- We need to get the pressure at the cell adjacent to the boundary.
   -- We aren't going to do any fancy reconstruction of the piston face pressure,
   -- we are simply going to take the adjacent cell pressure.
   cell = sampleFlow(blkId, args.i, args.j, args.k)
   p = cell.p
   flux = {}
   flux.mass = 0
   flux.momentum_x = p
   flux.momentum_y = 0.0
   flux.total_energy = 0
   return flux   
end
--]]

function ghostCells(args)
   -- A simple outflow boundary condition can be implemented by looking at
   -- the flow condition for the cell just inside the boundary.
   -- Notes:
   -- (1) This shows the use of some of the global data set in the Lua
   -- interpreter for this boundary condition, as well as the use of
   -- some of the elements of the args table that is given to this function.
   -- (2) Because we are using sampleFlow to read the cell state,
   -- we should apply these boundary conditions serially,
   -- to avoid a race condition.
   cell = sampleFlow(blkId, args.i, args.j, args.k)
   return cell, cell
end


function interface(args)
   -- We need to get the pressure at the cell adjacent to the boundary.
   -- We aren't going to do any fancy reconstruction of the piston face pressure,
   -- we are simply going to take the adjacent cell pressure.
   cell = sampleFlow(blkId, args.i, args.j, args.k)
   p = cell.p
   T = cell.T
   inter = {}
   inter.p = p
   inter.T = startT
   --inter.T = T
   inter.velx = u_max * math.cos( omega * args.t)   
   return inter   
end


