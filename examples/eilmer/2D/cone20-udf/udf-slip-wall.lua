-- udf-slip-wall.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedBC boundary condition.
--
-- This particular example is defining the slip-wall condition
-- for the cone20 test case.
--
-- RG & PJ 2015-03-14: ported from Eilme3 example

function reflectNormalVelocity(ux, vy, cosX, cosY)
   -- Copied from cns_bc.h.
   un = ux * cosX + vy * cosY;     -- Normal velocity
   vt = -ux * cosY + vy * cosX;    -- Tangential velocity
   un = -un;                       -- Reflect normal component
   ux = un * cosX - vt * cosY;     -- Back to Cartesian coords
   vy = un * cosY + vt * cosX;
   return ux, vy
end

function ghostCells(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   --
   -- args contains t, dt, timeStep, timeLevel, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   i = args.i; j = args.j; k = args.k
   cell1 = sampleFlow(blkId, i, j, k)
   cell1.velx, cell1.vely = reflectNormalVelocity(cell1.velx, cell1.vely, args.csX, args.csY)
   if whichBoundary == north then
      j = j - 1
   elseif whichBoundary == east then
      i = i - 1
   elseif whichBoundary == south then
      j = j + 1
   elseif whichBoundary == west then
      i = i + 1
   end
   cell2 = sampleFlow(blkBd, i, j, k)
   cell2.velx, cell2.vely = reflectNormalVelocity(cell2.velx, cell2.vely, args.csX, args.csY)
   return cell1, cell2
end


