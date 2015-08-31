-- udf-bc-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example contains multiple boundary conditions.
-- They are specialised by determining which boundary is in effect.
--
-- RG & PJ 2015-03-14 : ported from PJs eilmer3 example

function ghostCells_west(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, timeStep, timeLevel, dt, x, y, z, csX, csY, csZ, i, j, k
   -- but we don't happen to us any of them.
   --
   -- Supersonic inflow
   ghost = {}
   ghost.p = 95.84e3 -- pressure, Pa
   ghost.T = {1103.0}  -- temperatures, K (as a table)
   ghost.massf = {1.0} -- mass fractions to be provided as a table
   ghost.velx = 1000.0  -- x-velocity, m/s
   ghost.vely = 0.0     -- y-velocity, m/s
   ghost.velz = 0.0
   return ghost, ghost
end


function ghostCells_east(args)
   -- Because we are using sampleFlow to read the cell state at teh boundary,
   -- we have to apply these boundary conditions serially.
   cell = sampleFlow(blkId, args.i, args.j, args.k)
   return cell, cell
end

