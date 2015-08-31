-- udf-extrapolate-out.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the supersonic outflow
-- for the cone20 test case.
--
-- RG & PJ 2015-03-14: ported from Eilmer3 example

function ghostCells(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   --
   -- args contains t, dt, timeLevel, timeStep, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Sample the flow field at the current cell 
   -- which is beside the boundary.
   cell = sampleFlow(blkId, args.i, args.j, args.k)
   return cell, cell
end

