-- udf-bc.lua
-- Lua script for the user-defined boundary conditions.
--
-- The ghostCells_xxxx functions are defined at initialization and
-- later called by the UserDefinedGhostCell boundary condition at 
-- the start of each time step.
--
-- RG & PJ 2015-03-14 : ported from PJs eilmer3 example
--         2016-04-13 : reconfigure with a single ghostCells function
--                      as a delegator function

function ghostCells(args)
   -- Delegate to the appropriate sub-function
   if args.boundaryId == west then
      return ghostCells_west(args)
   elseif args.boundaryId == east then
      return ghostCells_east(args)
   else
      print("No ghost cell b.c. function available on boundary: ", args.boundaryId)
   end 
end

function ghostCells_west(args)
   -- For a supersonic inflow, just set the properties directly.
   ghost = {}
   ghost.p = 95.84e3 -- pressure, Pa
   ghost.T = 1103.0  -- temperatures, K (as a table)
   ghost.massf = {air=1.0} -- mass fractions to be provided as a table
   -- or, if working with a single-species gas model, omitted.
   ghost.velx = 1000.0  -- x-velocity, m/s
   ghost.vely = 0.0     -- y-velocity, m/s
   ghost.velz = 0.0
   return ghost, ghost
end

function ghostCells_east(args)
   -- A simple outflow boundary condition can be implemented by looking at
   -- the flow condition for the cell just inside the boundary.
   -- Notes:
   -- (1) This shows the use of some of the global data set in the Lua
   -- interpreter for this boundary condition, as well as the use of
   -- some of the elements of the args table that is given to this function.
   -- (2) Because we are using sampleFluidCell to read the cell state,
   -- we should apply these boundary conditions serially,
   -- to avoid a race condition.
   cell = sampleFluidCell(blkId, args.i, args.j, args.k)
   return cell, cell
end

