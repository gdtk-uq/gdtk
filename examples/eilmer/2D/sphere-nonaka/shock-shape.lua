-- shock-shape.lua
--
-- Author: Rowan J. Gollan
-- Date: 2018-01-04
-- 
-- In this shock-fitted solution, we'll take the edge
-- of the grid as the shock location.

local atan2 = math.atan2
local deg = math.deg
local sqrt = math.sqrt

print("Begin shock shape extraction.")

config.grid_motion = "shock_fitting"
jobName = "nonaka"
Db = 14.0e-3
R = Db/2.0
nyVtxsPerBlock = 31
-- Pick up flow solution at final time
fsol = FlowSolution:new{jobName=jobName, dir=".", tindx="last", nBlocks=4}
f = io.open("shock-shape.dat", 'w')
f:write("#  x                   y                    theta                d/R\n")
for ib=0,3 do
   jstart = 1
   if ib == 0 then
      jstart = 0
   end
   for j=jstart,nyVtxsPerBlock-1 do
      vtx = fsol:get_vtx{ib=ib, i=0, j=j}
      x = vtx.x - R
      y = vtx.y
      theta = deg(atan2(y, -x))
      d = sqrt(x*x + y*y) - R
      d_R = d/R
      f:write(string.format("%20.12e %20.12e %20.12e %20.12e\n", x, y, theta, d_R))
   end
end
f:close()

print("Done.")



