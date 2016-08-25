-- shock-detachnment.lua
--
-- Author: Rowan J. Gollan
-- Date: 2016-01-23
-- 
-- In this shock-fitted solution, we'll take the position
-- of the bottom left vertex in the grid as the shock
-- detachment distance. So, if we assume that block 0
-- is on the centreline, then all we need to do is get
-- the location of that bottom corner grid point.

config.grid_motion = "shock_fitting"
jobName = "lobb"
Db = 0.5 * 0.0254 -- diameter (in m) of ball bearing
-- Pick up flow solution at final time
fsol = FlowSolution:new{jobName=jobName, dir=".", tindx="last", nBlocks=4}
vtx = fsol:get_vtx{ib=0, i=0, j=0}
delta = -vtx.x
d_D = delta/Db
f = io.open("shock-detachment.txt", 'w')
f:write(string.format("%20.12e %20.12e\n", delta, d_D))
f:close()
print("shock-detachment= ", delta)
print("delta/D= ", d_D)


