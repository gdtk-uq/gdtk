-- shock-position.lua
-- Use the edge of the grid as the shock location.
-- $ e4shared --custom-post --script-file=shock-position.lua
-- PJ 2021-08-14
print("Begin to locate shock.")
fsol = FlowSolution:new{jobName="sco2", dir=".", tindx="last", nBlocks=4}
vtx = fsol:get_vtx{ib=0, i=0, j=0}
print("location=", vtx.x)
R = 0.005 -- metres
delta = -(R+vtx.x)
print("delta=", delta)
print("delta/R=", delta/R)
