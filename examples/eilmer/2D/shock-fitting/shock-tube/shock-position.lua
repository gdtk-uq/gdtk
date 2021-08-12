-- shock-position.lua
-- Use the edge of the grid as the shock location.
-- $ e4shared --custom-post --script-file=shock-position.lua
-- PJ 2021-08-12
print("Begin to locate shock.")
fsol = FlowSolution:new{jobName="sodsf", dir=".", tindx="last", nBlocks=2}
vtx = fsol:get_vtx{ib=0, i=0, j=0}
print("location=", vtx.x)
