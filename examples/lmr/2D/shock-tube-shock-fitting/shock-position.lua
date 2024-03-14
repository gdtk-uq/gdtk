-- shock-position.lua
-- Use the edge of the grid as the shock location.
-- $ lmr custom-script --job=shock-position.lua
-- PJ 2021-08-12, 2024-03-15
print("Begin to locate shock.")
fsol = FlowSolution:new{dir=".", snapshot="last", nBlocks=2}
vtx = fsol:get_vtx{ib=0, i=0, j=0}
print("location=", vtx.x)
