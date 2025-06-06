-- shock-standoff.lua
-- Invoke with the command line:
-- $ lmr custom-script --job=shock-standoff.lua
-- Author: Nick GIbbons
-- Date: 2025-06-05

R = 31.8e-3  -- radius of sphere, in metres
fsol = FlowSolution:new{dir=".", snapshot="last", nBlocks=4}

vtx0 = fsol:get_vtx{ib=0, i=0, j=0}
ss = -1.0*(vtx0.x + R)

print(string.format("Shock Standoff: %4.3f mm", ss*1000.0))



