-- grid.lua
-- Reacting flow over a 15 degree ramp
-- Kyle Damm & Nick Gibbons
-- 2024-04-08


-- Define the flow domain using an imported grid.
for blockfiles in io.popen("ls -1 ./su2grid | wc -l"):lines() do
   nblocks = tonumber(blockfiles)
end

grids = {}
for i=0,nblocks-1 do
   fileName = string.format("su2grid/block_%d_grid.su2", i)
   registerFluidGrid{grid=UnstructuredGrid:new{filename=fileName, fmt="su2text"}, fsTag="inflow"}
end

