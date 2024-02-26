
-- Read in grid:
for blockfiles in io.popen("ls -1 ./su2grid | wc -l"):lines() do
   nblocks = tonumber(blockfiles)
end

for i=0,nblocks-1 do
   fileName = string.format("su2grid/block_%d_grid.su2", i)
   grid = registerGrid{
      grid = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0},
      fsTag="inflow",
   }
end
