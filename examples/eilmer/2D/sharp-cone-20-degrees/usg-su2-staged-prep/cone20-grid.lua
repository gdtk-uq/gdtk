-- Unstructured Grid Example -- for use with Eilmer4
-- 2023-03-22 RJG, PJ, KAD
-- Process with command:
-- $ e4shared --prep-grid --job=cone20

config.dimensions = 2

-- Define the flow domain using an imported grid.
grids = {}
for i=0,3 do
   fileName = string.format("cone20_grid%d.su2", i)
   grids[i] = registerGrid{
      grid=UnstructuredGrid:new{filename=fileName, fmt="su2text"},
      fsTag="inflow-gas"
   }
end
