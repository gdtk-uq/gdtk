-- Author: Rowan J. Gollan
-- Date: 2016-04-22
--
-- History:
--   2018-02-14 : Use Fabian Zander's example with BC labels.
--   2023-11-28 : Update for use in staged preparation

gproGrid = "blk.tmp"
gproConn = "blk.tmp.conn"
gproPty = "blk.tmp.pty"

config.title = "Gridpro importer test"
config.dimensions = 3

-- Read in GridPro grids.
-- The GridPro grid file contains all the blocks in one file.
-- We'll separate these into individual pieces, and store
-- a list of grids.
scale = 1.0
grids = importGridproGrid(gproGrid, scale)

for i,g in ipairs(grids) do
   registerFluidGrid{grid=g, fsTag="initial"}
end

importGridproConnectivity(gproConn)
importGridproBCs(gproPty)

