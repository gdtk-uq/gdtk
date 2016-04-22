-- Author: Rowan J. Gollan
-- Date: 2016-04-22

gproGrid = "blk.tmp"
gproPty = "blk.tmp.pty"

config.title = "Gridpro importer test"
config.dimensions = 3

-- Set up a dummy flow condition for the test
setGasModel('ideal-air.lua')
inflow = FlowState:new{p=1.0e5, T=300, velx=100.0}

-- Read in GridPro grids.
-- The GridPro grid file contains all the blocks in one file.
-- We'll separate these into individual pieces, and store
-- a list of grids.
scale = 1.0
grids = importGridproGrid(gproGrid, scale)

-- Now we attach the grids to blocks
blks = {}
for i,g in ipairs(grids) do
   blks[#blks+1] = SBlock:new{grid=g, fillCondition=inflow}
end

print("Number of blocks constructed: ", #blks)

-- [TODO]: Apply BCs function not yet available.



