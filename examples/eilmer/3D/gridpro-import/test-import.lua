-- Author: Rowan J. Gollan
-- Date: 2016-04-22
--
-- History:
--   2018-02-14 : Use Fabian Zander's example with BC labels.

gproGrid = "blk.tmp"
gproConn = "blk.tmp.conn"
gproPty = "blk.tmp.pty"

config.title = "Gridpro importer test"
config.dimensions = 3

config.max_time = 2.0e-4
config.max_step = 10000000
config.dt_init = 1.0e-9
config.dt_plot = config.max_time/10.0

-- Set up a dummy flow condition for the test
setGasModel('ideal-air.lua')
initial = FlowState:new{p=1.0e4, T=300, velx=0.0}
inflow = FlowState:new{p=1.0e5, T=300, velx=-1000.0}

-- Read in GridPro grids.
-- The GridPro grid file contains all the blocks in one file.
-- We'll separate these into individual pieces, and store
-- a list of grids.
scale = 1.0
grids = importGridproGrid(gproGrid, scale)

-- Now we attach the grids to blocks
blks = {}
for i,g in ipairs(grids) do
   blks[#blks+1] = FluidBlock:new{grid=g, initialState=initial}
end

print("Number of blocks constructed: ", #blks)

applyGridproConnectivity(gproConn, blks)

bcMap = {SLIP_WALL=WallBC_WithSlip:new{},
   	 FIXED_T=WallBC_NoSlip_FixedT:new{Twall=30000},
   	 INFLOW_SUPERSONIC=InFlowBC_Supersonic:new{flowState=inflow},
   	 OUTFLOW_SIMPLE=OutFlowBC_Simple:new{}}

applyGridproBoundaryConditions(gproPty, blks, bcMap)
