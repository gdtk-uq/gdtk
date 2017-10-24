-- cone20.lua
-- Simple job-specification file for e4prep -- for use with Eilmer4
-- PJ & RG
-- 2015-02-24 -- adapted from the Python version of cone20
-- 2017-05-07 -- try import of Eilmer3 grids

-- We can set individual attributes of the global data object.
config.title = "Mach 1.5 flow over a 20 degree cone."
print(config.title)
config.dimensions = 2
config.axisymmetric = true

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}

grids = {}
for i=0, 1 do
   fileName = string.format("eilmer3-grids/cone20.grid.b%04d.t0000", i)
   grids[i] = StructuredGrid:new{file=fileName, fmt="text"}
end

-- Define the flow-solution blocks.
blk0 = FluidBlock:new{grid=grids[0], initialState=inflow}
blk1 = FluidBlock:new{grid=grids[1], initialState=initial}
-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_Supersonic:new{flowState=inflow}
blk1.bcList[east] = OutFlowBC_Simple:new{}

-- add history point 2/3 along length of cone surface
nx1 = grids[1]:get_niv()
setHistoryPoint{ib=1, i=math.floor(2*nx1/3), j=0}

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5
