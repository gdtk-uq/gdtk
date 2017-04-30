-- ramp.lua
-- A simple 3D simulation of flow over a ramp with 10-degree deflection.
-- PJ & RG
-- 2018-08-17 -- with su2 grid generated Kyle Damm
config.title = "Mach 1.5 flow over a 10-degree ramp."
print(config.title)
config.dimensions = 3
config.viscous = false
config.max_time = 5.0e-3
config.max_step = 3000
config.gasdynamic_update_scheme = "euler"
config.unstructured_limiter = "venkat"
config.interpolation_order = 1
config.dt_plot = 1.0e-3
config.dt_history = 1.0e-5

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0, vely=0.0}

-- Define the flow domain using an imported grid.
grid0 = UnstructuredGrid:new{filename="grid0.su2", fmt="su2text"}
grid1 = UnstructuredGrid:new{filename="grid1.su2", fmt="su2text"}
my_bcDict = {INFLOW=InFlowBC_Supersonic:new{flowCondition=inflow, label="inflow-boundary"},
	     OUTFLOW=OutFlowBC_Simple:new{label="outflow-boundary"},
	     RAMP_SURFACE=WallBC_WithSlip:new{},
	     INTERIOR=ExchangeBC_MappedCell:new{list_mapped_cells=true}}
blk0 = FluidBlock:new{grid=grid0, fillCondition=inflow, bcDict=my_bcDict}
blk1 = FluidBlock:new{grid=grid1, fillCondition=initial, bcDict=my_bcDict}
