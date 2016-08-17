-- ramp.lua
-- A simple 3D simulation of flow over a ramp with 10-degree deflection.
-- PJ & RG
-- 2018-08-17 -- with su2 grid generated Kyle Damm
config.title = "Mach 1.5 flow over a 10-degree ramp."
print(config.title)
config.dimensions = 3
config.viscous = false
config.max_time = 5.0e-3
config.max_step = 50
config.gasdynamic_update_scheme = "euler"
config.interpolation_order = 2
config.dt_plot = 1.0e-3
config.dt_history = 1.0e-5

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0, vely=0.0}

-- Define the flow domain using an imported grid.
mygrid = UnstructuredGrid:new{filename="simpleRamp.su2", fmt="su2text"}
myblk = UBlock:new{grid=mygrid, fillCondition=initial}
-- There are 12 boundaries within the single unstructured-grid
-- as an artifact of the original two blocks.
-- These are essentially the same as the two-block arrangement of the other examples,
-- so we need to join the interior boundaries.
-- Unspecified boundaries will be SlipWallBC.
myblk.bcList[10] = InFlowBC_Supersonic:new{flowCondition=inflow, label="inflow-boundary"} -- WEST0
myblk.bcList[11] = ExchangeBC_MappedCell:new{} -- WEST1
myblk.bcList[2] = ExchangeBC_MappedCell:new{} -- EAST0
myblk.bcList[3] = OutFlowBC_Simple:new{label="outflow-boundary"} -- EAST1

