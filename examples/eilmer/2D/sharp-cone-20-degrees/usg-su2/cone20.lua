-- cone20.lua
-- Unstructured Grid Example -- for use with Eilmer4
-- PJ & RG
-- 2015-11-08 -- adapted from cone20-simple

job_title = "Mach 1.5 flow over a 20 degree cone -- Unstructured Grid."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title
config.axisymmetric = true

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0, vely=0.0}

-- Demo: Verify Mach number of inflow.
Q = GasState:new{gm}
Q.p = 95.84e3
Q.T = {1103.0}
print("T", Q.T[1])
Q.massf = {air=1.0}
gm:updateSoundSpeed(Q)
print("Sound speed= ", Q.a)
print("Inflow Mach number= ", 1000.0/Q.a)

-- Define the flow domain using an imported grid.
mygrid = UnstructuredGrid:new{filename="cone20.su2", fmt="su2text"}
myblk = UBlock:new{grid=mygrid, fillCondition=initial}
-- There are 8 boundaries within the single unstructured-grid
-- These are essentially the same as the two-block arrangement of the other examples,
-- so we need to join the interior boundaries.
-- Unspecified boundaries will be SlipWallBC.
myblk.bcList[6] = InFlowBC_Supersonic:new{flowCondition=inflow, label="inflow-boundary"} -- WEST0
myblk.bcList[7] = ExchangeBC_MappedCell:new{} -- WEST1
myblk.bcList[0] = ExchangeBC_MappedCell:new{} -- EAST0
myblk.bcList[1] = OutFlowBC_Simple:new{label="outflow-boundary"} -- EAST1

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5

setHistoryPoint{x=1.0, y=0.2} -- nose of cone
setHistoryPoint{x=0.201, y=0.001} -- base of cone
