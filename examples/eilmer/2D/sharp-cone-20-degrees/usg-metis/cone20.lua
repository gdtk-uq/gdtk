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
Q = inflow:toTable()
print("T=", Q.T, "density=", Q.rho, "sound speed= ", Q.a)
print("inflow Mach number=", 1000.0/Q.a)
print("dynamic pressure q=", 1/2*Q.rho*1.0e6)

-- Define the flow domain using an imported grid
nblocks = 4
grids = {}
for i=0,nblocks-1 do
   fileName = string.format("block_%d_cone20.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

my_bcDict = {INFLOW=InFlowBC_Supersonic:new{flowState=inflow, label="inflow-boundary"},
	     OUTFLOW=OutFlowBC_Simple:new{label="outflow-boundary"},
	     SLIP_WALL=WallBC_WithSlip:new{},
	     METIS_INTERIOR=ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], initialState=inflow, bcDict=my_bcDict}
end

-- There are 8 boundaries within the single unstructured-grid
-- These are essentially the same as the two-block arrangement of the other examples,
-- so we need to join the interior boundaries.
-- Unspecified boundaries will be SlipWallBC.

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5

setHistoryPoint{x=1.0, y=0.2} -- nose of cone
setHistoryPoint{x=0.201, y=0.001} -- base of cone
