-- icon.lua
-- A fluidic oscillator to produce an unsteady jet. (long-body version)
-- PJ 2019-07-30, 2019-08-01

job_title = "Fluidic oscillator -- long-body."
print(job_title)
config.dimensions = 2
config.title = job_title

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=3.0e3, T=304.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=10.0e3, T=303.0, velx=400.0, vely=0.0}

print("T=", inflow.T, "density=", inflow.rho, "sound speed= ", inflow.a)
print("inflow Mach number=", 1000.0/inflow.a)
print("dynamic pressure q=", 1/2*inflow.rho*1.0e6)

-- Define the flow domain using an imported grid
nblocks = 4
grids = {}
for i=0,nblocks-1 do
   fileName = string.format("block_%d_fjl-gmsh.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

my_bcDict = {INFLOW=InFlowBC_Supersonic:new{flowState=inflow, label="inflow-boundary"},
	     OUTFLOW=OutFlowBC_Simple:new{label="outflow-boundary"},
	     SLIP_WALL=WallBC_WithSlip:new{},
	     CAVITY=WallBC_WithSlip:new{},
	     INJECTOR=WallBC_WithSlip:new{},
	     METIS_INTERIOR=ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], initialState=initial, bcDict=my_bcDict}
end

config.dt_init = 1.0e-9
config.max_time = 5.0e-3  -- seconds
config.max_step = 400000
config.dt_plot = 10.0e-6
