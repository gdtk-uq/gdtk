-- mach4nozzle.lua
-- Unstructured Grid Example -- for use with Eilmer4
-- PJ & RG
-- 2015-11-08 -- adapted from cone20-simple

job_title = "Mach 4 Nozzle Flow"
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title
config.axisymmetric = true
config.viscous = true
config.turbulence_model = "k_omega"
config.flux_calculator = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel("ideal-air-gas-model.lua")
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)

stagCond = FlowState:new{p=1.75e6, T=850.0, velx=0.0, tke=500.0, omega=1.0e6}
initial = FlowState:new{p=1.0e6, T=300.0, velx=300.0, tke=500.0, omega=1.0e6}

-- Define the flow domain using an imported grid
nblocks = 100
grids = {}
for i=0,nblocks-1 do
   fileName = string.format("GRID/block_%d_mach4nozzleStage1.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

my_bcDict = {Inflow=InFlowBC_FromStagnation:new{stagnationState=stagCond},
	          Outflow=OutFlowBC_FixedP:new{p_outside=500.0},
             Axis=WallBC_WithSlip:new{},
             WallFixedT=WallBC_NoSlip_FixedT:new{Twall=300.0},
	          METIS_INTERIOR=ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], initialState=initial, bcDict=my_bcDict}
end

-- Do a little more setting of global data.
config.max_time = 1.0e-3  -- seconds
config.max_step = 30000000
config.dt_init = 1.0e-10
config.dt_plot = config.max_time/10
config.dt_history = 10.0e-5
config.extrema_clipping = false

config.viscous_delay = 5.0e-3

config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"

config.cfl_schedule_times = {0.0, 1.0e-3, 2.0e-3,}
config.cfl_schedule_values = {0.5, 0.5, 0.8}

config.adjust_invalid_cell_data = true
config.max_invalid_cells = 20

-- Distribute the blocks
mpiTasks = mpiDistributeBlocks{ntasks=6, dist="load-balance"}