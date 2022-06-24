-- 3D blunt cone test simulation for parallel scaling tests
-- Notes: This file uses some odd values for the steady-state solver
--        settings to get consistent timing data.
-- @author: Nick Gibbons (02/06/22)

job_title = "Mach 5.0 flow over a 5 degree cone -- 3D."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 3
config.flux_calculator = "hanel"
config.viscous=true
config.unstructured_limiter = "park"
config.thermo_interpolator = "rhop"

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
inflow = FlowState:new{p=10e3, T=300.0, velz=-1735.9, vely=0.0}

-- Read in grid:
for blockfiles in io.popen("ls -1 ./su2grid | wc -l"):lines() do
   nblocks = tonumber(blockfiles)
end

grids = {}
for i=0,nblocks-1 do
   fileName = string.format("su2grid/block_%d_grid.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

my_bcDict = {top=InFlowBC_Supersonic:new{flowState=inflow, label="inflow-boundary"},
             bottom=WallBC_NoSlip_FixedT:new{Twall=300.0},
             east=OutFlowBC_Simple:new{label="outflow-boundary"},
             METIS_INTERIOR=ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], initialState=inflow, bcDict=my_bcDict}
end

-- Set up SSS config data

SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 10000,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   use_physicality_check = true,
   
   number_total_steps = 40,
   stop_on_relative_global_residual = 1.0e-8,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 10,
   max_restarts = 0,
   residual_based_cfl_scheduling = true,
   cfl_max = 5e1,

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 1e0,
   tau0 = 0.1,
   eta0 = 0.1,
   sigma0 = 1.0e-30,
   p0 = 0.8,

   -- Settings for inexact Newton phase
   cfl1 = 1e0,
   tau1 = 0.7,
   eta1 = 1e-6,
   sigma1 = 1.0e-30,
   p1 = 0.7,
   eta_strategy = "constant",

   -- Settings control write-out
   snapshots_count = 100000,
   number_total_snapshots = 100,
   write_diagnostics_count = 1,
   write_loads_count = 1000,
}
