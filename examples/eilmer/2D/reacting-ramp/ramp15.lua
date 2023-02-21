-- ramp15.lua
-- Reacting flow over a 15 degree ramp
-- Kyle Damm 
-- 2020-11-09
-- Updated by Nick Gibbons to test svan_albada (2023-02-11)

job_title = "Mach 3.4 reacting flow over a 15 degree ramp."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title
config.axisymmetric = false

config.unstructured_limiter = "svan_albada"
config.extrema_clipping=false                     -- On for safety, off for deep convergence
config.scale_species_after_reconstruction = false -- On for safety, off for deep convergence
config.enforce_species_density_positivity = false -- On for safety, off for deep convergence
config.smooth_limiter_coeff = 1e-3

config.flux_calculator = "ausmdv"
config.interpolation_order = 2

config.reacting = true
config.reactions_file = 'reac-file.lua'

nsp, nmodes, gmodel = setGasModel('gas-model.lua')
Mach = 3.4
-- inflow gas state
gs = GasState:new{gmodel}
gs.T = 700.0 -- Kelvin
gs.p = 101e3 -- Pa, 1 atm
gs.massf = {H2=0.005, O2=0.24, N2=0.755} -- general approximation of air
gamma = gmodel:gamma(gs)
R = gmodel:R(gs)
vel_inf = Mach*math.sqrt(gs.T*gamma*R)
inflow = FlowState:new{p=gs.p, T=gs.T, velx=vel_inf, massf=gs.massf}
initial = inflow -- FlowSolution:new{jobName="ramp15", dir="../explicit/", tindx=299, nBlocks=8}

-- Define the flow domain using an imported grid.
nblocks=4
grids = {}
for i=0,nblocks-1 do
   fileName = string.format("su2-grid/block_%d_ramp15.su2", i)
   --fileName = string.format("su2-grid/ramp15.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1}
end

my_bcDict = {north = OutFlowBC_Simple:new{},
             east  = OutFlowBC_Simple:new{},
             south = WallBC_WithSlip:new{},
             west  = InFlowBC_Supersonic:new{flowCondition=inflow},
             METIS_INTERIOR   = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, list_mapped_cells=false}}

blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], fillCondition=initial, bcDict=my_bcDict}
end


-- Solver settings
config.print_count = 20



SteadyStateSolver{
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   frozen_preconditioner_count = 40,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   use_physicality_check = true,
   
   number_total_steps = 2000,
   stop_on_relative_global_residual = 1.0e-8,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 40,
   max_restarts = 0,
   residual_based_cfl_scheduling = true,
   cfl_max = 1e6,

   -- Settings for start-up phase
   number_start_up_steps = 0,
   cfl0 = 1e0,
   tau0 = 0.5,
   eta0 = 0.01,
   sigma0 = 1.0e-30,
   p0 = 0.8,

   -- Settings for inexact Newton phase
   cfl1 = 2e-1,
   eta1 = 0.01,
   tau1 = 0.1,
   p1 = 0.8,
   sigma1 = 1.0e-30,
   eta_strategy = "constant",

   -- Settings control write-out
   snapshots_count = 20,
   number_total_snapshots = 5,
   write_diagnostics_count = 1,
   write_loads_count = 1000,
}
