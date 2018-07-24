-- ramp.lua
-- KD, 24-July-2018
-- Supersonic Ramp from Marques and Pereira paper 

config.title = "2D Supersonic Ramp, air at Mach 2.0"
print(config.title)
config.dimensions = 2

-- set gas model
setGasModel('ideal-air-gas-model.lua')

-- Define flow conditions
g = 1.4 -- gamma
R = 287.0 -- J/kg.K
Mach = 2.0
T_inf = 300.0 -- K
p_inf = 10000.0 -- Pa
a = math.sqrt(g*R*T_inf)
u_inf = Mach*a
inflow  = FlowState:new{p=p_inf, velx=u_inf, vely=0.0, T=T_inf}

-- Define grid
fileName = string.format("supersonic_ramp.su2")
grid = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}

-- set boundary conditions (METIS_INTERIOR unused for single block simulations)
my_bcDict = {INFLOW=InFlowBC_Supersonic:new{flowCondition=inflow},
	     OUTFLOW=OutFlowBC_Simple:new{},
	     PARALLEL=OutFlowBC_Simple:new{},
	     WALL=WallBC_WithSlip:new{}}

-- Define fluidBlocks
blk = UBlock:new{grid=grid, fillCondition=inflow, bcDict=my_bcDict}

-- Set simulation parameters
config.unstructured_limiter = "van_albada"
config.flux_calculator = "roe"
config.print_count = 40
-- Set reconstruction parameters
config.interpolation_order = 2
config.extrema_clipping = false
config.use_extended_stencil = false
config.apply_limiter = true
config.suppress_reconstruction_at_boundaries = true
config.thermo_interpolator = "rhop"
-- Set steady-state solver parameters
SteadyStateSolver{
   use_preconditioner = true,
   use_scaling = false,
   use_complex_matvec_eval = true,
   -- sigma = 1.0e-6, -- presently it's computed internally
   number_pre_steps = 0,
   number_total_steps = 2000,
   stop_on_relative_global_residual = 1.0e-12,
   -- Settings for FGMRES iterative solver
   max_outer_iterations = 30,
   --   number_inner_iterations = 10, -- not needed is preconditioning is false
   max_restarts = 10,
   -- Settings for start-up phase
   number_start_up_steps = 10,
   cfl0 = 1.0,
   eta0 = 0.5,
   tau0 = 0.8,
   sigma0 = 1.0e-50,
   -- Settings for inexact Newton phase
   cfl1 = 10.0,
   tau1 = 0.8,
   sigma1 = 1.0e-50,
   eta1 = 0.01,
   eta1_min = 0.01,
   eta_strategy = "constant",
   -- Settings control write-out
   snapshots_count = 10,
   number_total_snapshots = 5,
   write_diagnostics_count = 1
}
