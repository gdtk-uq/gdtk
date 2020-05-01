-- ramp.lua
-- KD, 29-Nov-2017
-- Supersonic Ramp from
--
-- Marques et al.
-- Comparison of Matrix-free Acceleration Techniques in Compressible Navier-Stokes Calculations
-- International Journal for Numerical Methods in Engineering
-- 2004, DOI: 10.1002/nme.1076
--
config.title = "Supersonic Ramp"
print(config.title)
config.dimensions = 2

-- set gas model
setGasModel('ideal-air-gas-model.lua')

-- inflow conditions
g = 1.4   -- gamma
R = 287.0 -- J/kg.K
Mach = 2.0
T_inf = 300.0   -- K
p_inf = 10000.0 -- Pa
a = math.sqrt(g*R*T_inf) -- speed of sound, m/s
u_inf = Mach*a -- m/s

inflow  = FlowState:new{p=p_inf, velx=u_inf, vely=0.0, T=T_inf}
initial  = inflow

nblocks = 4
grids = {}
for i=0,nblocks-1 do
   --fileName = string.format("su2-grid/ramp_grid.su2", i)
   fileName = string.format("su2-grid/block_%d_ramp_grid.su2", i) -- uncomment for multi-block simulations
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

-- set boundary conditions (METIS_INTERIOR unused for single block simulations)
my_bcDict = {INFLOW         = InFlowBC_Supersonic:new{flowCondition=inflow},
	     OUTFLOW        = OutFlowBC_Simple:new{},
	     PARALLEL       = OutFlowBC_Simple:new{},
	     WALL           = WallBC_WithSlip0:new{},
	     METIS_INTERIOR = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, list_mapped_cells=false}}

-- Define fluidBlocks
blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], fillCondition=initial, bcDict=my_bcDict}
end

-- output settings
config.print_count = 40

-- invsicid flux calc settings
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.extrema_clipping = false
config.apply_limiter = true
config.unstructured_limiter = "venkat"
config.venkat_K_value = 0.5
config.freeze_limiter_on_step = 700
config.thermo_interpolator = "rhop"

SteadyStateSolver{
   -- preconditioner
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 50,

   -- scaling
   use_scaling = true,
   use_complex_matvec_eval = true,

   -- start-up and convergence settings
   number_pre_steps = 10,
   number_total_steps = 2000,
   stop_on_relative_global_residual = 1.0e-10,

   -- GMRES settings
   max_outer_iterations = 5,
   max_restarts = 5,

   -- Settings for start-up phase
   number_start_up_steps = 150,
   cfl0 = 1.0,
   eta0 = 0.1,
   tau0 = 1.0,
   sigma0 = 1.0e-30,

   -- Settings for inexact Newton phase
   cfl1 = 1.0,
   tau1 = 1.0,
   sigma1 = 1.0e-30,
   eta1 = 0.1,
   eta_strategy = "constant",

   -- Settings control write-out
   snapshots_count = 10,
   number_total_snapshots = 5,
   write_diagnostics_count = 1
}
