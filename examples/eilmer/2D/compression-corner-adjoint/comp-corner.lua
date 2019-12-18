-- Author: KAD
-- Date: 2019-12-18
-- Mach 6 flow of ideal air over a 7.5-deg ramp.

config.title = string.format("Mach 6.0 air flowing over a 7.5 degree ramp.")
print(config.title)

-- domain parameters
config.axisymmetric = false
config.dimensions = 2

-- turbulence model
config.turbulence_model = "k_omega"

-- set gas model
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')

-- set inflow conditions
M_inf = 6.0
ReL = 1.0e07     -- m**(-1)
Tt_inf = 470.0   -- K
gamma = 1.4      -- for ideal air
R = 287.0        -- J/kg.K for ideal air

T_inf = Tt_inf * (1.0 / ((1.0 + M_inf*M_inf*(gamma-1.0)/2.0))) -- K
u_inf = M_inf * math.sqrt(gamma*R*T_inf)  -- m/s

gas_inf = FlowState:new{p=p_inf, T=T_inf}

mu_inf = gas_inf.mu  -- kg/m.s
p_inf = (ReL * R * T_inf * mu_inf)/u_inf -- Pa 

-- Turbulence quantities estimate
I_inf = 0.03  -- 3% turbulence intensity
mu_t_mu_inf = 10.0

tke_inf = 1.5 * (u_inf * I_inf)^2
mu_t_inf = mu_t_mu_inf * gas_inf.mu
omega_inf = gas_inf.rho * tke_inf / mu_t_inf

-- Set up inflow condition
inflow = FlowState:new{p=p_inf, T=T_inf, velx = u_inf, vely=0.0, tke=tke_inf, omega=omega_inf}
initial = inflow

-- Define the flow domain using an imported grid.
nblocks=4
grids = {}
for i=0,nblocks-1 do
   fileName = string.format("su2-grid-files/block_%d_ramp.su2", i)
   grids[i] = UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1.0}
end

-- set boundary conditions
my_bcDict = {INFLOW         = InFlowBC_Supersonic:new{flowCondition=inflow, group="inflow"},
	     OUTFLOW        = OutFlowBC_Simple:new{group="outflow"},
	     TOP            = WallBC_WithSlip0:new{group="top"},
	     RAMP           = WallBC_NoSlip_FixedT0:new{Twall=300.0, is_design_surface=true, num_cntrl_pts=4, group="objective_function_surface"},
	     WALL           = WallBC_NoSlip_FixedT0:new{Twall=300.0, group="wall"},
             METIS_INTERIOR = ExchangeBC_MappedCell:new{cell_mapping_from_file=true}}


-- create fluidblocks
blks = {}
for i=0,nblocks-1 do
   blks[i] = FluidBlock:new{grid=grids[i], initialState=initial, bcDict=my_bcDict}
end

-- general solver settings
config.print_count = 50

-- convective flux settings
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.apply_limiter = true
config.extrema_clipping = false
config.unstructured_limiter = "venkat"
config.venkat_K_value = 0.3
config.freeze_limiter_on_step = 1000
config.thermo_interpolator = "rhop"

-- viscous settings
config.viscous = true
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
config.diffuse_wall_bcs_on_init = true
config.number_init_passes = 10

-- flow solver settings
SteadyStateSolver{
   
   use_preconditioner = true,
   precondition_matrix_type = "ilu",
   ilu_fill = 0,
   frozen_preconditioner_count = 10,
   start_preconditioning = 1,
   
   use_scaling = true,
   use_complex_matvec_eval = true,
   
   number_pre_steps = 10,
   number_total_steps = 1200,
   stop_on_relative_global_residual = 1.0e-20,

   -- Settings for FGMRES iterative solver
   max_outer_iterations = 10,
   max_restarts = 10,

   -- Settings for start-up phase
   number_start_up_steps = 500,
   sigma0 = 1.0e-30,
   cfl0 = 0.01,
   p0 = 0.75,
   tau0 = 1.0,
   eta0 = 0.1,

   -- Settings for inexact Newton phase
   sigma1 = 1.0e-30,
   cfl1 = 1.0,
   p1 = 1.0,
   tau1 = 1.0,
   eta1 = 0.1,
   eta1_min = 1.0e-01,
   eta_ratio_per_step = 0.99,
   eta_strategy = "geometric",

   -- Settings control write-out
   snapshots_count = 100,
   number_total_snapshots = 5,
   write_diagnostics_count = 50
}

-- adjoint solver settings
ShapeSensitivityCalculator{
   read_frozen_limiter_values_from_file = true,
   pseudotime = false,
   cfl0 = 1.0e03,
   epsilon = 1.0e-30,
   eta = 1.0e-02,
   maxOuterIterations = 20,
   maxRestarts = 1000,
   stop_on_relative_global_residual = 1.0e-25,
   tol_bezier_curve_fit = 1.0e-06,
   max_steps_bezier_curve_fit = 10000
}
