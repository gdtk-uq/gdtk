--
-- 7-species air thermochemical nonequilibrium flow over a sphere.
--
-- Details (including sphere3.g grid file) are taken from the LARC hypersonic benchmark test cases website:
-- https://fun3d.larc.nasa.gov/chapter-2.html#hypersonic_benchmarks
--
-- @author: Kyle A. Damm (2024-08-28)
--

-- ==========================================================
-- General settings
-- ==========================================================
config.solver_mode  = "steady"
config.field_format = "rawbinary"
config.grid_format  = "rawbinary"

-- ==========================================================
-- Freestream conditions
-- ==========================================================
nsp, nmodes, gm = setGasModel('air-7sp-2T.gas')
gs = GasState:new{gm}

-- free stream is for flow at Mach 17.6 and Reynolds number of 376930.0 /m
gs.massf   = {N2=0.767,O2=0.233}
gs.rho     = 0.001  -- kg/m^3
gs.T       = 200.0  -- K
gs.T_modes = {gs.T} -- assume freestream is in thermodynamic equilibrium
V_inf      = 5000.0 -- m/s
gm:updateThermoFromRHOT(gs)

-- set inflow and initial condition
flowDict            = {}
inflow              = FlowState:new{p=gs.p, T=gs.T, T_modes=gs.T_modes, velx=V_inf, massf=gs.massf}
flowDict["inflow"]  = inflow
initial             = inflow
flowDict["initial"] = initial

-- ==========================================================
-- Define the flow domain using an imported grid
-- ==========================================================
grid_location       = "./su2_grid/"
config.dimensions   = 2
config.axisymmetric = true
nblocks             = 1

-- import grids
grids = {}
for i=0,nblocks-1 do
   fileName = string.format(grid_location.."block_%d_sphere.su2", i)
   grids[i] = registerFluidGrid{
      grid=UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1},
      fsTag="initial"
   }
end

-- define boundary conditions
T_wall = 500.0 -- K
bcDict = {
   inflow         = InFlowBC_Supersonic:new{flowCondition=inflow},
   outflow        = OutFlowBC_SimpleExtrapolate:new{},
   wall           = WallBC_NoSlip_FixedT0:new{catalytic_type="none", Twall=T_wall, group="wall"},
   symmetry       = WallBC_WithSlip:new{},
   METIS_INTERIOR = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, filename=grid_location.."mapped_cells"}
}

-- make fluidblocks
makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver settings
-- ==========================================================

-- reactions settings
config.reacting       = true
config.reactions_file = 'air-7sp-2T.chem'
config.enforce_species_density_positivity = false
config.scale_species_after_reconstruction = false

-- energy modes settings
config.energy_exchange_file = 'air-7sp-2T.exch'

-- diffusion model settings
config.mass_diffusion_model       = "ficks_first_law"
config.diffusion_coefficient_type = "binary_diffusion"

-- invsicid flux settings
config.flux_calculator             = "adaptive_ldfss0_ldfss2"
config.shock_detector              = "KAD"
config.compression_tolerance       = -1.0
config.interpolation_order         = 2
config.inviscid_least_squares_type = "unweighted_qr"
config.extrema_clipping            = false
config.thermo_interpolator         = "rhop"
config.unstructured_limiter        = "hvenkat"
config.smooth_limiter_coeff        = 1.0
config.apply_unstructured_limiter_stagnation_point_filter = true
config.apply_unstructured_limiter_min_pressure_filter     = true

-- viscous flux settings
config.viscous                    = true
config.spatial_deriv_calc         = "least_squares"
config.spatial_deriv_locn         = "cells"
config.viscous_least_squares_type = "weighted_qr"

-- Set temporal integration settings
NewtonKrylovGlobalConfig{
   -- phases
   number_of_phases = 3,
   max_steps_in_initial_phases = { 500, 500 },
   phase_changes_at_relative_residual = { 1.0e-03, 1.0e-06 },

   -- Preconditioner
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,

   -- nonlinear system solver settings
   max_newton_steps = 2000,
   max_consecutive_bad_steps = 10,
   stop_on_relative_residual = 1.0e-10,

   -- linear system solver settings
   frechet_derivative_perturbation = 1.0e-50,
   use_scaling = true,
   max_linear_solver_iterations = 100,
   max_linear_solver_restarts = 0,
   
   -- continuation settings
   inviscid_cfl_only = true,
   use_line_search = true,
   line_search_order = 3,
   use_physicality_check = true,
   allowable_relative_mass_change = 0.5,
   min_relaxation_factor_for_update = 0.1,
   min_relaxation_factor_for_cfl_growth = 0.5,
   use_residual_smoothing = true,

   -- output settings
   number_of_steps_for_setting_reference_residuals = 0,
   steps_between_status = 1,
   write_limiter_values = true,
   total_snapshots = 5,
   steps_between_snapshots = 50,
   steps_between_diagnostics = 1,
}

NewtonKrylovPhase:new{
   -- preconditioner settings
   frozen_preconditioner = true,
   use_adaptive_preconditioner = true,
   steps_between_preconditioner_update = 5,

   -- linear system solver settings
   linear_solve_tolerance = 0.01,
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 1,

   -- cfl settings
   use_auto_cfl = true,
   use_local_timestep = true,
   threshold_relative_residual_for_cfl_growth = 0.99,
   start_cfl = 1.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0,
   limit_on_cfl_decrease_ratio = 1.0
}

NewtonKrylovPhase:new{
   -- linear system solver settings
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = -1.0
}

NewtonKrylovPhase:new{
   frozen_shock_detector = true,
   frozen_limiter_for_residual = true,
   frozen_limiter_for_jacobian = true,
   start_cfl = -1
}
