--
-- Unsteady flow over a 10º/70º double cone, details are taken from ref. [1],
-- adapted from a script provided by H. G. Hornung 28-Jan-2020.
--
-- references:
--
-- [1] Unsteadiness Boundaries in Supersonic Flow over Double Cones
--     H.G. Hornung, R.J. Gollan, P.A. Jacobs,
--     Journal of Fluid Mechanics, 2021
--
-- author: Kyle A. Damm
-- date:   2025-09-12
--

-- ==========================================================
-- General settings
-- ==========================================================
config.solver_mode          = "dual_time_stepping"
config.field_format         = "rawbinary"
config.grid_format          = "rawbinary"

-- ==========================================================
-- Freestream conditions
-- ==========================================================
nsp, nmodes, gm = setGasModel('ideal-N2.gas')
gs = GasState:new{gm}

-- define flow state from ref. [1]
T_wall     = 300.0     -- degrees K
gs.T       = 300.0     -- degrees K
gs.p       = 100.0     -- Pa
V_inf      = 2720.0    -- m/s

-- set inflow and initial conditions
inflow = FlowState:new{p=gs.p, T=gs.T, velx=V_inf}
initial = FlowState:new{p=gs.p/5, T=gs.T, velx=0.0}

-- ==========================================================
-- Define the flow domain using an imported grid
-- ==========================================================
grid_location       = "./su2grid/"
config.dimensions   = 2
config.axisymmetric = true
nblocks             = 8

-- import grids
grids = {}
for i=0,nblocks-1 do
   fileName = string.format(grid_location.."block_%d_grid.su2", i)
   grids[i] = registerFluidGrid{
      grid=UnstructuredGrid:new{filename=fileName, fmt="su2text", scale=1},
      fsTag="initial"
   }
end

-- define boundary conditions
bcDict = {
   inflow         = InFlowBC_Supersonic:new{flowCondition=inflow},
   outflow        = OutFlowBC_SimpleExtrapolate:new{},
   wall           = WallBC_NoSlip_FixedT0:new{Twall=T_wall, group="wall"},
   slip           = WallBC_WithSlip:new{},
   METIS_INTERIOR = ExchangeBC_MappedCell:new{cell_mapping_from_file=true, filename=grid_location.."mapped_cells"}
}

-- make fluidblocks
flowDict            = {}
flowDict["initial"] = initial
flowDict["inflow"]  = inflow
makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver settings
-- ==========================================================

-- loads settings
config.write_loads               = true
config.boundary_groups_for_loads = "wall"

-- invsicid flux settings
config.flux_calculator             = "ldfss2"
config.interpolation_order         = 2
config.extrema_clipping            = false
config.thermo_interpolator         = "rhop"
config.unstructured_limiter        = "hvenkat"
config.smooth_limiter_coeff        = 1.0
config.inviscid_least_squares_type = "unweighted_qr"

-- viscous flux settings
config.viscous                    = true
config.spatial_deriv_calc         = "least_squares"
config.spatial_deriv_locn         = "cells"
config.viscous_least_squares_type = "weighted_qr"

-- Set temporal integration settings
config.print_count                    = 1
config.dualtimestepping_update_scheme = "bdf2"
config.fixed_time_step                = false
config.dt_init                        = 1.0e-08
config.dt_max                         = 1.0e-07
config.max_time                       = 1.0e-03
config.max_step                       = 1.0e+06
config.dt_plot                        = 1.0e-05

NewtonKrylovGlobalConfig {
   -- phases
   number_of_phases = 1,

   -- preconditioner settings
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,

   -- nonlinear system solver settings
   max_newton_steps = 100,
   max_consecutive_bad_steps = 10,
   stop_on_relative_residual = 1.0e-02,

   -- linear system solver settings
   frechet_derivative_perturbation = 1.0e-50,
   use_scaling = true,
   max_linear_solver_iterations = 100,
   max_linear_solver_restarts = 0,
   
   -- continuation settings
   inviscid_cfl_only = true,
   use_line_search = false,
   line_search_order = 3,
   use_physicality_check = true,
   allowable_relative_mass_change = 0.9,
   min_relaxation_factor_for_update = 0.1,
   min_relaxation_factor_for_cfl_growth = 0.5,

   -- output settings
   number_of_steps_for_setting_reference_residuals = 0,
   steps_between_status = 1,
   total_snapshots = 5,
   steps_between_snapshots = 1000,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   -- preconditioner settings
   frozen_preconditioner = true,
   use_adaptive_preconditioner = true,
   steps_between_preconditioner_update = 5,

   -- linear system solver settings
   linear_solve_tolerance = 0.1,
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 1,

   -- cfl settings
   use_auto_cfl = true,
   use_local_timestep = true,
   threshold_relative_residual_for_cfl_growth = 0.99,
   start_cfl = 10.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}
