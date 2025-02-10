-- Simulation of the phoebus capsule geometry from ref. [1].
--
-- References:
-- [1] D. Bianchi et al. “Numerical Analysis and Wind Tunnel Validation of 
--     Low-Temperature Ablators undergoing Shape Change”. In: International 
--     Journal of Heat and Mass Transfer 177 (2021), p. 121430. issn: 0017-9310.
--     doi: 10.1016/j.ijheatmasstransfer.2021.121430.
--
-- [2] A. L. Zibitsker et al. “Validation and analysis of a coupled 
--     fluid-ablation framework for modeling low-temperature ablator”. 
--     In: International Journal of Heat and Mass Transfer 218 (2024), 
--     p. 124728. issn: 0017-9310. doi: 10.1016/j.ijheatmasstransfer.2023.124728.
--
-- @author: Reece B. Otto (2024-11-19)

config.dimensions = 2
config.axisymmetric = true
config.solver_mode = "steady"

-- gas model
nsp, nmodes, gmodel = setGasModel("ideal-air.gas")

-- free-stream conditions (test 2 - 20bar from ref. [2])
u_inf = 927.0   -- [m/s]
p_inf = 1302.8  -- [Pa]
T_inf = 59.4    -- [K]
T_wall = 298.15 -- [K]

-- flow states
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0}
initial = inflow
flowDict = {
  initial=initial,
  inflow=inflow
}

-- boundary conditions
bcDict = {
  inflow = InFlowBC_ShockFitting:new{flowState=inflow},
  wall = WallBC_NoSlip_FixedT:new{group="wall", Twall=T_wall},
  outflow = OutFlowBC_Simple:new{},
  symmetry = WallBC_WithSlip:new{}
}

makeFluidBlocks(bcDict, flowDict)

-- ==========================================================
-- Solver configuration
-- ==========================================================
-- record surface loads
config.write_loads = true
config.boundary_groups_for_loads = "wall"

-- invsicid flux settings
config.flux_calculator = "ausmdv"
config.apply_entropy_fix = false
config.interpolation_order = 2
config.extrema_clipping = false
config.thermo_interpolator = "rhop"
config.apply_limiter = false

-- shock-fitting
config.shock_fitting_filter_velocity_scale = 1.0
config.grid_motion = "shock_fitting"

-- viscous flux settings
config.viscous = true

-- settings for steady-state solver
NewtonKrylovGlobalConfig{
  use_preconditioner = false, -- use fgmres instead
  use_fgmres = true, -- fgmres preconditions grid movement terms, ilu does not
  max_fgmres_preconditioning_iterations = 30,
  number_of_steps_for_setting_reference_residuals = 2,
  max_consecutive_bad_steps = 10,
  max_newton_steps = 3000,
  stop_on_relative_residual = 1.0e-9,
  inviscid_cfl_only = true,
  use_line_search = false,
  use_physicality_check = true,
  max_linear_solver_iterations = 400,
  max_linear_solver_restarts = 0,
  use_scaling = true,
  frechet_derivative_perturbation = 1.0e-250,
  preconditioner_perturbation = 1.0e-250,
  total_snapshots = 3,
  steps_between_snapshots = 10,
  steps_between_diagnostics = 1,
  steps_between_status = 10,
  write_snapshot_on_last_step = true,
  write_diagnostics_on_last_step = true,
  write_loads = true,
  write_residual_values = true,
  number_of_phases = 3,
  phase_changes_at_relative_residual = {5.0e-3, 5.0e-5},
  max_steps_in_initial_phases = {1000, 1000}
}

-- The grid is stationary in the first phase to give the
-- shock time to establish itself
NewtonKrylovPhase:new{
  residual_interpolation_order = 1,
  jacobian_interpolation_order = 1,
  linear_solve_tolerance = 1e-2,
  fgmres_preconditioning_solve_tolerance = 1e-2,
  use_auto_cfl = true,
  threshold_relative_residual_for_cfl_growth = 0.99,
  start_cfl = 0.5,
  max_cfl = 1.0e4,
  auto_cfl_exponent = 0.8,
  use_local_timestep = false
}

-- The grid begins moving with first order numerics
-- to fit the shock
NewtonKrylovPhase:new{
  start_cfl = -1,
  linear_solve_tolerance = 1e-2,
  fgmres_preconditioning_solve_tolerance = 1e-2,
  grid_motion_enabled = true,
  auto_cfl_exponent = 0.5
}

-- switch to second order numerics once the shock has been fit
NewtonKrylovPhase:new{
  residual_interpolation_order = 2,
  jacobian_interpolation_order = 2,
  linear_solve_tolerance = 1e-2,
  fgmres_preconditioning_solve_tolerance = 1e-2,
  start_cfl = -1
}
