-- Simulation of RAMCII: Coarse viscous simulation
--
-- @author: Robert Watt
-- 2021-05-01

job_title = "RAMCII coarse viscous"
config.solver_mode = "steady"
config.dimensions = 2
config.axisymmetric = true 

config.max_step = 200

config.flux_calculator = "hanel"

nsp, nmodes, gmodel = setGasModel('gas.gas')
config.reacting = true
config.reactions_file = 'chem.chem'
config.energy_exchange_file = 'kinetics.kin'

dofile('conditions.lua')

-- Compute full inflow state from T,rho, and u
Q = GasState:new{gmodel}
Q.T = T_inf;
Q.rho = rho_inf;
Q.massf = mass_fraction
Q.T_modes = {T_inf, T_inf}
gmodel:updateThermoFromRHOT(Q);
p_inf = Q.p

inflow = FlowState:new{
   p=p_inf, 
   T=T_inf, 
   velx=u_inf,
   vely=0.0, 
   massf={N2=0.767, O2=0.233},
   T_modes={T_inf, T_inf}
}
initial = FlowState:new{
   p=p_inf, 
   T=T_inf, 
   velx=u_inf,
   vely=0.0, 
   massf={N2=0.767, O2=0.233},
   T_modes={T_inf, T_inf}
}

bcDict = {
   outflow = OutFlowBC_Simple:new{},
   inflow = InFlowBC_Supersonic:new{flowState=inflow},
   wall = WallBC_WithSlip:new{group="wall"},
   symm = WallBC_WithSlip0:new{},
}

fsDict = {
   initial = initial
}

makeFluidBlocks(bcDict, fsDict)
mpiDistributeBlocks{ntasks=8}

config.thermo_interpolator = "rhop"
config.apply_heuristic_pressure_based_limiting = true
config.epsilon_van_albada = 1e-300
config.scale_species_after_reconstruction = false
config.extrema_clipping = false

config.T_frozen = 1
config.T_frozen_energy = 1

config.boundary_groups_for_loads = "wall"


config.electron_pressure_convection_term = true
config.electric_field_work = true

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 20000,
   stop_on_relative_residual = 1e-8,
   number_of_phases = 2,
   write_residual_values = true,
   max_steps_in_initial_phases = {10000},
   phase_changes_at_relative_residual = {
      1e-3,
   },
   inviscid_cfl_only = false,
   use_line_search = false,
   use_physicality_check = true,
   max_linear_solver_iterations = 150,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1e-150,
   preconditioner_perturbation = 1e-150,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 10,
   steps_between_snapshots = 10,
   steps_between_diagnostics = 1,

   write_loads = true,
   steps_between_loads_update = 10000,
   write_loads_on_last_step = true,
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 10,
   linear_solve_tolerance = 1e-2,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.5,
   start_cfl = 0.25,
   max_cfl = 1e6,
   auto_cfl_exponent = 0.8,

   -- use_local_timestep = false,
   -- grid_motion_enabled = false,
   -- shock_fitting_scale_factor = 1.0,
}

-- NewtonKrylovPhase:new{
--    residual_interpolation_order = 1,
--    jacobian_interpolation_order = 1,
--    use_adaptive_preconditioner = false,
--    steps_between_preconditioner_update = 5,
--    use_auto_cfl = true,
--    threshold_relative_residual_for_cfl_growth = 0.1,
--    start_cfl = 50.0,
--    linear_solve_tolerance = 1e-2,
--    frozen_limiter_for_jacobian = false,
--    frozen_limiter_for_residual = false,
--    auto_cfl_exponent = 0.6,
--    -- grid_motion_enabled = true,
-- }


NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 5,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.01,
   start_cfl = 5.0,
   linear_solve_tolerance = 1e-2,
   frozen_limiter_for_jacobian = false,
   frozen_limiter_for_residual = false,
   auto_cfl_exponent = 0.6,
}

