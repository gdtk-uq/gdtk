-- steady.lua
print("Set up steady-state solve for Hakkinen et al's 1959 experiment.")
--
-- 0. Assume that a previous processing step has set up the grids.
--
-- 1. Domain type, gas model and flow states
config.solver_mode = "steady"
config.dimensions = 2
--
-- Flow conditions to match those of Figure 6: pf/p0=1.4, Re_shock=2.96e5
p_inf = 6205.0 -- Pa
u_inf = 514.0 -- m/s
T_inf = 164.4 -- degree K
-- End of plate, and of the whole flow domain.
mm = 1.0e-3 -- metres per mm
L2 = 90.0*mm

nsp, nmodes = setGasModel('ideal-air.gas')
print('GasModel set to ideal air. nsp= ', nsp, ' nmodes= ', nmodes)
inflow = FlowState:new{p=p_inf, velx=u_inf, T=T_inf}
flowDict = {initial=inflow, inflow=inflow}
--
-- 2. boundary conditions attached to fluid blocks.
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_FixedPT:new{p_outside=p_inf, T_outside=T_inf},
   noslipwall=WallBC_NoSlip_Adiabatic:new{},
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- 3. Simulation parameters.
config.flux_calculator= "ausmdv"
config.interpolation_order = 2
config.viscous = true
config.extrema_clipping = false
config.apply_limiter = false
config.apply_heuristic_pressure_based_limiting = true
-- config.spatial_deriv_locn = 'vertices'
-- config.spatial_deriv_calc = 'divergence'

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-6,
   number_of_phases = 2,
   max_steps_in_initial_phases = { 15 },
   use_physicality_check = true,
   max_linear_solver_iterations = 20,
   total_snapshots = 3,
   steps_between_status = 1,
   steps_between_snapshots = 5,
   steps_between_diagnostics = 1
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.9,
   start_cfl = 2.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 0.9
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 10.0
}
