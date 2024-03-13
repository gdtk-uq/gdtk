-- steady.lua
-- 2024-03-09 PJ, RG and KD
print("Model of Mohammadian's convex-ramp experiment with thermal nonequilibrium.")
config.solver_mode = 'steady'
config.dimensions = 2

p_inf = 66.43 -- Pa
u_inf = 1589.8 -- m/s
T_inf = 41.92 -- degree K
T_vib = 1000.0 -- freestream has frozen vibrational energy
T_wall = 296.0 -- degree K, assumed cold-wall temperature
--
nsp, nmodes = setGasModel('air-5sp-2T.gas')
print('5-species, 2T air model: nsp= ', nsp, ' nmodes= ', nmodes)
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_vib,}, velx=u_inf,
                       massf={N2=0.767,O2=0.233}}
flowDict = {
   initial=inflow,
   inflow=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_FixedPT:new{p_outside=p_inf/5, T_outside=T_inf},
   noslipwall=WallBC_NoSlip_FixedT:new{Twall=T_wall, group='loads'},
}
--
makeFluidBlocks(bcDict, flowDict)
mpiDistributeBlocks{ntasks=8}
--
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.viscous = true
config.thermo_interpolator = "rhop"
config.extrema_clipping = false
config.spatial_deriv_locn = "cells"
config.spatial_deriv_calc = "least_squares"

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   stop_on_relative_residual = 1.0e-6,
   number_of_phases = 2,
   phase_changes_at_steps = { 170 },
   inviscid_cfl_only = true,
   use_physicality_check = true,
   max_linear_solver_iterations = 50,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 3,
   steps_between_status = 5,
   steps_between_snapshots = 10,
   steps_between_diagnostics = 1,
   write_loads = true,
   steps_between_loads_update = 20,
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   frozen_preconditioner = true,
   steps_between_preconditioner_update = 5,
   use_adaptive_preconditioner = false,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.9,
   start_cfl = 5.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 0.9
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 5.0
}
