-- steady.lua
-- 2024-03-09 PJ, RG and KD
print("Model of Mohammadian's convex-ramp experiment with thermal nonequilibrium.")
config.solver_mode = 'steady'
config.dimensions = 2
print("Set up steady-state solve for Hakkinen et al's 1959 experiment.")
--
nsp, nmodes = setGasModel('air-5sp-2T.gas')
print('5-species, 2T air model: nsp= ', nsp, ' nmodes= ', nmodes)
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_vib,}, velx=u_inf,
                       massf={N2=0.767,O2=0.233}}
initial = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf,}, velx=0,
                        massf={N2=0.767,O2=0.233}}
flowDict = {
   initial=initial,
   inflow=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_FixedPT:new{p_outside=p_inf/5, T_outside=T_inf},
   noslipwall=WallBC_NoSlip_FixedT:new{Twall=T_wall},
}
--
makeFluidBlocks(bcDict, flowDict)
--
config.flux_calculator= "ausmdv"
config.interpolation_order = 2
config.viscous = true
-- config.spatial_deriv_locn = 'vertices'
-- config.spatial_deriv_calc = 'divergence'

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-6,
   number_of_phases = 2,
   phase_changes_at_steps = { 15 },
   use_physicality_check = true,
   max_linear_solver_iterations = 10,
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
