print("Circular-arc foil in transonic flow, transient simulation.")
-- PJ, 2024-03-26, adpated from sharp-cone example

config.solver_mode = "steady"

setGasModel('ideal-air.gas')
p_stag = 101.3e3 -- Pa
T_stag = 300.0 -- degree K
stagGas = FlowState:new{p=p_stag, T=T_stag}
M_inf = 0.84 -- for supercritical flow over foil
p_inf = p_stag/idealgasflow.p0_p(M_inf)
T_inf = T_stag/idealgasflow.T0_T(M_inf)
a_inf = math.sqrt(1.4*287.1*T_inf)
velx_inf = M_inf * a_inf
print("Free-stream: p=", p_inf, "T=", T_inf, "velx=", velx_inf)
initGas = FlowState:new{p=p_inf, T=T_inf, velx=velx_inf}

flowDict = {
   initial=initGas,
   inflow=stagGas
}
bcDict = {
   inflow=InFlowBC_FromStagnation:new{stagnationState=stagGas},
   outflow=OutFlowBC_FixedP:new{p_outside=p_inf}
}
makeFluidBlocks(bcDict, flowDict)

config.flux_calculator= "ausmdv"
config.interpolation_order = 2

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 3,
   stop_on_relative_residual = 1.0e-6,
   number_of_phases = 2,
   max_steps_in_initial_phases = { 15 },
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
