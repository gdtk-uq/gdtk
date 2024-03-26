print("Circular-arc foil in transonic flow, transient simulation.")
-- PJ, 2017-10-27, adapted from channel-with-bump example
--     2024-03-26, adpated for Eilmer5

nsp, nmodes = setGasModel('ideal-air.gas')
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

config.solver_mode = "transient"
config.flux_calculator = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 0.8
config.max_time = 0.400 -- 100.0*L/velx_inf=0.366
config.max_step = 50000
config.dt_init = 1.0e-6
config.dt_plot = config.max_time/40
