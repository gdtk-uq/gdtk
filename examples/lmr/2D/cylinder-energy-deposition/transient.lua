-- transient.lua
print("Riggins drag reduction calculation -- set up simulation")
dofile("parameters.lua")

-- Gas model and flow states
nsp, nmodes, gmodel = setGasModel("ideal-air.gas")
inflow = FlowState:new{p=p_inf, T=T_inf, velx=V_inf}
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0}
flowDict = {
   initial=initial,
   inflow=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   body=WallBC_NoSlip_Adiabatic:new{group="loads"},
   outflow=OutFlowBC_Simple:new{}
}
--
makeFluidBlocks(bcDict, flowDict)

-- Now, set some simulation options.
config.viscous = true
body_flow_length = R/V_inf
t_final = 20 * body_flow_length -- allow time to establish
config.flux_calculator = "adaptive_hanel_ausmdv"
config.interpolation_delay = 10 * body_flow_length
config.max_time = t_final
config.max_step = 40000
config.dt_init = 1.0e-9
config.cfl_value = 0.4
config.dt_plot = config.max_time/10.0
config.write_loads = true
config.boundary_groups_for_loads = "loads"
config.dt_loads = config.dt_plot

config.udf_source_terms = true
config.udf_source_terms_file = "energy-deposition.lua"
