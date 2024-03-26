-- Rowan J. Gollan and Peter J: 2024-03-26, Ported to Eilmer5
--
-- This script is used to setup a simlulation of Lehr's
-- hemispherical projectile fired into a detonable gas.
--
-- Reference:
--   Lehr, H. (1972)
--   Acta Astronautica, 17, pp.589--597
--
print("Lehr experiment M=3.55 -- set up grid")
R = 7.5e-3 -- nose radius, metres
config.axisymmetric = true

-- Free stream conditions
-- taken from Table 1, Entry 5.
nsp, nmodes, gmodel = setGasModel('h2o2.gas')
p_inf = 186.0/760.0*101325.0 -- Pa
u_inf = 1892 -- m/s
T_inf = 292 -- K
molef_inf = {H2=2/3, O2=1/3}
massf_inf = gmodel:molef2massf(molef_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, massf=massf_inf}
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0, massf=massf_inf}
flowDict = {
   initial=initial,
   inflow=inflow
}
bcDict = {
   inflow_sf=InFlowBC_ShockFitting:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
--
makeFluidBlocks(bcDict, flowDict)

-- Now, set some configuration options.
body_flow_length = R/u_inf
t_final = 20 * body_flow_length -- allow time to establish
config.reacting = true
config.reactions_file = 'h2o2.chem'
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 3 * body_flow_length
config.interpolation_delay = 10 * body_flow_length
config.max_time = t_final
config.max_step = 800000
config.dt_init = 1.0e-10
config.cfl_value = 0.4
config.dt_plot = config.max_time/100.0
