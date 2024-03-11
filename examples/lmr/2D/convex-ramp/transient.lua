-- transient.lua
-- 2024-03-09 PJ, RG and KD
print("Model of Mohammadian's convex-ramp experiment with thermal nonequilibrium.")
config.solver_mode = 'transient'
config.dimensions = 2

p_inf = 66.43 -- Pa
u_inf = 1589.8 -- m/s
T_inf = 41.92 -- degree K
T_vib = 1000.0 -- freestream has frozen vibrational energy
T_wall = 296.0 -- degree K, assumed cold-wall temperature

nsp, nmodes = setGasModel('air-5sp-2T.gas')
print('5-species, 2T air model: nsp= ', nsp, ' nmodes= ', nmodes)
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_vib,}, velx=u_inf,
                       massf={N2=0.767,O2=0.233}}
initial = FlowState:new{p=p_inf/5, T=T_inf, T_modes={T_inf,}, velx=0,
                        massf={N2=0.767,O2=0.233}}
flowDict = {
   initial=initial,
   inflow=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_FixedPT:new{p_outside=p_inf/5, T_outside=T_inf},
   noslipwall=WallBC_NoSlip_FixedT:new{Twall=T_wall, group='loads'},
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- mpiDistributeBlocks{ntasks=4}
config.gasdynamic_update_scheme = 'classic-rk3'
config.flux_calculator = 'adaptive'
config.viscous = true
config.spatial_deriv_locn = 'vertices'
config.spatial_deriv_calc = 'divergence'
config.cfl_value = 1.0
L2 = 0.254 -- End of ramp, and of the whole flow domain.
config.max_time = 5.0*L2/u_inf -- time in flow lengths
config.max_step = 200000
config.dt_init = 1.0e-8
config.dt_plot = config.max_time/10
config.write_loads = true
config.dt_loads = config.max_time/50
