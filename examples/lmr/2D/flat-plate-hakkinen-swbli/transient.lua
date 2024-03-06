-- transient.lua
-- 2024-03-06 PJ and RG
print("Set up a transient simulation of Hakkinen et al's 1959 experiment.")
config.solver_mode = 'transient'
config.dimensions = 2

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
flowDict = {
   initial=inflow,
   inflow=inflow
}
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_FixedPT:new{p_outside=p_inf, T_outside=T_inf},
   noslipwall=WallBC_NoSlip_Adiabatic:new{},
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
config.max_time = 5.0*L2/u_inf -- time in flow lengths
config.max_step = 200000
config.dt_init = 1.0e-8
config.dt_plot = config.max_time/10
