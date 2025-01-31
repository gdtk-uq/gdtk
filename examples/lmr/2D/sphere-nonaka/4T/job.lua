print("Sphere fired into air.")
gproGrid = "gridpro/sphere-orig.grd"
--
no_flow_times = 5.0
njb = 4
Db = 14.0e-3
Rc = Db/2
setGasModel('air-5sp-4T.gas')
p_inf = 4850.0 -- Pa
T_inf = 293.0 -- K
u_inf = 3490.0 -- m/s
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf, T_inf, T_inf}, velx=u_inf,
                       massf={N2=0.767,O2=0.233}}
flowDict = {
   initial=inflow,
   inflow=inflow
}
--
body_flow_time = Db/u_inf
config.solver_mode = 'transient'
config.reacting = true
config.reactions_file = 'air-5sp-6r-4T.chem'
config.energy_exchange_file = 'air-energy-exchange-4T.exch'
config.dimensions = 2
config.axisymmetric = true
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.interpolation_delay = 3*body_flow_time
config.grid_motion = "shock_fitting"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.shock_fitting_delay = body_flow_time
config.max_invalid_cells = 20
config.adjust_invalid_cell_data = true
config.max_time = no_flow_times*body_flow_time
config.max_step = 80000
config.dt_init = 1.0e-10
config.cfl_value = 0.75
config.dt_plot = config.max_time/no_flow_times
config.allow_reconstruction_for_species = false
config.allow_reconstruction_for_energy_modes = false
config.apply_limiter = false
config.extrema_clipping = false
--
grids = importGridproGrid(gproGrid, 1.0)
registerFluidGridArray{
   grid=grids[1],
   nib=1, njb=njb,
   fsTag='initial',
   shock_fitting=true,
   bcTags={west='inflow_sf', north='outflow', east='wall'}
}
--
bcDict = {
   inflow_sf=InFlowBC_ShockFitting:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{},
   wall=WallBC_WithSlip0:new{}
}
--
makeFluidBlocks(bcDict, flowDict)
--
