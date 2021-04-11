config.title = "Sphere fired into air."
print(config.title)
gproGrid = "gridpro/sphere-orig.grd"
no_flow_times = 10.0

njb = 4
Db = 14.0e-3
Rc = Db/2
setGasModel('air-5sp-2T.gas')
p_inf = 4850.0 -- Pa
T_inf = 293.0 -- K
u_inf = 3490.0 -- m/s
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf}, velx=u_inf,
                       massf={N2=0.767,O2=0.233}}

body_flow_time = Db/u_inf
config.reacting = true
config.reactions_file = 'air-5sp-6r-2T.chem'
config.energy_exchange_file = 'air-VT.exch'
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
config.dt_plot = config.max_time/100

grids = importGridproGrid(gproGrid, 1.0)
blk = FBArray:new{grid=grids[1], initialState=inflow, label='blk',
                  bcList={west=InFlowBC_ShockFitting:new{flowState=inflow},
                          north=OutFlowBC_Simple:new{}},
                  nib=1, njb=njb}

-- Set some recording points to capture flow data
setHistoryPoint{x=0.0, y=0.0}
setHistoryPoint{x=Rc, y=Rc}
setHistoryPoint{ib=0, i=0, j=0}
config.dt_history = 1.0e-7
