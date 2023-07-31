-- t5ns.lua
-- Reacting 2-Temperature air flow T5 nozzle starting process
-- HH 2023-05-30, modified with transient throat condition

config.title = "T5 nozzle starting process shot 2946. Transient throat."
print(config.title)
config.dimensions = 2
config.axisymmetric = true 
config.use_viscosity_from_cells=true
config.adjust_invalid_cell_data = true
config.max_invalid_cells = 10
config.viscous = true
config.turbulence_model = "spalart_allmaras"
config.flux_calculator = "adaptive"
config.spatial_deriv_calc = "divergence"
config.spatial_deriv_locn = "faces"
config.gasdynamic_update_scheme = "classic-rk3"
config.viscous_signal_factor = 0.4

T_wall = 300  --K

setGasModel('gas.lua')
config.reactions_file = 'chem.lua'
config.energy_exchange_file= 'kinetics.lua'
config.interpolation_order = 2

--
initial = FlowState:new{p=1.0, T=300., T_modes={300.0}, velx=0.0,  massf={N2=0.767, O2=0.233}}
--  inflow = InFlowBC_Transient:new{fileName = throat-file.dat}

config.reacting = true

print("Second, set up the simulation grid and the initial flow state.")
-- Estimate turbulence quantities for the inflow.
-- turb_intensity = 0.05
-- turb_to_laminar_visc_ratio = 100.0
-- throat_tke = 1.5 * (1800. * turb_intensity)^2
-- throat_mu_t = turb_to_laminar_visc_ratio * 3.0e-5
-- throat_omega = 4.0 * throat_tke / throat_mu_t
-- print("Inflow turbulence: tke=", throat_tke, "omega=", throat_omega)

-- Supersonic expansion is fully defined by its north edge
-- nozzle_profile = ArcLengthParameterizedPath:new{underlying_path=Bezier:new{points=bez_points}}
nozzle_profile = Spline2:new{filename="t5-profile.dat"}
exp_region = NozzleExpansionPatch:new{north=nozzle_profile}
-- Define structured grids for both regions.
x_cf_throat = RobertsFunction:new{end0=true, end1=true, beta=1.53}
x_cf = RobertsFunction:new{end0=true, end1=true, beta=1.2} --1.1
y_cf = RobertsFunction:new{end0=false, end1=true, beta=1.01}
-- throat_grid = StructuredGrid:new{psurface=throat_region, niv=11, njv=25,
--                                 cfList={west=y_cf, east=y_cf,
--                                         south=x_cf_throat, north=x_cf_throat}}
exp_grid = StructuredGrid:new{psurface=exp_region, niv=161, njv=31,
                              cfList={west=y_cf, east=y_cf,
                                      south=x_cf, north=x_cf}}
-- Divide the full domain into many blocks, mainly in the axial direction.
-- throat_fba = FBArray:new{grid=throat_grid, initialState=initial, label="throat",
--                         bcList={west=InFlowBC_Transient:new{fileName = "throat-file.dat"},
--                                 north=WallBC_NoSlip_FixedT:new{Twall=300.0, wall_function=true}},
--                         nib=1, njb=2}
exp_fba = FBArray:new{grid=exp_grid, initialState=initial, label="expansion",
                      bcList={west=InFlowBC_Transient:new{fileName = "throat-file.dat"},
		              north=WallBC_NoSlip_FixedT:new{Twall=300.0, wall_function=true},
                              east=OutFlowBC_Simple:new{}},
                      nib=16, njb=3}
identifyBlockConnections()
x_tr = 0.425; x1 = 9.856209152470000179e-01 
y0 = 0.0; y1 = 1.570611508610000040e-01
turbZone_botLeftCorner = Vector3:new{x=x_tr, y=y0}
turbZone_topRightCorner = Vector3:new{x=x1, y=y1}
TurbulentZone:new{p0=turbZone_botLeftCorner, p1=turbZone_topRightCorner}

-- Set a few more config options
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.max_time = 4.0e-3
config.max_step = 40000000
config.cfl_value = 0.7
config.dt_init = 1.0e-10
config.dt_plot = 1e-4
setHistoryPoint{x=99.0,y=0.1}
setHistoryPoint{x=99.0,y=0.15}
setHistoryPoint{x=99.0,y=0.05}
config.dt_history = 10.0e-06
