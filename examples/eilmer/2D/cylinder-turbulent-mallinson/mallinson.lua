-- mallinson.lua : Mach 8.8 Turbulent flow over a cylinder
-- Dimir Y.X. Pot, Samuel J. Stennett, Wilson Y.K. Chan, 2018-03-14
-- Ported from Eilmer3, 2018-10-27
-- Note that the grid has been set up to reflect the 
-- "600x130cells-yplus0.2-ar577" configuration.
-- Reference:
--     Mallinson SG et al. (2000), Shock Waves, v.10, pp.313-322
--     Boyce RR & Hillier R (2000), AIAA 2000-2226

config.title = "Mallinson's Mach 9 flow over a hollow cylinder (k-omega)"
print(config.title)
config.dimensions = 2
config.axisymmetric = true     -- Axisymmetric calculations
config.turbulence_model = "k_omega"
config.viscous = true
config.flux_calculator = 'ausmdv'
config.gasdynamic_update_scheme = "classic-rk3"
--
config.max_time = 5.0e-3
config.dt_plot =  0.5e-3
config.dt_history = 1.0e-3
config.max_step = 30000000
--
config.cfl_value = 0.4
config.cfl_count = 3
config.stringent_cfl = false  -- true is more robust
config.dt_init = 1.0e-14      

-- Gas model and flow conditions to match Mabey's data set 74021802
nsp, nmodes, gm = setGasModel('ideal-N2-gas-model.lua')
p_inf = 3.3e3   -- Pa 
u_inf = 1498.0  -- m/s
T_inf = 69.7    -- K
T_wall = 295.0  -- cylinder surface temperature, K

-- Set inflow turbulence quantities to almost-laminar values so that
-- we can use the TurbulenceZone later to fix a transition location
-- on the surface of the cylinder.
tke_inf = 1e-12
omega_inf = 1.0

-- Set up flow state
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, tke=tke_inf, omega=omega_inf}
print("Inflow Check\n", inflow)
                       
-- It seems that everyone has been starting the flow field with 
-- inflow conditions throughout and then letting the boundary layer 
-- grow out into the main stream, instead of having an initial state 
-- that is more representative of the actual tunnel conditions (which
-- unfortunately seemed to wreck havoc with our turbulence model in e3).
--initial_flag = "inflow"   # "inflow" or "actual"
--if initial_flag == "inflow":
--    initial = inflow
--elif initial_flag == "actual":
--    initial = FlowState:new{p=0.1*p_inf, u=0.0, v=0.0, T=296.0,
--                            tke=tke_inf/100.0, omega=omega_inf/10.0}

-- Define geometry (all dimensions in metres)
x0 = 0.0; x1 = 0.85
y0 = 0.0; y1 = 0.0375; y2 = 0.057; y3 = 0.2 

-- Define approximately where boundary layer transitions to turbulence.
-- We set the location of the start of the TurbulentZone at 0.047 m to
-- get the transition location to match that of the experiment (0.085 m
-- downtream of the leading edge.
x_tr = 0.047    
turbZone_botLeftCorner = Vector3:new{x=x_tr, y=y1}
turbZone_topRightCorner = Vector3:new{x=x1, y=y3}
TurbulentZone:new{p0=turbZone_botLeftCorner, p1=turbZone_topRightCorner}

-- Start building grid
a = Vector3:new{x=x0, y=y1}; b = Vector3:new{x=x1, y=y1}
c = Vector3:new{x=x1, y=y3}; d = Vector3:new{x=x0, y=y2}
patch = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d}
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
cflist = {north=cfx, east=RobertsFunction:new{end0=true,end1=false,beta=1.0002},
          south=cfx, west=RobertsFunction:new{end0=true,end1=false,beta=1.0020}}
grd = StructuredGrid:new{psurface=patch, niv=601, njv=131, cfList=cflist}

-- Assemble the block from the grid and boundary data.
blks = FBArray:new{grid=grd, nib=4, njb=1, fillCondition=inflow, 
                       bcList = {north=UserDefinedBC:new{fileName='udf-supersonic-in.lua'}, 
                                 east=OutFlowBC_Simple:new{}, 
                                 south=WallBC_NoSlip_FixedT:new{Twall=T_wall, wall_function=false}, 
                                 west=UserDefinedBC:new{fileName='udf-supersonic-in.lua'}}}

identifyBlockConnections()
