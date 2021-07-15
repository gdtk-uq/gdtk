-- cyl50.lua
-- PJ, updated from Tcl script, 14-Aug-2006
--     Eilmer3 port, July 2008
--     Eilmer4 port, May 2015
--     Change to using implicit update, July 2021

config.title = "Mach 2 flow along the axis of a 5mm cylinder."
print(config.title)
config.axisymmetric = true

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
inflow = FlowState:new{p=257.3, T=222.0, velx=597.3, vely=0.0}

-- Set up a quadrilateral grid in the (x,y)-plane.
--  y     c
--  ^   / |
--  | /   |
--  d     |
--  a-----b
--  0------------> x
a = Vector3:new{x=0.0,y=0.005}; b = Vector3:new{x=1.0,y=0.005}
c = Vector3:new{x=1.0,y=0.7}; d = Vector3:new{x=0.0,y=0.06}
clusterx = RobertsFunction:new{end0=true,end1=false,beta=1.1}
clustery = RobertsFunction:new{end0=true,end1=false,beta=1.01}
myCFList = {north=clusterx, east=clustery, south=clusterx, west=clustery}
grd = StructuredGrid:new{psurface=AOPatch:new{p00=a, p10=b, p11=c, p01=d},
			 cfList=myCFList, niv=51, njv=51}

-- Assemble the blocks from the grid and boundary data.
blks = FBArray:new{grid=grd, nib=2, njb=2,
                   initialState=inflow,
                   bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
                           east=OutFlowBC_Simple:new{},
                           south=WallBC_NoSlip_FixedT:new{Twall=222.0},
                           west=InFlowBC_Supersonic:new{flowState=inflow}}}

-- Do a little more setting of the simulation configuration data.
config.viscous = true
config.spatial_deriv_locn = "vertices"
config.spatial_deriv_calc = "divergence"
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "backward_euler"
config.max_time = 8.0e-3  -- seconds
config.max_step = 230000
config.dt_init = 3.0e-8
config.cfl_schedule_times = {0.0, 10.0e-6, 100.0e-6}
config.cfl_schedule_values = {0.5, 5.0, 50.0}
config.dt_plot = 4.0e-3
