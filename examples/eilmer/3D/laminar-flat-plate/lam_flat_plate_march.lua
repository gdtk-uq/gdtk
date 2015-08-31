-- lam_flat_plate_march.lua
-- Laminar flow over a flat plate 
--
-- This script creates and runs a simulation of a worked
-- example (Example 5-11 on page 155) in Joseph Schetz's
-- book "Boundary Layer Analysis".
--
-- The example involves a laminar Mach 4.0 flow over a
-- 1.0 m long flat plate. The boundary layer will grow on
-- the NORTH boundary. Simulation results will be compared 
-- with those from the CLBL code from Schetz's book.
--
-- Peter Jacobs & Wilson Chan, 
-- 18 Jan 2010 Python version for eilmer3
-- 07 Jun 2015 Lua version for Eilmer4                           

config.title = "Schetz's Mach 4 laminar flow over a flat plate 3D"
print(config.title)
config.dimensions = 3
config.viscous = true
config.flux_calc = "ausmdv"
config.gasdynamic_update_scheme = "euler"
config.max_time = 2.4e-3  -- will allow 3 flow lengths   
config.dt_plot =  config.max_time
config.dt_history = 1.0e-5
config.max_step = 300000
config.cfl_target = 0.4
config.cfl_count = 3
config.dt_init = 1.0e-9  -- only an initial guess - the simulation will take over

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)

-- Define flow conditions to match Schetz's worked example 5-1
p_inf = 1.013e3  -- Pa
u_inf = 1390.0   -- m/s
T_inf = 300.0    -- degrees K

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0}

-- Geometry of plate and flow domain.
-- Have set up the flat plate 0.1 m longer than the actual plate
-- so that we don't have to use the simulation profiles right on 
-- the boundary.
L = 1.1       -- Length of flat plate (in metres)
H = 0.4 * L   -- Height of flow domain (in metres)
dz = 0.010    -- depth in z-direction (also k-index direction)

-- In 3D, look down on the top surface
--           wall
--       p01--------p11
-- flow=> |         |
--       p00        |
--          -\-     |
--    flow=>    -\- |
--        0        p10 ----> x
-- 
p000 = Vector3:new{0.0, 3.0*H/4.0, 0.0}
p100 = Vector3:new{L, 0.0, 0.0}
p110 = Vector3:new{L, H, 0.0}
p010 = Vector3:new{0.0, H, 0.0}
p001 = Vector3:new{0.0, 3.0*H/4.0, dz}
p101 = Vector3:new{L, 0.0, dz}
p111 = Vector3:new{L, H, dz}
p011 = Vector3:new{0.0, H, dz}

clusterx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
clustery_e = RobertsFunction:new{end0=false,end1=true,beta=1.016}
clustery_w = RobertsFunction:new{end0=false,end1=true,beta=1.05}
clusterz = LinearFunction:new{}
cflist = {clusterx, clustery_e, clusterx, clustery_w,
	  clusterx, clustery_e, clusterx, clustery_w,
	  clusterz, clusterz, clusterz, clusterz}
vol = TFIVolume:new{vertices={p000,p100,p110,p010,p001,p101,p111,p011}}
grd = StructuredGrid:new{pvolume=vol, cfList=cflist, niv=221, njv=193, nkv=3}

-- Assemble the block from the grid and boundary data.
bcList={north=FixedTWallBC:new{Twall=300.0},
	east=ExtrapolateOutBC:new{},
	south=SupInBC:new{flowCondition=inflow},
	west=SupInBC:new{flowCondition=inflow},
	bottom=SlipWallBC:new{},
	top=SlipWallBC:new{}}
blks = SBlockArray{grid=grd, nib=22, njb=2, nkb=1,
		   bcList=bcList, fillCondition=inflow}

config.block_marching = true
config.nib = 22
config.njb = 2
config.nkb = 1
