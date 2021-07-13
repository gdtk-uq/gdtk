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
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "euler"
config.max_time = 10.0e-3  -- will allow lots more than 3 flow lengths
config.dt_plot =  config.max_time
config.dt_history = 1.0e-5
config.max_step = 300000
config.cfl_value = 0.4
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
p000 = Vector3:new{x=0.0, y=3.0*H/4.0, z=0.0}
p100 = Vector3:new{x=L,   y=0.0,       z=0.0}
p110 = Vector3:new{x=L,   y=H,         z=0.0}
p010 = Vector3:new{x=0.0, y=H,         z=0.0}
p001 = Vector3:new{x=0.0, y=3.0*H/4.0, z=dz}
p101 = Vector3:new{x=L,   y=0.0,       z=dz}
p111 = Vector3:new{x=L,   y=H,         z=dz}
p011 = Vector3:new{x=0.0, y=H,         z=dz}

clusterx = RobertsFunction:new{end0=true, end1=false, beta=1.05}
clustery_e = RobertsFunction:new{end0=false, end1=true, beta=1.016}
clustery_w = RobertsFunction:new{end0=false, end1=true, beta=1.05}
clusterz = LinearFunction:new{}
cflist = {edge01=clusterx, edge12=clustery_e, edge32=clusterx, edge03=clustery_w,
	  edge45=clusterx, edge56=clustery_e, edge76=clusterx, edge47=clustery_w,
	  edge04=clusterz, edge15=clusterz,   edge26=clusterz, edge37=clusterz}
vol = TFIVolume:new{vertices={p000,p100,p110,p010,p001,p101,p111,p011}}
grd = StructuredGrid:new{pvolume=vol, cfList=cflist, niv=221, njv=193, nkv=5}

-- Assemble the block from the grid and boundary data.
bcList={north=WallBC_NoSlip_FixedT:new{Twall=300.0},
	east=OutFlowBC_Simple:new{},
	south=InFlowBC_Supersonic:new{flowState=inflow},
	west=InFlowBC_Supersonic:new{flowState=inflow},
	bottom=WallBC_WithSlip:new{},
	top=WallBC_WithSlip:new{}}
blks = FBArray:new{grid=grd, nib=22, njb=4, nkb=1,
                   bcList=bcList, initialState=inflow}

config.block_marching = true
config.nib = 22
config.njb = 4
config.nkb = 1
