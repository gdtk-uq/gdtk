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

config.title = "Schetz's Mach 4 laminar flow over a flat plate"
print(config.title)
config.viscous = true
config.include_ghost_cells_in_spatial_deriv_clouds = true
config.spatial_deriv_calc = "least_squares"
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "euler"
config.max_time = 2.4e-3  -- will allow 3 flow lengths   
config.dt_plot =  0.4e-3
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

--           wall
--        c---------b
-- flow=> |         |
--        d         |
--          -\-     |
--    flow=>    -\- |
--        0         a ----> x
-- 
a = Vector3:new{x=L}; b = Vector3:new{x=L, y=H}
c = Vector3:new{y=H}; d = Vector3:new{y=3.0*H/4.0}
patch = makePatch{north=Line:new{p0=c, p1=b}, east=Line:new{p0=a, p1=b},
		  south=Line:new{p0=d, p1=a}, west=Line:new{p0=d, p1=c}}
clusterx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
clustery_e = RobertsFunction:new{end0=false,end1=true,beta=1.016}
clustery_w = RobertsFunction:new{end0=false,end1=true,beta=1.05}
grd = StructuredGrid:new{psurface=patch,
			 cfList = {north=clusterx, east=clustery_e,
				   south=clusterx, west=clustery_w},
			 niv=221, njv=193}

-- create a special boundary condition for the no_slip_fixed_T wall that doesn't reference KOmegaWall
LaminarWallBC = WallBC_NoSlip_FixedT:new{Twall=300.0}
table.remove(LaminarWallBC.preSpatialDerivAction, 5)

-- Assemble the block from the grid and boundary data.
blks = FluidBlockArray{grid=grd, nib=22, njb=2, 
		       initialState=inflow,
		       bcList={north=LaminarWallBC,
			       east=OutFlowBC_Simple:new{},
			       south=InFlowBC_Supersonic:new{flowState=inflow},
			       west=InFlowBC_Supersonic:new{flowState=inflow}}}

-- convert structured-grid blocks to unstructured-grid blocks
for i=1,44 do
   SBlock2UBlock(blocks[i])
end

