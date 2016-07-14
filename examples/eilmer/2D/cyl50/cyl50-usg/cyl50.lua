-- cyl50.lua
-- PJ, updated from Tcl script, 14-Aug-2006
--     Eilmer3 port, July 2008
--     Eilmer4 port, May 2015

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
sth = Line:new{p0=a,p1=b}; nth = Line:new{p0=d,p1=c};
wst = Line:new{p0=a,p1=d}; est = Line:new{p0=b,p1=c}
clusterx = RobertsFunction:new{end0=true,end1=false,beta=1.1}
clustery = RobertsFunction:new{end0=true,end1=false,beta=1.01}
grd = StructuredGrid:new{psurface=makePatch{north=nth, east=est, south=sth, west=wst, gridType="ao"},
			 cfList = {north=clusterx, east=clustery, south=clusterx, west=clustery},
			 niv=51, njv=51}

-- create special boundary condition for the no_slip_fixed_T wall
LaminarWallBC = WallBC_NoSlip_FixedT:new{Twall=222.0}
table.remove(LaminarWallBC.preSpatialDerivAction, 5)
-- Assemble the block from the grid and boundary data.
blks = SBlockArray{grid=grd, nib=2, njb=2, 
		   fillCondition=inflow,
		   bcList={north=InFlowBC_Supersonic:new{flowCondition=inflow},
			   east=OutFlowBC_Simple:new{},
			   south=LaminarWallBC,
			   west=InFlowBC_Supersonic:new{flowCondition=inflow}}}
-- convert structured blocks to unstructured blocks
for i=1,4 do
   SBlock2UBlock(blocks[i])
end

-- Do a little more setting of the simulation configuration data.
config.viscous = true
config.interpolation_order = 2
config.flux_calculator = "adaptive"
config.include_ghost_cells_in_spatial_deriv_clouds = true
config.spatial_deriv_calc = "least_squares"
config.gasdynamic_update_scheme = "euler"
config.max_time = 8.0e-3  -- seconds
config.max_step = 230000
config.dt_init = 3.0e-8
config.dt_plot = 4.0e-3/20.0
