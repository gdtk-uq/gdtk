-- ==========================================================
-- Title / Description
-- ==========================================================

-- Simulation builds on work in examples/3D/moving-grid/flying-cube/
-- Displays a free flying rotating square in M6 freestream conditions

config.title = "Moving Square Projectile"
config.dimensions = 2
config.axisymmetric = false

-- ==========================================================
-- Boundary Condition Pre-Calcs
-- ==========================================================

nsp,nmodes,gasModel = setGasModel("ideal-air-gas-model.lua")

init_state = FlowState:new{p = 500.0, T =  300.0}
inflow_state = FlowState:new{p = 760.0, T = 71.0, velx = 1005.0} -- nominal TUSQ M6 nozzle conditions

InOutBC = InOutFlowBC_Ambient:new{flowState = inflow_state} -- Boundary condition applies an inflow condition
															-- for flow coming in to the domain, and a simple
															-- 1D flux extrapolation for flow leaving the domain.
															-- Subsonic outflow can cause issues with this BC.

WallBC  = WallBC_WithSlip1:new{group = "walls"} -- Using this BC is important as it doesn't rely on
												-- ghost cells. The ghost cell implementation
												-- does not account for the wall velocity when reflecting
												-- the velocity vector.

-- ========================================================
-- Geometry, grid and block setup.
-- ========================================================

-- Geometry is a square with corner radius r_corner.
-- Apologies for the messy grid setup.

-- Geometry constants.

Rx = math.cos(math.pi/4) 
Ry = math.sin(math.pi/4)

a_sq = 0.01 -- 10 mm

SF_o = 4 -- a_outter = SF_o * a_sq
SF_c = 1 -- cell scaling factor

ncells1 = math.floor(SF_c * 21)
ncells2 = math.floor(SF_c * 24)
ncells3 = math.floor(SF_c * 3)

r_corner = a_sq/10

C1 = Vector3:new{x=-a_sq,y=-a_sq}
C2 = Vector3:new{x=-a_sq,y=a_sq}
C3 = Vector3:new{x=a_sq,y=a_sq}
C4 = Vector3:new{x=a_sq,y=-a_sq}

c1p1 = C1 + Vector3:new{x=r_corner,y=0}
c1p2 = C1 + Vector3:new{x=0,y=r_corner} 

c2p1 = C2 + Vector3:new{x=0,y=-r_corner}
c2p2 = C2 + Vector3:new{x=r_corner,y=0}

c3p1 = C3 + Vector3:new{x=-r_corner,y=0}
c3p2 = C3 + Vector3:new{x=0,y=-r_corner}

c4p1 = C4 + Vector3:new{x=0,y=r_corner}
c4p2 = C4 + Vector3:new{x=-r_corner,y=0}

a1 = SF_o * c1p2
b1 = c1p2
c1 = c1p1
d1 = SF_o * c1p1

a2 = SF_o * c2p1
b2 = c2p1
c2 = c1p2
d2 = a1

a3 = SF_o * c2p2
b3 = c2p2
c3 = c2p1
d3 = a2

a4 = a3
b4 = SF_o * c3p1
c4 = c3p1
d4 = b3

a5 = b4
b5 = SF_o * c3p2
c5 = c3p2
d5 = c4

a6 = b5
b6 = SF_o * c4p1
c6 = c4p1
d6 = c5

a7 = b6
b7 = SF_o * c4p2
c7 = c4p2
d7 = c6

a8 = b7
b8 = d1
c8 = c1p1
d8 = c4p2

a1d1 = Arc:new{p0 = d1, p1 = a1, centre = Vector3:new{x=0,y=0,z=0}} -- W
a1b1 = Line:new{p0 = a1, p1 = b1} -- N
c1b1 = Arc:new{p0 = c1, p1 = b1, centre = C1 + Vector3:new{x=r_corner,y=r_corner}} -- E
d1c1 = Line:new{p0 = d1, p1 = c1} --S

a2d2 = Arc:new{p0 = d2, p1 = a2, centre = Vector3:new{x=0,y=0,z=0}} -- W
a2b2 = Line:new{p0 = a2, p1 = b2} -- N
c2b2 = Line:new{p0 = c2, p1 = b2} -- E
d2c2 = Line:new{p0 = d2, p1 = c2} --S

a3d3 = Arc:new{p0 = d3, p1 = a3, centre = Vector3:new{x=0,y=0,z=0}} -- W
a3b3 = Line:new{p0 = a3, p1 = b3} -- N
c3b3 = Arc:new{p0 = c3, p1 = b3, centre = C2 + Vector3:new{x=r_corner,y=-r_corner}} -- E
d3c3 = Line:new{p0 = d3, p1 = c3} --S

d4a4 = Line:new{p0 = d4, p1 = a4} -- W
a4b4 = Arc:new{p0 = a4, p1 = b4, centre = Vector3:new{x=0,y=0,z=0}} -- N
d4c4 = Line:new{p0 = d4, p1 = c4} -- S
c4b4 = Line:new{p0 = c4, p1 = b4} -- E

d5a5 = Line:new{p0 = d5, p1 = a5} -- W
a5b5 = Arc:new{p0 = a5, p1 = b5, centre = Vector3:new{x=0,y=0,z=0}} -- N
d5c5 = Arc:new{p0 = d5, p1 = c5, centre = C3 + Vector3:new{x=-r_corner,y=-r_corner}} -- S
c5b5 = Line:new{p0 = c5, p1 = b5} -- E

d6a6 = Line:new{p0 = d6, p1 = a6} -- N
b6a6 = Arc:new{p0 = b6, p1 = a6, centre = Vector3:new{x=0,y=0,z=0}} -- W
c6b6 = Line:new{p0 = c6, p1 = b6} -- S
c6d6 = Line:new{p0 = c6, p1 = d6} -- E

d7a7 = Line:new{p0 = d7, p1 = a7} -- N
b7a7 = Arc:new{p0 = b7, p1 = a7, centre = Vector3:new{x=0,y=0,z=0}} -- W
c7b7 = Line:new{p0 = c7, p1 = b7} -- S
c7d7 = Arc:new{p0 = c7, p1 = d7, centre = C4 + Vector3:new{x=-r_corner,y=r_corner}} -- E

d8a8 = Line:new{p0 = d8, p1 = a8} -- N
b8a8 = Arc:new{p0 = b8, p1 = a8, centre = Vector3:new{x=0,y=0,z=0}} -- E
c8b8 = Line:new{p0 = c8, p1 = b8} -- S
c8d8 = Line:new{p0 = c8, p1 = d8} -- W

grid1 = StructuredGrid:new{
	psurface = makePatch{north = a1b1, east = c1b1, south = d1c1, west = a1d1},
	cfList = {north = RobertsFunction:new{end0=false, end1=true, beta=1.15},
	          south = RobertsFunction:new{end0=false, end1=true, beta=1.15}},
	niv = ncells1 + 1, njv = ncells3 + 1}

grid2 = StructuredGrid:new{
	psurface = makePatch{north = a2b2, east = c2b2, south = d2c2, west = a2d2},
	cfList = {north = RobertsFunction:new{end0=false, end1=true, beta=1.15},
	          south = RobertsFunction:new{end0=false, end1=true, beta=1.15}},
	niv = ncells1 + 1, njv = ncells2 + 1}
	
grid3 = StructuredGrid:new{
	psurface = makePatch{north = a3b3, east = c3b3, south = d3c3, west = a3d3},
	cfList = {north = RobertsFunction:new{end0=false, end1=true, beta=1.15},
	          south = RobertsFunction:new{end0=false, end1=true, beta=1.15}},
	niv = ncells1 + 1, njv = ncells3 + 1}	
	
grid4 = StructuredGrid:new{
	psurface = makePatch{north = a4b4, east = c4b4, south = d4c4, west = d4a4},
	cfList = {east = RobertsFunction:new{end0=true, end1=false, beta=1.15},
	          west = RobertsFunction:new{end0=true, end1=false, beta=1.15}},
	niv = ncells2 + 1, njv = ncells1 + 1}

grid5 = StructuredGrid:new{
	psurface = makePatch{north = a5b5, east = c5b5, south = d5c5, west = d5a5},
	cfList = {east = RobertsFunction:new{end0=true, end1=false, beta=1.15},
	          west = RobertsFunction:new{end0=true, end1=false, beta=1.15}},
	niv = ncells3 + 1, njv = ncells1 + 1}

grid6 = StructuredGrid:new{
	psurface = makePatch{north = d6a6, east = b6a6, south = c6b6, west = c6d6},
	cfList = {north = RobertsFunction:new{end0=true, end1=false, beta=1.15},
	          south = RobertsFunction:new{end0=true, end1=false, beta=1.15}},
	niv = ncells1 + 1, njv = ncells2 + 1}

grid7 = StructuredGrid:new{
	psurface = makePatch{north = d7a7, east = b7a7, south = c7b7, west = c7d7},
	cfList = {north = RobertsFunction:new{end0=true, end1=false, beta=1.15},
	          south = RobertsFunction:new{end0=true, end1=false, beta=1.15}},
	niv = ncells1 + 1, njv = ncells3 + 1}

grid8 = StructuredGrid:new{
	psurface = makePatch{north = d8a8, east = b8a8, south = c8b8, west = c8d8},
	cfList = {north = RobertsFunction:new{end0=true, end1=false, beta=1.15},
	          south = RobertsFunction:new{end0=true, end1=false, beta=1.15}},
	niv = ncells1 + 1, njv = ncells2 + 1}

-- Now create our blocks

blk1 = FluidBlock:new{grid = grid1, initialState = init_state, 
	                  bcList = {west = InOutBC, east = WallBC}}

blk2 = FluidBlock:new{grid = grid2, initialState = init_state, 
	                  bcList = {west = InOutBC, east = WallBC}}

blk3 = FluidBlock:new{grid = grid3, initialState = init_state, 
	                  bcList = {west = InOutBC, east = WallBC}}

blk4 = FluidBlock:new{grid = grid4, initialState = init_state,
                      bcList = {north = InOutBC, south = WallBC}}

blk5 = FluidBlock:new{grid = grid5, initialState = init_state,
                      bcList = {north = InOutBC, south = WallBC}}

blk6 = FluidBlock:new{grid = grid6, initialState = init_state,
	                  bcList = {east = InOutBC, west = WallBC}}

blk7 = FluidBlock:new{grid = grid7, initialState = init_state,
	                  bcList = {east = InOutBC, west = WallBC}}

blk8 = FluidBlock:new{grid = grid8, initialState = init_state,
                      bcList = {east = InOutBC, west = WallBC}}

identifyBlockConnections()

-- ==========================================================
-- Solver settings
-- ==========================================================

config.interpolation_order = 2
config.dt_init = 1.0e-7
config.fixed_time_step = false
config.max_time = 2.0e-3
config.max_step = 500000
config.dt_plot = config.max_time / 20
config.cfl_schedule = {{0.0, 0.15},{config.max_time/5.0, 0.20},{config.max_time/2.0, 0.25}}

-- Flux Calculator

-- The main challenge with this calculation is pressure instabilities from the strong boundary shock.
-- These instabilities affect the pressure distribution over the bodies, and hence aerodynamic forces.

-- Flux blending trigger settings.
config.flux_calculator = "adaptive_hanel_ausmdv"
config.compression_tolerance = -0.25
--config.apply_heuristic_pressure_based_limiting = true -- setting requires some work to implement properly
config.shear_tolerance = 0.25

config.adjust_invalid_cell_data = true
config.max_invalid_cells = 10

config.write_loads = true
config.dt_loads = 1e-4
config.boundary_groups_for_loads = "walls"

-- Moving mesh related parameters

config.gasdynamic_update_scheme = "moving-grid-2-stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "grid-motion.lua" -- lua function for prescribing grid motion (body dynamics)
config.udf_supervisor_file = 'udf-supervisor.lua' -- file for calculating projectile dynamics

run_time_loads = {{group = "walls", moment_centre = Vector3:new{x=0,y=0,z=0}}}
config.compute_run_time_loads = true

-- Setup the parameters required for our dynamics model
-- Required states are angle and angular velocity.

config.user_pad_length = 3

thetaz_init = 0 -- rad
rb_init = 2000 -- rad/s
step_init = 0

user_pad_data = {thetaz_init, rb_init, 
                 step_init}

