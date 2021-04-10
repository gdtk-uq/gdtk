-- Author: Rowan J. Gollan
-- Date: 2017-01-04
--
-- Modified & Adapted for new test case
-- Author: Fabian Zander
-- Date: 2019-05-02
--
-- This Lua script is used to prepare an Eilmer4 simulation of a projectile.
-- This simulation makes use of the user-defined moving grid and user-defined
-- boundary condition capabilities in Eilmer4.
--
-- This script has been adapted as a validation case against the
-- projectile case in the 16th AFMC paper (2007)
-- "Development of Casbar: a Two-phase Flow Code for the Interior Ballistics
-- Problem", R.J Gollan, I.A. Johnston, B.T. O'Flaherty and P.A. Jacobs

config.title = "Flying projectile."
config.dimensions = 2
config.axisymmetric = true

-- Read in some config parameters as global variables
dofile('sim-config.lua')

-- Gas model and flow state setup
setGasModel("ideal-air-gas-model.lua")
pInit = 1.0e5
TInit = 348.43
initialLeft = FlowState:new{p=pInit, T=TInit}
initialRight = FlowState:new{p=1.0e-10, T=TInit}

-- Geometry, grid and block setup.
nxcells = 100
nycells = 2

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=L1, y=0.0}
c = Vector3:new{x=L1, y=H}
d = Vector3:new{x=0.0, y=H}

e = Vector3:new{x=L2, y=0.0}
f = Vector3:new{x=L3, y=0.0}
g = Vector3:new{x=L3, y=H}
h = Vector3:new{x=L2, y=H}

ab = Line:new{p0=a, p1=b}
bc = Line:new{p0=b, p1=c}
ad = Line:new{p0=a, p1=d}
dc = Line:new{p0=d, p1=c}

ef = Line:new{p0=e, p1=f}
fg = Line:new{p0=f, p1=g}
eh = Line:new{p0=e, p1=h}
hg = Line:new{p0=h, p1=g}

grid0 = StructuredGrid:new{psurface=makePatch{north=dc, east=bc, south=ab, west=ad},
        niv=nxcells+1, njv=nycells+1}
grid1 = StructuredGrid:new{psurface=makePatch{north=hg, east=fg, south=ef, west=eh},
        niv=nxcells+1, njv=nycells+1}
-- Now we create our blocks. We are using the WallBC_WithSlip1 for the
-- moving mesh. This should account for the contraction/expansion of flow domains
blk0 = FBArray:new{grid=grid0, initialState=initialLeft, nib=1, njb=1, bcList={north=WallBC_WithSlip1:new{}, east=WallBC_WithSlip1:new{group='pistonUpstream'}, south=WallBC_WithSlip1:new{}, west=WallBC_WithSlip1:new{}}}
blk1 = FBArray:new{grid=grid1, initialState=initialRight, nib=1, njb=1, bcList={north=WallBC_WithSlip1:new{}, east=WallBC_WithSlip1:new{}, south= WallBC_WithSlip1:new{}, west=WallBC_WithSlip1:new{group='pistonDownStream'}}}
-- Set up the run-time loads for computing the piston movement
run_time_loads={ {group="pistonUpstream",
                              moment_centre=Vector3:new{x= 0., y= 0.}},
                 {group="pistonDownstream",
                              moment_centre=Vector3:new{x= 0., y= 0.}}
                }
-- Setting up our simulation
config.interpolation_order = 2
config.dt_init = 1.0e-7
config.fixed_time_step = false
config.max_time = 4.0e-2
config.max_step = 30000
config.dt_plot = config.max_time/10.0
config.cfl_value = 0.1
config.flux_calculator = "ausmdv"
-- The following lines are the important selections related
-- to a moving grid simulation. We need to:
-- 1. Select a gasdynamic update scheme compatible with
--    a moving grid simulation. The available choices are:
--    + "moving_grid_1_stage"
--    + "moving_grid_2_stage"
-- 2. Select a mode of grid motion. The available choices are:
--    + "none" : no grid motion.
--    + "user_defined" : user controls vertex velocities with Lua script
--    + "shock_fitting" : to be used in the special case of a shock-fitting
--                        simulation
-- 3. Since we are doing "user_defined" grid motion, we need to supply
--    the name of the Lua script that contains the motion definition.
config.gasdynamic_update_scheme = "moving_grid_1_stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "grid-motion.lua"
-- We require a udf in which we calculate the projectile dynamics and also
-- do some data manipulation, storage and writing
config.udf_supervisor_file='udf-process.lua'
-- We are using the built-in computation of loads on surfaces to calculate the
-- forces on the projectile. Here we specify that the calculation should be
-- done, and the step frequency with which to compute the forces
config.compute_run_time_loads = true
config.run_time_loads_count = 1
-- Set up the userPad for storing our piston data, i.e. position and velocity
config.user_pad_length = 2
-- The user pad values are automatically populated with 0.0, however, here it
-- is done explicitly to demonstrate how it is done
user_pad_data = {0.0, 0.0}