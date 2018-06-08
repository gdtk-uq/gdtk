-- Author: Ingo HJ Jahn
-- Date: 2018-06-07
--
-- This Lua script is used to prepare an Eilmer4
-- simulation with two-way coupling between the 
-- the fluid and the projectile that is accelerated
-- by the fluid pressure. This simulation makes use 
-- of the user defined supervisory file and user 
-- moving grid capabilities.
--


config.title = "Projectile being propelled by pressurised reservoir"
config.dimensions = 2

-- Read in some config parameters as global variables
dofile('sim-config.lua')

maxTime = 50e-3

config.udf_supervisor_file='udf-process.lua'   -- manages calculation of force acting on projectile and acceleration 
config.interpolation_order = 2
config.dt_init = 1.0e-8
config.fixed_time_step = false
config.max_time = maxTime
config.max_step = 3000000
config.dt_plot = maxTime/100.0
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
config.udf_grid_motion_file = "move-grid.lua"

config.viscous = false
--config.suppress_reconstruction_at_boundaries = true

-- Gas model and flow state setup
setGasModel("ideal-air-gas-model.lua")
initial = FlowState:new{p=startP, T=startT}

-- Geometry, grid and block setup.
nxcells = 50
nycells = 3

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=startL, y=0.0}
c = Vector3:new{x=0.0, y=H}
d = Vector3:new{x=startL, y=H}

ab = Line:new{p0=a, p1=b}
ac = Line:new{p0=a, p1=c}
bd = Line:new{p0=b, p1=d}
cd = Line:new{p0=c, p1=d}

grid = StructuredGrid:new{psurface=makePatch{north=cd, east=bd, south=ab, west=ac},
			  niv=nxcells+1, njv=nycells+1} 

blk = FluidBlock:new{grid=grid, initialState=initial}

-- Boundary condition uses development version of WallBC to correctly include
-- work done by fluid as projectile is pushed.
blk.bcList[east] = WallBC_WithSlip2:new{}


