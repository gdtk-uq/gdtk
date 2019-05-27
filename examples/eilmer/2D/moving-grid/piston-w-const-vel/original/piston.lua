-- Author: Rowan J. Gollan
-- Date: 2017-01-04
--
-- This Lua script is used to prepare an Eilmer4
-- simulation of piston motion and, most importantly,
-- its compressive effect on a slug of gas. This
-- simulation makes use of the user-defined moving
-- grid and user-defined boundary condition capabilities
-- in Eilmer4.
--
-- In this simulation, the piston is driven with a constant
-- velocity. Its speed is such that it drives a shock wave
-- through the gas in front of the piston. This simulation
-- mimics the evolution of the gas in the right-half of the
-- classic Sod shock tube problem.
--

config.title = "Fast-moving piston in duct."
config.dimensions = 2

-- Read in some config parameters as global variables
dofile('sim-config.lua')

maxTime = 0.6e-3

config.interpolation_order = 2
config.dt_init = 1.0e-7
config.fixed_time_step = false
config.max_time = maxTime
config.max_step = 30000
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
config.udf_grid_motion_file = "grid-motion.lua"

-- Gas model and flow state setup
setGasModel("ideal-air-gas-model.lua")
pInit = 1.0e4
TInit = 278.8
initial = FlowState:new{p=pInit, T=TInit}

-- Geometry, grid and block setup.
nxcells = 50
nycells = 4

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=L, y=0.0}
c = Vector3:new{x=0.0, y=H}
d = Vector3:new{x=L, y=H}

ab = Line:new{p0=a, p1=b}
ac = Line:new{p0=a, p1=c}
bd = Line:new{p0=b, p1=d}
cd = Line:new{p0=c, p1=d}

grid = StructuredGrid:new{psurface=makePatch{north=cd, east=bd, south=ab, west=ac},
			  niv=nxcells+1, njv=nycells+1}
blk = FluidBlock:new{grid=grid, initialState=initial}
-- We now set a special boundary condition for the boundary that acts as the piston face (WEST face).
blk.bcList[west] = UserDefinedFluxBC:new{fileName='piston-bc.lua', funcName='convective_flux'}






