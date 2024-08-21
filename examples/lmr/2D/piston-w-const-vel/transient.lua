-- transient.lua
print("Fast-moving piston with constant velocity.")
--
-- 0. Assume that a previous processing has step set up the grids.
config.solver_mode = "transient"
config.dimensions = 2
pSpeed = 293.5 -- m/s
--
-- 1. Gas model and flow state setup
setGasModel('ideal-air.gas')
pInit = 1.0e4; TInit = 278.8
initial = FlowState:new{p=pInit, T=TInit}
--
-- 2. Fluid blocks, with initial flow states and boundary conditions.
-- Block boundaries that are not otherwise assigned a boundary condition
-- are initialized as WallBC_WithSlip.
flowDict = {initial=initial}
bcDict = {piston_face=WallBC_WithSlip1:new{}}
makeFluidBlocks(bcDict, flowDict)
--
-- 3. Simulation parameters.
config.gasdynamic_update_scheme = "moving_grid_1_stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "grid-motion.lua"

config.max_time = 0.6e-3
config.max_step = 3000
config.dt_plot = 0.1e-3
config.cfl_value = 0.1
