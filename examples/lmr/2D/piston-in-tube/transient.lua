-- transient.lua
print("Piston in tube, 2 FluidBlocks.")
--
-- 0. Assume that a previous processing has step set up the grids.
config.solver_mode = "transient"
config.dimensions = 2
config.axisymmetric = true
--
-- 1. Gas model and flow state setup
setGasModel('ideal-air.gas')
pInit = 100.0e3 -- Pa
TInit = 348.43 -- K
-- density is 1.0 kg/m**3
initialLeft = FlowState:new{p=pInit, T=TInit}
initialRight = FlowState:new{p=1.0e-10, T=TInit}
--
-- 2. Fluid blocks, with initial flow states and boundary conditions.
-- Block boundaries that are not otherwise assigned a boundary condition
-- are initialized as WallBC_WithSlip.
flowDict = {
   initialLeft=initialLeft,
   initialRight=initialRight
}
bcDict = {
   piston_upstream_face=WallBC_WithSlip1:new{group="pistonUpstream"},
   piston_downstream_face=WallBC_WithSlip1:new{group="pistonDownstream"}
}
makeFluidBlocks(bcDict, flowDict)
--
-- 3. Simulation parameters.
config.max_time = 40.0e-3
config.max_step = 3000
config.dt_plot = config.max_time/8

config.gasdynamic_update_scheme = "moving_grid_1_stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "grid-motion.lua"

-- Calculate the projectile dynamics in user-defined functions
-- but using loads computed by the flow solver.
config.udf_supervisor_file='udf-supervisor.lua'
config.compute_run_time_loads = true
config.run_time_loads_count = 1
run_time_loads={
   {group="pistonUpstream", moment_centre=Vector3:new{x=0, y=0}},
   {group="pistonDownstream", moment_centre=Vector3:new{x=0, y=0}}
}
-- Dimension userPad for storing piston position and velocity.
config.user_pad_length = 2
user_pad_data = {0, 0}
