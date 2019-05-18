-- pit2.lua
-- Authors: Rowan G., FabZ, Ingo J., Peter J. 
-- Date: 2017-01-04 -- 2019-05-18
--
-- Dimensions match verification case 2 in the 16th AFMC paper (2007)
-- "Development of Casbar: a Two-phase Flow Code for the Interior Ballistics
-- Problem", R.J Gollan, I.A. Johnston, B.T. O'Flaherty and P.A. Jacobs

config.title = "Piston in tube, 2 FluidBlocks."
config.dimensions = 2
config.axisymmetric = true
dofile('sim-config.lua')

-- Gas model and initial gas states.
setGasModel("ideal-air-gas-model.lua")
pInit = 100.0e3 -- Pa
TInit = 348.43 -- K
-- density is 1.0 kg/m**3
initialLeft = FlowState:new{p=pInit, T=TInit}
initialRight = FlowState:new{p=1.0e-10, T=TInit}

-- Geometry, grid and block setup.
-- Gas region that drives piston.
driver_patch = CoonsPatch:new{p00=Vector3:new{x=0, y=0},
                              p10=Vector3:new{x=L1, y=0},
                              p11=Vector3:new{x=L1, y=H},
                              p01=Vector3:new{x=0, y=H}}
-- Gas region that is compressed by piston.
driven_patch = CoonsPatch:new{p00=Vector3:new{x=L2, y=0.0},
                              p10=Vector3:new{x=L3, y=0.0},
                              p11=Vector3:new{x=L3, y=H},
                              p01=Vector3:new{x=L2, y=H}}
nxc = 100; nyc = 2
grid0 = StructuredGrid:new{psurface=driver_patch, niv=nxc+1, njv=nyc+1}
grid1 = StructuredGrid:new{psurface=driven_patch, niv=nxc+1, njv=nyc+1}

blk0 = FluidBlock:new{
   grid=grid0, initialState=initialLeft,
   bcList={east=WallBC_WithSlip1:new{group='pistonUpstream'}}
}
blk1 = FluidBlock:new{
   grid=grid1, initialState=initialRight,
   bcList={west=WallBC_WithSlip1:new{group='pistonDownStream'}}
}

config.interpolation_order = 2
config.flux_calculator = "ausmdv"
config.cfl_value = 0.1
config.dt_init = 1.0e-7
config.max_time = 40.0e-3
config.max_step = 30000
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
