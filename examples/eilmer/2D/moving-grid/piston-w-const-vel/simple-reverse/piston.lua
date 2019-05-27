-- Author: Rowan J. Gollan & Peter J.
-- Date: 2017-01-04 -- 2019-05-21
--
-- This simulation mimics the evolution of the gas in
-- the right-half of the classic Sod shock tube problem,
-- and exercises the moving-grid capability of Eilmer4.
-- The piston has a constant velocity, high enough to drive
-- a shock wave through the gas in front of it.
--
config.title = "Fast-moving piston with constant velocity, reversed."
config.dimensions = 2
pSpeed = 293.5 -- m/s

-- Gas model and flow state setup
setGasModel("ideal-air-gas-model.lua")
pInit = 1.0e4; TInit = 278.8
initial = FlowState:new{p=pInit, T=TInit}

-- Geometry, grid and block setup.
L = 0.5; H = 0.1
-- Gas region that drives piston.
patch0 = CoonsPatch:new{p00=Vector3:new{x=0, y=0},
                        p10=Vector3:new{x=L, y=0},
                        p11=Vector3:new{x=L, y=H},
                        p01=Vector3:new{x=0, y=H}}
grid0 = StructuredGrid:new{psurface=patch0, niv=51, njv=3}
blk0 = FluidBlock:new{grid=grid0, initialState=initial,
                      bcList={east=WallBC_WithSlip1:new{}}
}

-- Simulation control parameters
config.gasdynamic_update_scheme = "moving_grid_1_stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "grid-motion.lua"

config.max_time = 0.6e-3
config.max_step = 3000
config.dt_plot = 0.1e-3
config.cfl_value = 0.1






