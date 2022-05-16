-- Author: Ingo HJ Jahn
-- Date: 2018-10-10
--
-- This Lua script is used to prepare an Eilmer4
-- simulation to simulate Stokes Second Flow
-- See:
--  Schlichting, H (1979), Boudary-Layer Theory, McGraw-Hill, New York
--


config.title = "Stokes Second Flow"
config.dimensions = 2

-- Read in some config parameters as global variables
dofile('sim-config.lua')

maxTime = 50e-3

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

config.viscous = true
--config.suppress_reconstruction_at_boundaries = true

-- Gas model and flow state setup
setGasModel("ideal-air-gas-model.lua")
initial = FlowState:new{p=startP, T=startT}

-- Geometry, grid and block setup.

--    d----N----c
--    |         |
--    |         |
--    W         E
--    |         |
--    |         |
--    a----S----b


nxcells = 1*N_refine
nycells = 5*N_refine


a = Vector3:new{x=-L_x/2., y=0.0}
b = Vector3:new{x= L_x/2., y=0.0}
c = Vector3:new{x= L_x/2., y=L_y}
d = Vector3:new{x=-L_x/2., y=L_y}

cf = RobertsFunction:new{end0=true,end1=false,beta=1.05}

grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=a, p10=b, p01=d, p11=c},
			  niv=nxcells+1, njv=nycells+1,
              cfList={east=cf, west=cf} }

blk = FluidBlock:new{grid=grid, initialState=initial}

-- Boundary conditions:
--blk.bcList['south'] = WallBC_NoSlip_FixedT:new{Twall=startT}
-- blk.bcList['south'] = WallBC_NoSlip_Adiabatic:new{}
blk.bcList['south'] = UserDefinedBC:new{fileName='southwall-bc.lua'}

--blk.bcList['south'] = WallBC_TranslatingSurface_FixedT:new{Twall=startT, v_trans={x=0,y=0,z=0}}
--blk.bcList['south'] = WallBC_TranslatingSurface_Adiabatic:new{v_trans={x=0,y=0,z=0}}

-- east and west are set as extraploation as the flow is invariant in the x-direction
blk.bcList['east'] = OutFlowBC_SimpleExtrapolate:new{xOrder=0}
blk.bcList['west'] = OutFlowBC_SimpleExtrapolate:new{xOrder=0}
