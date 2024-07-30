-- lmr(5) script for Langley Expansion Tube simulation
--
-- let.lua
-- Author: Peter J. and Rowan G.
-- Date: 2024-07-29
--
-- This script is used to emulate Case 1 described in:
--
-- Jacobs, P.A. (1994)
-- Numerical Simulation of Transient Hypervelocity Flow in an Expansion Tube
-- Computers & Fluids, 23:1 pp. 77--101
--
print("Langley Expansion Tube, Case 1")
config.solver_mode = 'transient'
config.dimensions = 2
config.axisymmetric = true
config.grid_format = 'rawbinary'
config.field_format = 'rawbinary'
dofile('./geom.lua')
dofile('./blkIds.lua')
--
-- gas model and initial conditions
setGasModel("ideal-helium.lua")
driver_gas = FlowState:new{p=33e6, T=330}
test_gas = FlowState:new{p=3.45e3, T=300}
acceleration_gas = FlowState:new{p=16.0, T=300}
--
flowDict = {
   driver_gas=driver_gas,
   test_gas=test_gas,
   acceleration_gas=acceleration_gas
}
bcDict = {
   outFlow=OutFlowBC_Simple:new{},
   coldWall=WallBC_NoSlip_FixedT:new{Twall=300.0},
   warmWall=WallBC_NoSlip_FixedT:new{Twall=330.0},
   upstream_face_of_diaphragm = ExchangeBC_FullFacePlusUDF:new{
      otherBlock=downstreamBlk, otherFace='west', fileName="diaphragm.lua"},
   downstream_face_of_diaphragm = ExchangeBC_FullFacePlusUDF:new{
      otherBlock=upstreamBlk, otherFace='east', fileName="diaphragm.lua"}
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- simulation settings
config.viscous = true
config.spatial_deriv_locn = 'vertices'
config.spatial_deriv_calc = 'divergence'
config.viscous_signal_factor = 0.1
config.max_time = 6.0e-3 -- s
config.max_step = 1000000
config.gasdynamic_update_scheme = "classic-rk3"
config.dt_init = 1.0e-9
config.cfl_value = 0.60
config.flux_calculator = "ausmdv"
config.cfl_count = 3
config.suppress_radial_reconstruction_at_xaxis = true
config.max_invalid_cells = 100
config.adjust_invalid_cell_data = true
--
config.dt_plot = 0.5e-4 -- s (at 0.1 ms)
--
-- A place to record the state of the diaphragm.
-- 0 == diaphragm closed
-- 1 == diaphragm open (ruptured)
config.user_pad_length = 1
user_pad_data = {0}
-- We want the intermediate-tube right-most block that sets the rupture state
-- to be on the MPI master task 0.  Its user_pad_data is broadcast.
mpiDistributeBlocks{ntasks=7, dist="load-balance",
                    preassign={[upstreamBlk]=0}}
-- The function that sets the diaphragm state is also a user-defined function.
config.udf_supervisor_file='supervisor.lua'
--
-- Configure history points
config.dt_history = 1.0e-6
-- at secondary diaphragm
setHistoryPoint{x=C0.x, y=C0.y}
setHistoryPoint{x=C1.x, y=C1.y}
-- near end of acceleration tube for plotting Pitot pressures
accEnd = 24
setHistoryPoint{x=accEnd, y=0.0}
setHistoryPoint{x=accEnd, y=2.12e-3}
setHistoryPoint{x=accEnd, y=18.8-3}
setHistoryPoint{x=accEnd, y=34e-3}
setHistoryPoint{x=accEnd, y=52.5e-3}
setHistoryPoint{x=accEnd, y=71e-3}
setHistoryPoint{x=accEnd, y=76e-3}
