-- het.lua
-- lmr(5) script for Hypervelocity Expansion Tube (HET) simulation
-- https://www.austin.caltech.edu/hypersonics/het.html
--
-- This script is used to emulate Condition 3 described in:
-- Tal Schwartz et.al. (2023)
-- Laser absorption sensor tergeting potassium for hypersonic velocity,
-- temperature and enthalpy measurements.
-- AIAA Journal Vol 61(8): 3287-3297
--
-- Peter J. 2025-02-09
--
print("Hypervslocity Expansion Tube, Condition 3")
config.solver_mode = 'transient'
config.dimensions = 2
config.axisymmetric = true
config.grid_format = 'rawbinary'
config.field_format = 'rawbinary'
dofile('./blkIds.lua')
--
-- gas model and initial conditions
-- Because the 2T model is giving grief, we'll work with a 1T model.
nsp, nmodes, gmodel = setGasModel("co2-he-1T.gas")
print("number_species = ",nsp)
print("number_modes = ",nmodes)
print("gas_model = ",gmodel)
config.reactions_file = 'co2-he-1T-chemistry.chem'
config.reacting = true
T_init = 289.85 -- 16.5 degrees C
molef_driver = {He=1.0}
molef_driven = {CO2=1.0}
molef_accel = {CO2=1.0}
driver_gas = FlowState:new{p=3.080e6, T=T_init, massf=gmodel:molef2massf(molef_driver)}
test_gas = FlowState:new{p=1.2e3, T=T_init, massf=gmodel:molef2massf(molef_driven)}
acceleration_gas = FlowState:new{p=15.7, T=T_init, massf=gmodel:molef2massf(molef_driven)}
--
flowDict = {
   driver_gas=driver_gas,
   test_gas=test_gas,
   acceleration_gas=acceleration_gas
}
bcDict = {
   outFlow=OutFlowBC_Simple:new{},
   coldWall=WallBC_NoSlip_FixedT:new{Twall=T_init},
   warmWall=WallBC_NoSlip_FixedT:new{Twall=T_init},
   upstream_face_of_diaphragm = ExchangeBC_FullFacePlusUDF:new{
      otherBlock=downstreamBlk, otherFace='west', fileName="diaphragm.lua"},
   downstream_face_of_diaphragm = ExchangeBC_FullFacePlusUDF:new{
      otherBlock=upstreamBlk, otherFace='east', fileName="diaphragm.lua"}
}
--
makeFluidBlocks(bcDict, flowDict)
--
-- simulation settings
-- config.viscous = true
config.viscous = false -- TEMPORARY
config.spatial_deriv_locn = 'vertices'
config.spatial_deriv_calc = 'divergence'
config.viscous_signal_factor = 0.1
config.max_time = 5.0e-3 -- s
config.max_step = 1000000
config.gasdynamic_update_scheme = "predictor-corrector"
config.dt_init = 1.0e-9
config.cfl_value = 0.5
config.flux_calculator = "ausmdv"
config.cfl_count = 3
config.suppress_radial_reconstruction_at_xaxis = true
config.max_invalid_cells = 100
config.adjust_invalid_cell_data = true
--
config.dt_plot = 0.5e-4 -- s (at 0.05 ms)
--
-- A place to record the state of the diaphragm.
-- 0 == diaphragm closed
-- 1 == diaphragm open (ruptured)
config.user_pad_length = 1
user_pad_data = {0}
-- We want the intermediate-tube right-most block that sets the rupture state
-- to be on the MPI master task 0.  Its user_pad_data is broadcast.
mpiDistributeBlocks{ntasks=14, dist="load-balance", preassign={[upstreamBlk]=0}}
-- The function that sets the diaphragm state is also a user-defined function.
config.udf_supervisor_file='supervisor.lua'
--
-- Configure history points
config.dt_history = 1.0e-6
setHistoryPoint{x=5.179, y=0.0} -- at secondary diaphragm
setHistoryPoint{x=9.139, y=0.0} -- near end of acceleration tube
