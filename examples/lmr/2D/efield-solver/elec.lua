-- elec.lua
-- Test file for electric field solver
-- NNG
-- 2021-01-20

config.solver_mode = 'transient'
config.dimensions = 2
config.solve_electric_field = true
config.electric_field_count=1
config.conductivity_model_name = "test"

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel('ideal-air.lua')
initial = FlowState:new{p=101.35e3, T=300.0, velx=0.0}

flowDict = {
   initial=initial
}

bcDict = {
   north=WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
   south=WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
   east=WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
   west=WallBC_WithSlip:new{field_bc=FixedField_Test:new{}},
}


makeFluidBlocks(bcDict, flowDict)

-- No boundary conditions needed, just yet
identifyBlockConnections()

-- global data.
config.max_step = 1
