-- cst.lua
-- High-performance shock tube with helium driving air
-- in a constant-diameter tube.  The temperatures in the air
-- are high enough to induce strong thermochemical effects.
--
-- Adapted/ported from 3D/sod-shock-tube and eilmer3 classic shock tube.
-- PJ 2019-07-15

config.title = "High-performance shock tube with helium driving air."
print(config.title)

nsp, nmodes = setGasModel('lut-plus-helium-gas-model.lua')
print("GasModel nsp= ", nsp, " nmodes= ", nmodes)
helium = FlowState:new{p=30.0e6, T=3000, massf={ideal=1.0}}
air = FlowState:new{p=30.0e3, T=300.0, massf={lut=1.0}}

L = 1.0; R = 0.1
a0 = {x=0.0, y=0.0}; a1 = {x=0.0, y=R}
b0 = {x=L/2, y=0.0}; b1 = {x=L/2, y=R}
c0 = {x=L, y=0.0}; c1 = {x=L, y=R}
quad_driver = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
quad_driven = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}

grid_driver = StructuredGrid:new{psurface=quad_driver, niv=401, njv=3}
grid_driven = StructuredGrid:new{psurface=quad_driven, niv=401, njv=3}

blk_driver = FluidBlock:new{grid=grid_driver, initialState=helium}
blk_driven = FluidBlock:new{grid=grid_driven, initialState=air}
identifyBlockConnections()

config.max_time = 100.0e-6  -- seconds
config.max_step = 8000

