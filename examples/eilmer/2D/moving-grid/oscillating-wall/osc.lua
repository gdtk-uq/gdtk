-- osc.lua
-- Example to exercise the setting of vertex velocities along a boundary.
-- We make a flow domain where part of the wall (b0-c0) deflects periodically.
--
-- Peter J. 2022-02-26
--
config.title = "Oscillating wall"
config.dimensions = 2
setGasModel("ideal-air-gas-model.lua")
initialGas = FlowState:new{p=100.0e3, T=300.0, velx=700.0}

-- Geometry, grid and block setup.
--  a1----b1----c1---------d1
--  |     |     |          |
--  | [0] | [1] |    [2]   |
--  |     |     |          |
--  a0----b0----c0---------d0
L1 = math.pi
H = math.pi
a0 = Vector3:new{x=-L1, y=0}; a1 = Vector3:new{x=-L1, y=H}
b0 = Vector3:new{x=0, y=0}; b1 = Vector3:new{x=0, y=H}
c0 = Vector3:new{x=L1, y=0}; c1 = Vector3:new{x=L1, y=H}
d0 = Vector3:new{x=3*L1, y=0}; d1 = Vector3:new{x=3*L1, y=H}

blk0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
blk1 = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
blk2 = CoonsPatch:new{p00=c0, p10=d0, p11=d1, p01=c1}

niv = 31; njv = 31
grid0 = StructuredGrid:new{psurface=blk0, niv=niv, njv=njv}
grid1 = StructuredGrid:new{psurface=blk1, niv=niv, njv=njv}
grid2 = StructuredGrid:new{psurface=blk2, niv=2*niv, njv=njv}

blk0 = FluidBlock:new{
   grid=grid0, initialState=initialGas,
   bcList={west=InFlowBC_Supersonic:new{flowState=initialGas}}
}
blk1 = FluidBlock:new{
   grid=grid1, initialState=initialGas,
}
blk2 = FluidBlock:new{
   grid=grid2, initialState=initialGas,
   bcList={east=OutFlowBC_Simple:new{}}
}
identifyBlockConnections()

config.max_time = 80.0e-3
config.max_step = 3000
config.dt_plot = config.max_time/40

config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "user_defined"
config.udf_grid_motion_file = "grid-motion.lua"
