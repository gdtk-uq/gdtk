-- =====
-- grid-setup:start
--
-- G1. Geometry
a0 = {x=0.0, y=0.0};     a1 = {x=0.0, y=1.0}
b0 = {x=0.2, y=0.0};     b1 = {x=0.2, y=1.0}
c0 = {x=1.0, y=0.29118}; c1 = {x=1.0, y=1.0}
--
quad0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
quad1 = AOPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
--
-- G2. Grids
nx0 = 10; nx1 = 30; ny = 40
print("registering grid 0")
grid0 = registerFluidGrid{
   grid=StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1},
   fsTag="inflow",
   bcTags={west="inflow"}
}
grid1 = registerFluidGrid{
   grid=StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1},
   fsTag="initial",
   bcTags={east="outflow"}
}
identifyGridConnections()
-- grid-setup:end
-- =====
-- sim-setup:start
--
-- S0. Assume that a previous processing has step set up the grids.
--
-- S1. Domain type, gas model and flow states
config.solver_mode = "transient"
config.axisymmetric = true
setGasModel('ideal-air.gas')
initial = FlowState:new{p=5955.0, T=304.0} -- Pa, degrees K
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
flowDict = {initial=initial, inflow=inflow}
--
-- S2. Fluid blocks, with initial flow states and boundary conditions.
-- Block boundaries that are not otherwise assigned a boundary condition
-- are initialized as WallBC_WithSlip.
bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
--
makeFluidBlocks(bcDict, flowDict)
--
--
-- S3. Simulation parameters.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_plot = 1.5e-3
config.extrema_clipping = false
--
-- S4. Add some history recording at 1/3 and 2/3 length along cone
setHistoryPoint{x=b0.x + (1/3)*(c0.x-b0.x), y=c0.y/3}
setHistoryPoint{ib=1, i=math.floor(2*nx1/3), j=0}
config.dt_history = 0.1e-3
-- sim-setup:end
