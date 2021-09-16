config.title = "Mach 1.5 flow over a 20 degree cone."
print(config.title)
--
-- 1. Flow domain and grids, dimensions in metres.
config.axisymmetric = true
--
a0 = {x=0.0, y=0.0};     a1 = {x=0.0, y=1.0}
b0 = {x=0.2, y=0.0};     b1 = {x=0.2, y=1.0}
c0 = {x=1.0, y=0.29118}; c1 = {x=1.0, y=1.0}
--
quad0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
quad1 = AOPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
--
grid0 = StructuredGrid:new{psurface=quad0, niv=11, njv=41}
grid1 = StructuredGrid:new{psurface=quad1, niv=31, njv=41}
--
-- 2. Gas model and flow states.  SI units.
setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=5955.0, T=304.0} -- Pa, degrees K
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}
--
-- 3. Fluid blocks, with initial flow states and boundary conditions.
blk0 = FluidBlock:new{grid=grid0, initialState=inflow}
blk1 = FluidBlock:new{grid=grid1, initialState=initial}
identifyBlockConnections()
blk0.bcList['west'] = InFlowBC_Supersonic:new{flowState=inflow}
blk1.bcList['east'] = OutFlowBC_Simple:new{}
--
-- 4. Simulation parameters.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_plot = 1.5e-3
config.extrema_clipping = false
