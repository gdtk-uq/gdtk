-- corner.lua

config.title = "Mach 5.09 shock wave moving past a backward facing step"
print(config.title)
config.dimensions = 2

-- select flux calculator for comparison

config.flux_calculator = "ausmdv"
--config.flux_calculator = "ausm_plus_up"
--config.flux_calculator = "efm"
--config.flux_calculator = "hlle"
--config.flux_calculator = "adaptive"

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=1000.0, T=298.0, velx=0.0, vely=0.0}
-- the shock wave, not the flow, is at Mach 5.09
inflow = FlowState:new{p=30059.45, T=1781.08, velx=1411.09, vely=0.0} 

-- points
ii = 1.0e-1; jj = 1.0e-1
a = Vector3:new{x = 0, y = 32*jj}
b = Vector3:new{x = 20*ii, y = 32*jj}
c = Vector3:new{x = 20*ii, y = 0}
d = Vector3:new{x = 70*ii, y = 0}
e = Vector3:new{x = 70*ii, y = 32*jj}
f = Vector3:new{x = 70*ii, y = 70*jj}
g = Vector3:new{x = 20*ii, y = 70*jj}
h = Vector3:new{x = 0, y = 70*jj}

-- lines
ab = Line:new{p0 = a, p1 = b}
bg = Line:new{p0 = b, p1 = g}
ah = Line:new{p0 = a, p1 = h}
hg = Line:new{p0 = h, p1 = g}
cd = Line:new{p0 = c, p1 = d}
de = Line:new{p0 = d, p1 = e}
be = Line:new{p0 = b, p1 = e}
cb = Line:new{p0 = c, p1 = b}
ef = Line:new{p0 = e, p1 = f}
gf = Line:new{p0 = g, p1 = f}

-- Mesh the patches, with particular discretisation.
dx = 1; dy = 1
ny0 = 38*dy; ny2 = 32*dy
nx0 = 20*dx; nx1 = 50*dx; nx3 = 10*dx
quad0 = makePatch{north=hg, east=bg, south=ab, west=ah}
quad1 = makePatch{north=gf, east=ef, south=be, west=bg}
quad2 = makePatch{north=be, east=de, south=cd, west=cb}

grid0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny0+1}
grid1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny0+1}
grid2 = StructuredGrid:new{psurface=quad2, niv=nx1+1, njv=ny2+1}

-- Define the flow-solution blocks.
blk0 = FluidBlock:new{grid=grid0, initialState=initial}
blk1 = FluidBlock:new{grid=grid1, initialState=initial}
blk2 = FluidBlock:new{grid=grid2, initialState=initial}

-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_Supersonic:new{flowState=inflow}
blk1.bcList[east] = OutFlowBC_Simple:new{}
blk1.bcList[north] = OutFlowBC_Simple:new{}
blk2.bcList[east] = OutFlowBC_Simple:new{}
blk2.bcList[south] = OutFlowBC_Simple:new{}

config.max_time = 10.0e-3  -- seconds
config.max_step = 1000
config.dt_init = 1.0e-7
config.cfl_value = 0.5
config.dt_plot = 10.0e-5
config.dt_history = 1.0e-5

dofile("sketch-domain.lua")
