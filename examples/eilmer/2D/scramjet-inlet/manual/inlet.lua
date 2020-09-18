-- inlet.lua
-- Optimize the profile of an scramjet inlet for MECH7101 assignment.
--
-- Peter J. 2020-08-01
--
config.title = "Flow into a scramjet inlet."
print(config.title)

-- Gas model and flow states for simulation.
nsp, nmodes, gmodel = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp)

-- The low_pressure_gas is an arbitrary fill condition for the blocks.
-- It will be swept away.
-- The external_stream will provide an environment for the rocket's exhaust gas.
low_pressure_gas = FlowState:new{p=100.0, T=250.0}
speed = 4.0 * low_pressure_gas.a -- set speed from Mach number
print("speed=", speed)
external_stream = FlowState:new{p=10.0e3, T=250.0, velx=speed}

-- Define the geometry of the scramjet ramp and a bit of duct
-- following that ramp.
-- Remember that Eilmer length units are metres.
L0 = 1.00
L1 = 0.30
L2 = 0.20
h0 = 0.05 -- adjust height of duct
h1 = 0.05 -- adjust height of b on ramp
h2 = 0.15 -- adjust height of c on ramp
h3 = 0.30
theta = math.rad(10.0) -- adjust dip of cowl
r = 0.15
Htop = 1.5
pnts = {
   -- Points defining the compression-ramp surface.
   a = {x=0.0, y=0.0},
   b = {x=L0/3, y=h1},
   c = {x=2*L0/3, y=h2},
   d = {x=L0, y=h3},
   d2 = {x=L0+0.5*L1, y=h3},
   e = {x=L0+L1, y=h3},
   f = {x=L0+L1+L2, y=0.95*h3},
   g = {x=0.0, y=h3},
   h = {x=L0-r*math.cos(theta), y=h3+h0-r*math.sin(theta)},
   i = {x=L0, y=h3+h0},
   i2 = {x=L0+0.5*L1, y=h3+h0},
   j = {x=L0+L1, y=h3+h0},
   k = {x=L0+L1+L2, y=h3+h0},
   l = {x=0.0, y=Htop},
   m = {x=L0-r, y=Htop},
   m2 = {x=L0+L1, y=Htop},
   n = {x=L0+L1+L2, y=Htop}
}
ramp = Bezier:new{points={pnts.a, pnts.b, pnts.c, pnts.d}}
t_truncate = 0.87
pnts.d1 = ramp(t_truncate) -- We are going to end block 0 here

lines = {}
patches = {}
-- Block leading to the throat.
lines.w0 = Line:new{p0=pnts.a, p1=pnts.g}
lines.s0 = SubRangedPath:new{underlying_path=ramp, t0=0.0, t1=t_truncate}
lines.n0 = Line:new{p0=pnts.g, p1=pnts.h}
lines.e0 = Line:new{p0=pnts.d1, p1=pnts.h}
patches[0] = AOPatch:new{north=lines.n0, south=lines.s0,
                         west=lines.w0, east=lines.e0}
-- 2 Blocks inside the cowl surface.
lines.n1 = Polyline:new{segments={Bezier:new{points={pnts.h, pnts.i, pnts.i2}},
                                  Line:new{p0=pnts.i2, p1=pnts.j}}}
lines.e1 = Line:new{p0=pnts.e, p1=pnts.j}
lines.s1 = Polyline:new{segments={Bezier:new{points={pnts.d1, pnts.d, pnts.d2}},
                                  Line:new{p0=pnts.d2, p1=pnts.e}}}
lines.n2 = Line:new{p0=pnts.j, p1=pnts.k}
lines.s2 = Line:new{p0=pnts.e, p1=pnts.f}
patches[1] = AOPatch:new{north=lines.n1, south=lines.s1,
                         west=lines.e0, east=lines.e1}
patches[2] = CoonsPatch:new{p00=pnts.e, p10=pnts.f, p11=pnts.k, p01=pnts.j}
-- Block upstream outer flow.
patches[3] = CoonsPatch:new{p00=pnts.g, p10=pnts.h, p11=pnts.m, p01=pnts.l}
-- 2 Blocks downstream outer flow.
lines.n4 = Line:new{p0=pnts.m, p1=pnts.m2}
lines.w4 = Line:new{p0=pnts.h, p1=pnts.m}
lines.e4 = Line:new{p0=pnts.j, p1=pnts.m2}
patches[4] = CoonsPatch:new{north=lines.n4, south=lines.n1,
                            west=lines.w4, east=lines.e4}
patches[5] = CoonsPatch:new{p00=pnts.j, p10=pnts.k, p11=pnts.n, p01=pnts.m2}

-- Mesh the patches, with particular discretisation.
cfy = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grids = {}
grids[0] = StructuredGrid:new{psurface=patches[0],
                              niv=61, njv=31}
grids[1] = StructuredGrid:new{psurface=patches[1],
                              niv=45, njv=31}
grids[2] = StructuredGrid:new{psurface=patches[2],
                              niv=21, njv=31}
grids[3] = StructuredGrid:new{psurface=patches[3],
                              cfList={east=cfy, west=cfy},
                              niv=61, njv=61}
grids[4] = StructuredGrid:new{psurface=patches[4],
                              cfList={east=cfy, west=cfy},
                              niv=45, njv=61}
grids[5] = StructuredGrid:new{psurface=patches[5],
                              cfList={east=cfy, west=cfy},
                              niv=21, njv=61}

-- Define the flow-solution blocks.
blks = {}
for ib=0,5 do
   blks[ib] = FluidBlock:new{grid=grids[ib], initialState=low_pressure_gas}
end

-- Set boundary conditions, first, connections
connectBlocks(blks[0], east, blks[1], west)
connectBlocks(blks[0], north, blks[3], south)
connectBlocks(blks[1], east, blks[2], west)
connectBlocks(blks[3], east, blks[4], west)
connectBlocks(blks[4], east, blks[5], west)
-- then, directly specify the other boundary conditions.
blks[0].bcList[west] = InFlowBC_Supersonic:new{flowCondition=external_stream}
blks[3].bcList[west] = InFlowBC_Supersonic:new{flowCondition=external_stream}
blks[2].bcList[east] = OutFlowBC_Simple:new{}
blks[5].bcList[east] = OutFlowBC_Simple:new{}

-- Some more simulation configuration.
config.max_time = 5.0e-3  -- seconds
config.max_step = 30000
config.dt_init = 1.0e-7
config.dt_plot = config.max_time/10.0
config.dt_loads = config.max_time/10.0
-- For difficult cells.
config.cfl_value = 0.4
config.cfl_count = 3
config.max_invalid_cells = 20
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
config.suppress_reconstruction_at_shocks = true

dofile("sketch-domain.lua")
