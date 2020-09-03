-- spike.lua
-- Optimize the profile of an aerospike nozzle for ENGG7601 assignment.
--
-- Modelled on the Aerospike nozzles described in
-- G.V.R. Rao "Spike nozzle contour for optimal thrust."
-- Planetary and Space Science vol.4, Jan 1961, pages 92-101.
--
-- Peter J. 2019-09-03
--
config.title = "Flow over a toroidal aerospike nozzle."
config.axisymmetric = true
print(config.title)

-- Gas model and flow states for simulation.
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)

-- The stagnation gas represents a reservoir condition inside the rocket motor.
-- The low_pressure_gas is an arbitrary fill condition for the blocks.
-- It will be swept away.
-- The external_stream will provide an environment for the rocket's exhaust gas.
stagnation_gas = FlowState:new{p=5.0e6, T=3000.0}
low_pressure_gas = FlowState:new{p=300.0, T=300.0}
external_stream = FlowState:new{p=80.0e3, T=300.0, velx=2.0e3}

-- Define the geometry of the throat region and aerospike nozzle surface.
-- Use some of the notation from Rao's paper.
-- Remember that Eilmer length units are metres.
cm = 1.0e-2 -- centimetres per metre
yE = 10*cm -- outer radial position of throat
L = 10*cm -- length of aerospike
L2 = 5*cm -- downstream extension length for flow domain
H = 15*cm -- height/outer-radial-position of flow domain
yD = 2*cm -- radial position of D
alpha1 = math.rad(45) -- angle to control-point 1
alpha2 = math.rad(10) -- angle to control point 2
sqrt2 = math.sqrt(2)
pnts = {
   -- Points defining the block with the nozzle surface.
   E = Vector3:new{x=0, y=yE},
   T = Vector3:new{x=-sqrt2*cm, y=yE-sqrt2*cm},
   D = Vector3:new{x=L, y=yD},
   TDb1 = Vector3:new{x=-sqrt2*cm+0.4*L*math.cos(alpha1),
                      y=yE-sqrt2*cm-0.4*L*math.sin(alpha1)},
   TDb2 = Vector3:new{x=L-0.4*L*math.cos(alpha2),
                      y=yD+0.4*L*math.sin(alpha2)},
   F = Vector3:new{x=L, y=7*cm},
   EFb1 = Vector3:new{x=3.0*cm, y=8*cm},
   EFb2 = Vector3:new{x=6.5*cm, y=7*cm},
   -- Points defining the combustion chamber and throat.
   b0 = Vector3:new{x=-8*cm, y=8*cm},
   b1 = Vector3:new{x=-4*cm, y=8*cm},
   b2 = Vector3:new{x=-3*cm, y=10*cm},
   c0 = Vector3:new{x=-8*cm, y=12*cm},
   c1 = Vector3:new{x=-4*cm, y=12*cm},
   c2 = Vector3:new{x=-2*cm, y=12*cm},
   -- Points defining the exterior surface of the combustor.
   d0 = Vector3:new{x=-8*cm, y=12.2*cm},
   d1 = Vector3:new{x=-4*cm, y=12.2*cm},
   d2 = Vector3:new{x=-2*cm, y=12.2*cm},
   -- Points around the exterior boundary of the flow domain.
   e0 = Vector3:new{x=-8*cm, y=H},
   e1 = Vector3:new{x=2*cm, y=H},
   e2 = Vector3:new{x=L, y=H},
   e3 = Vector3:new{x=L+L2, y=H},
   e4 = Vector3:new{x=L+L2, y=7*cm},
   e5 = Vector3:new{x=L+L2, y=yD},
   Ee1b1 = Vector3:new{x=2*cm, y=12*cm},
   Ee1b2 = Vector3:new{x=2*cm, y=13*cm}
}

lines = {}
patches = {}
-- Block leading to the throat.
lines.w0 = Line:new{p0=pnts.b0, p1=pnts.c0}
lines.s0 = Bezier:new{points={pnts.b0, pnts.b1, pnts.b2, pnts.T}}
lines.n0 = Bezier:new{points={pnts.c0, pnts.c1, pnts.c2, pnts.E}}
lines.TE = Line:new{p0=pnts.T, p1=pnts.E}
patches[0] = CoonsPatch:new{north=lines.n0, south=lines.s0, west=lines.w0, east=lines.TE}
-- Block adjacent to the nozzle surface.
lines.TD = Bezier:new{points={pnts.T, pnts.TDb1, pnts.TDb2, pnts.D}}
lines.EF = Bezier:new{points={pnts.E, pnts.EFb1, pnts.EFb2, pnts.F}}
lines.DF = Line:new{p0=pnts.D, p1=pnts.F}
patches[1] = CoonsPatch:new{north=lines.EF, south=lines.TD, west=lines.TE, east=lines.DF}
-- Block downstream of nozzle.
patches[2] = CoonsPatch:new{p00=pnts.D, p10=pnts.e5, p11=pnts.e4, p01=pnts.F}
-- Block outside combustor.
lines.n3 = Line:new{p0=pnts.e0, p1=pnts.e1}
lines.s3 = Bezier:new{points={pnts.d0, pnts.d1, pnts.d2, pnts.E}}
lines.w3 = Line:new{p0=pnts.d0, p1=pnts.e0}
lines.e3 = Bezier:new{points={pnts.E, pnts.Ee1b1, pnts.Ee1b2, pnts.e1}}
patches[3] = CoonsPatch:new{north=lines.n3, south=lines.s3, west=lines.w3, east=lines.e3}
-- Continue downstream, outer part of flow domain.
lines.n4 = Line:new{p0=pnts.e1, p1=pnts.e2}
lines.e4 = Line:new{p0=pnts.F, p1=pnts.e2}
patches[4] = CoonsPatch:new{north=lines.n4, south=lines.EF, west=lines.e3, east=lines.e4}
patches[5] = CoonsPatch:new{p00=pnts.F, p10=pnts.e4, p11=pnts.e3, p01=pnts.e2}

-- Mesh the patches, with particular discretisation.
ny = 30
cleft = RobertsFunction:new{end0=true, end1=false, beta=1.2}
cright = RobertsFunction:new{end0=false, end1=true, beta=1.2}
grids = {}
grids[0] = StructuredGrid:new{psurface=patches[0],
                              cfList={south=cright, north=cright},
                              niv=61, njv=ny+1}
grids[1] = StructuredGrid:new{psurface=patches[1],
                              cfList={north=cleft, south=cleft},
                              niv=81, njv=ny+1}
grids[2] = StructuredGrid:new{psurface=patches[2],
                              niv=31, njv=ny+1}
grids[3] = StructuredGrid:new{psurface=patches[3],
                              cfList={east=cleft}, niv=61, njv=ny+1}
grids[4] = StructuredGrid:new{psurface=patches[4],
                              cfList={west=cleft, east=cleft, south=cleft},
                              niv=81, njv=ny+1}
grids[5] = StructuredGrid:new{psurface=patches[5],
                              cfList={west=cleft, east=cleft},
                              niv=31, njv=ny+1}

-- Define the flow-solution blocks.
blks = {}
for ib=0,5 do
   blks[ib] = FluidBlock:new{grid=grids[ib], initialState=low_pressure_gas}
end

-- Set boundary conditions.
identifyBlockConnections()
-- then, directly specify the stagnation conditions for the subsonic inflow.
blks[0].bcList['west'] = InFlowBC_FromStagnation:new{stagCondition=stagnation_gas}
-- to get loads on thrust surface, add that boundary condition to the group
blks[1].bcList['south'] = WallBC_WithSlip:new{group="loads"}
blks[2].bcList['east'] = OutFlowBC_Simple:new{}
blks[5].bcList['east'] = OutFlowBC_Simple:new{}
blks[3].bcList['west'] = InFlowBC_Supersonic:new{flowCondition=external_stream}

-- Some more simulation configuration.
config.max_time = 1.0e-3  -- seconds
config.max_step = 30000
config.dt_init = 1.0e-7
config.dt_plot = config.max_time/10.0
config.dt_loads = config.max_time/10.0
config.write_loads = true
-- We don't want to see the noise from starting the sharp corner.
config.max_invalid_cells = 10
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
-- Also need to deal with the difficult hypersonic patch
-- on the outside of the combustor, just upstream of the throat.
config.suppress_reconstruction_at_shocks = true

dofile("sketch-domain.lua")
