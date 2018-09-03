-- sbo.lua
-- Optimize the profile of a payload fairing for ENGG7601 assignment.
-- Peter J. 2018-09-03
--
config.title = "Supersonic flow over a sharp 2D-axisymmetric body"
config.axisymmetric = true
print(config.title)

-- Flight conditions.
T_inf = 300.0
p_inf = 5000.0
M_inf = 3.0

-- Gas model and flow states for simulation.
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=0.1*p_inf, T=T_inf}
a_inf = initial.a -- sound-speed same as for inflow
v_inf = M_inf * a_inf
print("a_inf=", a_inf, "v_inf=", v_inf, "m/s")
inflow = FlowState:new{p=p_inf, T=T_inf, velx=v_inf}
print("density=", inflow.rho, "kg/m^3")
print("dynamic_pressure=", 0.5*inflow.rho*v_inf*v_inf, "Pa")

-- Geometry of flow domain; leave these fixed.
L0 = 1.0 -- metres
L1 = 5.0
L2 = 1.0
Rbase = 1.5
H1 = 3*Rbase
H2 = 5*Rbase

-- The shape of the fairing is defined by the Bezier curve
-- through the points b0 through b4.
-- Of these, the end points must remain fixed but the interior
-- points b1, b2 and b3 may have their y-values adjusted.
-- Below, these values are set, scaled on the base radius
-- of the fairing.
-- The values 0.25, 0.50 and 0.75 will define the Bezier
-- as a straight line and the fairing will be a cone.
pnts = {
   a = Vector3:new{x=-L0, y=0.0},
   b0 = Vector3:new{x=0.0, y=0.0},
   b1 = Vector3:new{x=0.25*L1, y=0.80*Rbase}, -- y is adjustable
   b2 = Vector3:new{x=0.50*L1, y=1.40*Rbase}, -- y is adjustable
   b3 = Vector3:new{x=0.75*L1, y=1.20*Rbase}, -- y is adjustable
   b4 = Vector3:new{x=L1, y=Rbase},
   c = Vector3:new{x=L1+L2, y=0.5*Rbase},
   d = Vector3:new{x=-L0, y=H1},
   e = Vector3:new{x=0.0, y=H1},
   f = Vector3:new{x=0.5*L1, y=H2},
   g = Vector3:new{x=L1, y=H2},
   h = Vector3:new{x=L1+L2, y=H2}
}

lines = {}
patches = {}
-- inflow and outflow blocks are quadrilaterals in xy-plane
patches[0] = CoonsPatch:new{p00=pnts.a, p10=pnts.b0, p11=pnts.e, p01=pnts.d}
patches[2] = CoonsPatch:new{p00=pnts.b4, p10=pnts.c, p11=pnts.h, p01=pnts.g}
-- lower boundary of middle block is body surface
lines.body = Bezier:new{points={pnts.b0, pnts.b1, pnts.b2, pnts.b3, pnts.b4}} 
-- upper boundary of middle block is a low-order Bezier
lines.upper = Bezier:new{points={pnts.e, pnts.f, pnts.g}}
-- vertical lines to provide west and east edges
lines.b0e = Line:new{p0=pnts.b0, p1=pnts.e}
lines.b4g = Line:new{p0=pnts.b4, p1=pnts.g}
patches[1] = CoonsPatch:new{north=lines.upper, south=lines.body,
                            west=lines.b0e, east=lines.b4g}

-- Mesh the patches, with particular discretisation.
ny = 60
clustery = RobertsFunction:new{end0=true, end1=false, beta=1.2}
grids = {}
grids[0] = StructuredGrid:new{psurface=patches[0],
                              cfList={west=clustery, east=clustery},
                              niv=17, njv=ny+1}
grids[1] = StructuredGrid:new{psurface=patches[1],
                              cfList={west=clustery, east=clustery},
                              niv=81, njv=ny+1}
grids[2] = StructuredGrid:new{psurface=patches[2],
                              cfList={west=clustery, east=clustery},
                              niv=17, njv=ny+1}

-- Define the flow-solution blocks.
blks = {}
blks[0] = FluidBlock:new{grid=grids[0], initialState=inflow}
blks[1] = FluidBlock:new{grid=grids[1], initialState=initial}
blks[2] = FluidBlock:new{grid=grids[2], initialState=initial}

-- Set boundary conditions.
identifyBlockConnections()
inBC = InFlowBC_Supersonic:new{flowState=inflow}
blks[0].bcList[west] = inBC
blks[0].bcList[north] = inBC
blks[1].bcList[north] = inBC
blks[1].bcList[south] = WallBC_WithSlip:new{group="loads"}
blks[2].bcList[east] = OutFlowBC_Simple:new{}

-- Some more configuration.
config.max_time = 30.0e-3  -- seconds
config.max_step = 2500
config.dt_init = 1.0e-6
config.dt_plot = config.max_time/10.0
config.dt_loads = config.max_time/10.0

dofile("sketch-domain.lua")
