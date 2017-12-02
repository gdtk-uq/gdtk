-- conepe.lua
-- Parametric, extended setup for sharp-cone simulation.
-- PJ & RG
-- 2016-09-23 -- adapted from conep.lua

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.axisymmetric = true

-- The gas model is defined via a gas-model file.
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0}
-- Compute inflow from Mach number.
inflow_gas = FlowState:new{p=95.84e3, T=1103.0}
M = 1.5
Vx = M * inflow_gas.a
print("inflow velocity Vx=", Vx)
print("dynamic pressure q=", 1/2*inflow_gas.rho*Vx*Vx)
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=Vx}
print("T=", inflow.T, "density=", inflow.rho, "sound speed= ", inflow.a)

-- Parameters defining cone and flow domain.
theta = 32 -- cone half-angle, degrees
L = 0.8    -- axial length of cone, metres
rbase = L * math.tan(math.pi*theta/180.0)
x0 = 0.2   -- upstream distance to cone tip
H = 2.0    -- height of flow domain, metres
config.title = string.format("Mach %.1f flow over a %.1f-degree cone.",
			     M, theta)
print(config.title)

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=x0, y=0.0}
c = Vector3:new{x=x0+L, y=rbase}
d = Vector3:new{x=x0+L, y=H}
e = Vector3:new{x=x0, y=H}
f = Vector3:new{x=0.0, y=H}
ab = Line:new{p0=a, p1=b} -- lower boundary, axis
bc = Line:new{p0=b, p1=c} -- lower boundary, cone surface
fe = Line:new{p0=f, p1=e}; ed = Line:new{p0=e, p1=d} -- upper boundary
af = Line:new{p0=a, p1=f} -- vertical line, inflow
be = Line:new{p0=b, p1=e} -- vertical line, between quads
cd = Line:new{p0=c, p1=d} -- vertical line, outflow
quad0 = makePatch{north=fe, east=be, south=ab, west=af}
quad1 = makePatch{north=ed, east=cd, south=bc, west=be, gridType="ao"}
-- extend the flow domain
xend = x0 + 2*L
quad2 = CoonsPatch:new{p00=c, p10=Vector3:new{x=xend, y=rbase/2},
		       p11=Vector3:new{x=xend, y=H}, p01=d}
-- Mesh the patches, with particular discretisation.
dx = 1.0/40
nx0 = math.floor(x0/dx); nx1 = math.floor(L/dx); ny = math.floor(H/dx)
grid0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1}
grid2 = StructuredGrid:new{psurface=quad2, niv=nx1+1, njv=ny+1}
-- Define the flow-solution blocks.
blk0 = FluidBlock:new{grid=grid0, initialState=inflow}
blk1 = FluidBlock:new{grid=grid1, initialState=initial}
blk2 = FluidBlock:new{grid=grid2, initialState=initial}
-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = InFlowBC_Supersonic:new{flowState=inflow}
blk2.bcList[east] = OutFlowBC_Simple:new{}

-- add history point 1/3 along length of cone surface
setHistoryPoint{x=2*b.x/3+c.x/3, y=2*b.y/3+c.y/3}
-- add history point 2/3 along length of cone surface
setHistoryPoint{ib=1, i=math.floor(2*nx1/3), j=0}

-- Do a little more setting of global data.
config.max_time = 30.0e-3  -- seconds
config.max_step = 15000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_history = 10.0e-5

dofile("sketch-domain-extended.lua")
