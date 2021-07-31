-- cone.lua
-- Conical/source flow example.
-- PJ, 2021-07-31

config.title = "Conical/source flow."
print(config.title)
config.dimensions = 2
config.axisymmetric = true

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=5955.0, T=304.0, velx=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, velx=1000.0}

print("inflow=", inflow)
print("T=", inflow.T, "density=", inflow.rho, "sound speed= ", inflow.a)
print("inflow Mach number=", 1000.0/inflow.a)

r_a = 1.0; r_b = 1.5; r_c = 2.0
theta = 7.0*math.pi/180.0; snth = math.sin(theta); csth = math.cos(theta)
a0 = {x=r_a, y=0.0}
b0 = {x=r_b, y=0.0}; b1 = {x=r_b*csth, y=r_b*snth}
c0 = {x=r_c, y=0.0}; c1 = {x=r_c*csth, y=r_c*snth}
origin = {x=0.0, y=0.0}

if false then
   -- Circular arc inflow boundary
   a1 = {x=r_a*csth, y=r_a*snth}
   a0a1 = Arc:new{p0=a0, p1=a1, centre=origin}
else
   -- Straight-line (vertical) inflow boundary
   a1 = {x=r_a, y=r_a*math.tan(theta)}
   a0a1 = Line:new{p0=a0, p1=a1}
end
a1b1 = Line:new{p0=a1, p1=b1} -- upper boundary, conical surface
a0b0 = Line:new{p0=a0, p1=b0} -- lower boundary, axis
b0c0 = Line:new{p0=b0, p1=c0}
b1c1 = Line:new{p0=b1, p1=c1}
b0b1 = Arc:new{p0=b0, p1=b1, centre=origin}
c0c1 = Arc:new{p0=c0, p1=c1, centre=origin} -- exit

quad0 = CoonsPatch:new{north=a1b1, east=b0b1, south=a0b0, west=a0a1}
quad1 = CoonsPatch:new{north=b1c1, east=c0c1, south=b0c0, west=b0b1}
grid0 = StructuredGrid:new{psurface=quad0, niv=21, njv=11}
grid1 = StructuredGrid:new{psurface=quad1, niv=21, njv=11}
blk0 = FluidBlock:new{grid=grid0, initialState=initial}
blk1 = FluidBlock:new{grid=grid1, initialState=initial}

-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList['west'] = InFlowBC_ConstFlux:new{flowState=inflow, x0=0.0, y0=0.0, r=r_a}
blk1.bcList['east'] = OutFlowBC_Simple:new{}

config.max_time = 15.0e-3
config.max_step = 3000
config.dt_plot = 5.0e-3
config.dt_init = 1.0e-6
config.extrema_clipping = false
