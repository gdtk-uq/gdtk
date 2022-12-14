# cone20.py
# PJ 2022-12-14 modelled off the Eilmer4 example
#
config.title = "Mach 1.5 flow over a 20 degree cone."
print(config.title)
#
init_gas_model("ideal-air-gas-model.lua")
initial = FlowState(p=5955.0, T=304.0)
inflow = FlowState(p=95.84e3, T=1103.0, velx=1000.0)
# Demo: Verify Mach number of inflow and compute dynamic pressure.
print("inflow=", inflow)
print("T=", inflow.gas.T, "density=", inflow.gas.rho, "sound speed= ", inflow.gas.a)
print("inflow Mach number=", 1000.0/inflow.gas.a)
print("dynamic pressure q=", 1/2*inflow.gas.rho*1.0e6)
#
#
# Set up two quadrilaterals in the (x,y)-plane by first defining
# the corner nodes, then the lines between those corners.
a = Vector3(0.0, 0.0)
b = Vector3(0.2, 0.0)
c = Vector3(1.0, 0.29118)
d = Vector3(1.0, 1.0)
e = Vector3(0.2, 1.0)
f = Vector3(0.0, y=1.0)
from gdtk.geom.path import Line
ab = Line(p0=a, p1=b) # lower boundary, axis
bc = Line(p0=b, p1=c) # lower boundary, cone surface
fe = Line(p0=f, p1=e); ed = Line(p0=e, p1=d) # upper boundary
af = Line(p0=a, p1=f) # vertical line, inflow
be = Line(p0=b, p1=e) # vertical line, between quads
cd = Line(p0=c, p1=d) # vertical line, outflow
quad0 = CoonsPatch(north=fe, east=be, south=ab, west=af)
quad1 = CoonsPatch(north=ed, east=cd, south=bc, west=be)
#
# Mesh the patches, with particular discretisation.
#
nx0 = 10; nx1 = 30; ny = 40
grd0 = StructuredGrid(psurf=quad0, niv=nx0+1, njv=ny+1)
grd1 = StructuredGrid(psurf=quad1, niv=nx1+1, njv=ny+1)
#
blk0 = FluidBlock(i=0, grid=grd0, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow),'iplus':ExchangeBC()})
blk1 = FluidBlock(i=1, grid=grd1, initialState=initial,
                  bcs={'iminus':ExchangeBC(),'iplus':OutflowBC()})
#
config.axisymmetric = True
config.flux_calc = FluxCalc.ausmdv_plus_hanel
config.max_time = 5.0e-3
config.max_step = 3000
config.plot_dt = 1.5e-3
