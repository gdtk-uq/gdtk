# csbli.py
#
# A simulation of the crossing-shocks+boundary-layer interaction
# for Anand V. and David Mee's Discovery Project. 
#
# Peter J. 2023-07-29
#
config.title = "Crossing shocks interacting with a boundary-layer over a short plate."
config.viscous = True
config.max_time = 1.5e-3
config.max_step = 300000
config.dt_init = 1.0e-8
add_cfl_value(0.0, 0.5)
add_dt_plot(0.0, 0.3e-3)
#
inflow = FlowState(p=3.5315e3, T=340.2, velx=2158.2)
#
mm = 1.0e-3
z1 = 50*mm
x1 = 50*mm
L = 120*mm
from math import radians, cos, sin
theta = radians(10)
x2 = x1 + L*cos(theta)
dz2 = L*sin(theta)
z2 = z1 - dz2
th2 = radians(30)
L2 = dz2/sin(th2)
x3 = x2 + L2*cos(th2)
x4 = x3 + 30*mm
H = 50*mm

vol0 = TFIVolume(p000=Vector3(0.0,0.0,z1), p100=Vector3(x1,0.0,z1),
                 p110=Vector3(x1,0.0,0.0), p010=Vector3(0.0,0.0,0.0),
                 p001=Vector3(0.0,H,z1), p101=Vector3(x1,H,z1),
                 p111=Vector3(x1,H,0.0), p011=Vector3(0.0,H,0.0))
vol1 = TFIVolume(p000=Vector3(x1,0.0,z1), p100=Vector3(x2,0.0,z2),
                 p110=Vector3(x2,0.0,0.0), p010=Vector3(x1,0.0,0.0),
                 p001=Vector3(x1,H,z1), p101=Vector3(x2,H,z2),
                 p111=Vector3(x2,H,0.0), p011=Vector3(x1,H,0.0))
vol2 = TFIVolume(p000=Vector3(x2,0.0,z2), p100=Vector3(x3,0.0,z1),
                 p110=Vector3(x3,0.0,0.0), p010=Vector3(x2,0.0,0.0),
                 p001=Vector3(x2,H,z2), p101=Vector3(x3,H,z1),
                 p111=Vector3(x3,H,0.0), p011=Vector3(x2,H,0.0))
vol3 = TFIVolume(p000=Vector3(x3,0.0,z1), p100=Vector3(x4,0.0,z1),
                 p110=Vector3(x4,0.0,0.0), p010=Vector3(x3,0.0,0.0),
                 p001=Vector3(x3,H,z1), p101=Vector3(x4,H,z1),
                 p111=Vector3(x4,H,0.0), p011=Vector3(x3,H,0.0))
from gdtk.geom.cluster import RobertsFunction
cf_j = RobertsFunction(True, False, 1.05)
cf_k = RobertsFunction(True, False, 1.05)
grd0 = StructuredGrid(pvolume=vol0, niv=61, njv=61, nkv=61, cf_list=[None, cf_j, cf_k])
b0 = FluidBlock(i=0, grid=grd0, initialState=inflow,
                bcs={'iminus':InflowFunctionBC('laminar_boundary_layer'),
                     'iplus':ExchangeBC(),
                     'kminus':WallNoSlipFixedTBC(300.0)})
grd1 = StructuredGrid(pvolume=vol1, niv=121, njv=61, nkv=61, cf_list=[None, cf_j, cf_k])
b1 = FluidBlock(i=1, grid=grd1, initialState=inflow,
                bcs={'iminus':ExchangeBC(),
                     'iplus':ExchangeBC(),
                     'jminus':WallNoSlipFixedTBC(300.0),
                     'kminus':WallNoSlipFixedTBC(300.0)})
grd2 = StructuredGrid(pvolume=vol2, niv=41, njv=61, nkv=61, cf_list=[None, cf_j, cf_k])
b2 = FluidBlock(i=2, grid=grd2, initialState=inflow,
                bcs={'iminus':ExchangeBC(),
                     'iplus':ExchangeBC(),
                     'jminus':WallNoSlipFixedTBC(300.0),
                     'kminus':WallNoSlipFixedTBC(300.0)})
grd3 = StructuredGrid(pvolume=vol3, niv=41, njv=61, nkv=61, cf_list=[None, cf_j, cf_k])
b3 = FluidBlock(i=3, grid=grd3, initialState=inflow,
                bcs={'iminus':ExchangeBC(),
                     'iplus':OutflowBC(),
                     'jminus':WallNoSlipFixedTBC(300.0),
                     'kminus':WallNoSlipFixedTBC(300.0)})

