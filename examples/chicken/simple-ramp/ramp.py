# ramp.py
# A simple 3D simulation of flow over a ramp with 10-degree deflection.
# PJ 2022-10-03 adapted from the Eilmer example.
# We'll model the domain as a pair of simple boxes.
Z = 0.0
L1 = 0.2
L2 = 1.0
W = 0.1
import math
H1 = 0.8*math.tan(math.radians(10.0))
H2 = 1.0

vol0 = TFIVolume(p000=Vector3(Z,Z,Z),   p100=Vector3(L1,Z,Z),
                 p110=Vector3(L1,W,Z),  p010=Vector3(Z,W,Z),
                 p001=Vector3(Z,Z,H2),  p101=Vector3(L1,Z,H2),
                 p111=Vector3(L1,W,H2),  p011=Vector3(Z,W,H2))
grd0 = StructuredGrid(pvolume=vol0, niv=11, njv=5, nkv=41)
vol1 = TFIVolume(p000=Vector3(L1,Z,Z),  p100=Vector3(L2,Z,H1),
                 p110=Vector3(L2,W,H1), p010=Vector3(L1,W,Z),
                 p001=Vector3(L1,Z,H2), p101=Vector3(L2,Z,H2),
                 p111=Vector3(L2,W,H2), p011=Vector3(L1,W,H2))
grd1 = StructuredGrid(pvolume=vol1, niv=31, njv=5, nkv=41)
#
initial = FlowState(p=5595.0, T=304.0)
inflow = FlowState(p=95.84e3, T=1103.0, velx=1000.0)
#
b0 = FluidBlock(i=0, grid=grd0, initialState=initial,
                bcs={'iminus':InflowBC(inflow),'iplus':ExchangeBC()})
b1 = FluidBlock(i=1, grid=grd1, initialState=initial,
                bcs={'iminus':ExchangeBC(),'iplus':OutflowBC()})
#
config.max_time = 5.0e-3
config.max_step = 1000
