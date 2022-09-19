# sod.py
# Simulate the Sod shock tube with chicken.
# PJ 2022-09-19
# We'll model the shock tube as a pair of simple boxes.
Z = 0.0
L = 1.0
H = 0.1
vol0 = TFIVolume(p000=Vector3(Z,Z,Z),   p100=Vector3(L/2,Z,Z),
                 p110=Vector3(L/2,H,Z), p010=Vector3(Z,H,Z),
                 p001=Vector3(Z,Z,H),   p101=Vector3(L/2,Z,H),
                 p111=Vector3(L/2,H,H), p011=Vector3(Z,H,H))
grd0 = StructuredGrid(pvolume=vol0, niv=11, njv=5, nkv=5)
vol1 = TFIVolume(p000=Vector3(L/2,Z,Z), p100=Vector3(L,Z,Z),
                 p110=Vector3(L,H,Z),   p010=Vector3(L/2,H,Z),
                 p001=Vector3(L/2,Z,H), p101=Vector3(L,Z,H),
                 p111=Vector3(L,H,H),   p011=Vector3(L/2,H,H))
grd1 = StructuredGrid(pvolume=vol1, niv=11, njv=5, nkv=5)
#
high_pressure = FlowState(p=100.0e3, T=300.0)
low_pressure = FlowState(p=20.0e3, T=270.0)
#
b0 = FluidBlock(i=0, grid=grd0, initialState=high_pressure)
b1 = FluidBlock(i=1, grid=grd1, initialState=low_pressure)
