# sod_j.py
# Simulate the Sod shock tube with chicken -- tube extends in j-index direction.
# PJ 2022-10-19
# We'll model the shock tube as a pair of simple boxes.
Z = 0.0
L = 1.0
H = 0.1
vol0 = TFIVolume(p000=Vector3(Z,Z,Z),   p100=Vector3(H,Z,Z),
                 p110=Vector3(H,L/2,Z), p010=Vector3(Z,L/2,Z),
                 p001=Vector3(Z,Z,H),   p101=Vector3(H,Z,H),
                 p111=Vector3(H,L/2,H), p011=Vector3(Z,L/2,H))
grd0 = StructuredGrid(pvolume=vol0, niv=3, njv=51, nkv=3)
vol1 = TFIVolume(p000=Vector3(Z,L/2,Z), p100=Vector3(H,L/2,Z),
                 p110=Vector3(H,L,Z),   p010=Vector3(Z,L,Z),
                 p001=Vector3(Z,L/2,H), p101=Vector3(H,L/2,H),
                 p111=Vector3(H,L,H),   p011=Vector3(Z,L,H))
grd1 = StructuredGrid(pvolume=vol1, niv=3, njv=51, nkv=3)
#
high_pressure = FlowState(p=100.0e3, T=348.4)
low_pressure = FlowState(p=10.0e3, T=278.8)
#
b0 = FluidBlock(j=0, grid=grd0, initialState=high_pressure, bcs={'jplus':ExchangeBC()})
b1 = FluidBlock(j=1, grid=grd1, initialState=low_pressure, bcs={'jminus':ExchangeBC()})
#
config.max_time = 0.6e-3
config.max_step = 1000
