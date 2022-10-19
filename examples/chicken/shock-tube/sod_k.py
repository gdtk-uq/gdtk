# sod_k.py
# Simulate the Sod shock tube with chicken -- tube extends in k-index direction.
# PJ 2022-10-19
# We'll model the shock tube as a pair of simple boxes.
Z = 0.0
L = 1.0
H = 0.1
vol0 = TFIVolume(p000=Vector3(Z,Z,Z),   p100=Vector3(H,Z,Z),
                 p110=Vector3(H,H,Z),   p010=Vector3(Z,H,Z),
                 p001=Vector3(Z,Z,L/2), p101=Vector3(H,Z,L/2),
                 p111=Vector3(H,H,L/2), p011=Vector3(Z,H,L/2))
grd0 = StructuredGrid(pvolume=vol0, niv=3, njv=3, nkv=51)
vol1 = TFIVolume(p000=Vector3(Z,Z,L/2), p100=Vector3(H,Z,L/2),
                 p110=Vector3(H,H,L/2), p010=Vector3(Z,H,L/2),
                 p001=Vector3(Z,Z,L),   p101=Vector3(H,Z,L),
                 p111=Vector3(H,H,L),   p011=Vector3(Z,H,L))
grd1 = StructuredGrid(pvolume=vol1, niv=3, njv=3, nkv=51)
#
high_pressure = FlowState(p=100.0e3, T=348.4)
low_pressure = FlowState(p=10.0e3, T=278.8)
#
b0 = FluidBlock(k=0, grid=grd0, initialState=high_pressure, bcs={'kplus':ExchangeBC()})
b1 = FluidBlock(k=1, grid=grd1, initialState=low_pressure, bcs={'kminus':ExchangeBC()})
#
config.max_time = 0.6e-3
config.max_step = 1000
