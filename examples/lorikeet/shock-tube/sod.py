# sod.py
# Simulate the Sod shock tube with lorikeet.
# PJ 2022-12-11
# We'll model the shock tube as a pair of simple boxes.
Z = 0.0
L = 1.0
H = 0.1
box0 = CoonsPatch(p00=Vector3(Z,Z),   p10=Vector3(L/2,Z),
                  p11=Vector3(L/2,H), p01=Vector3(Z,H))
grd0 = StructuredGrid(psurf=box0, niv=51, njv=3)
box1 = CoonsPatch(p00=Vector3(L/2,Z), p10=Vector3(L,Z),
                  p11=Vector3(L,H),   p01=Vector3(L/2,H))
grd1 = StructuredGrid(psurf=box1, niv=51, njv=3)
#
init_gas_model("ideal-air-gas-model.lua")
high_pressure = FlowState(p=100.0e3, T=348.4)
low_pressure = FlowState(p=10.0e3, T=278.8)
#
b0 = FluidBlock(i=0, grid=grd0, initialState=high_pressure, bcs={'iplus':ExchangeBC()})
b1 = FluidBlock(i=1, grid=grd1, initialState=low_pressure, bcs={'iminus':ExchangeBC()})
#
config.max_time = 0.6e-3
config.max_step = 1000
