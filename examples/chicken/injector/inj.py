# inj.py
#
# This script creates and runs a simulation of transverse injection
# into a Mach 4 stream.  The boundary-layer inflow definition comes
# from the baoundary-layer short-plate example.
# The main flow duct is aligned with the i-index direction and
# the boundary layer grows on the 'kminus' boundary.
#
# PJ 2023-05-15 Adapted from short-plate boundary-layer example.
#
config.title = "Transverse injection into a supersonic flow."
config.viscous = True
config.max_time = 0.5e-3
config.max_step = 300000
config.dt_init = 1.0e-8
add_cfl_value(0.0, 0.5)
add_dt_plot(0.0, 0.1e-3)
#
inflow = FlowState(p=1.013e3, T=300.0, velx=1390.0)
injector = FlowState(p=10.13e3, T=300.0, vely=350.0)
#
d = 0.005  # Scale everything on the injector hole size
X0 = 0.0; X1 = 10*d; X2 = X1 + d; X3 = X2 + 16*d
Z0 = 0.0; Z1 = d/2; Z2 = Z1 + 5*d
Y0 = 0.0; Y1 = 7*d
from gdtk.geom.cluster import RobertsFunction
# Inflow to injector
vol00 = TFIVolume(p000=Vector3(X0,Y0,Z2), p100=Vector3(X1,Y0,Z2),
                  p110=Vector3(X1,Y0,Z1), p010=Vector3(X0,Y0,Z1),
                  p001=Vector3(X0,Y1,Z2), p101=Vector3(X1,Y1,Z2),
                  p111=Vector3(X1,Y1,Z1), p011=Vector3(X0,Y1,Z1))
cf00_i = RobertsFunction(False, True, 1.20)
cf00_j = RobertsFunction(False, True, 1.20)
cf00_k = RobertsFunction(True, False, 1.05)
grd00 = StructuredGrid(pvolume=vol00, niv=101, njv=51, nkv=71,
                       cf_list=[cf00_i,cf00_j,cf00_k])
blk00 = FluidBlock(i=0, j=0, grid=grd00, initialState=inflow,
                  bcs={'iminus':InflowFunctionBC('laminar_boundary_layer'),
                       'iplus':ExchangeBC(),
                       'jplus':ExchangeBC(),
                       'kminus':WallNoSlipFixedTBC(300.0)})
vol01 = TFIVolume(p000=Vector3(X0,Y0,Z1), p100=Vector3(X1,Y0,Z1),
                  p110=Vector3(X1,Y0,Z0), p010=Vector3(X0,Y0,Z0),
                  p001=Vector3(X0,Y1,Z1), p101=Vector3(X1,Y1,Z1),
                  p111=Vector3(X1,Y1,Z0), p011=Vector3(X0,Y1,Z0))
grd01 = StructuredGrid(pvolume=vol01, niv=101, njv=11, nkv=71,
                       cf_list=[cf00_i,None,cf00_k])
blk01 = FluidBlock(i=0, j=1, grid=grd01, initialState=inflow,
                  bcs={'iminus':InflowFunctionBC('laminar_boundary_layer'),
                       'iplus':ExchangeBC(),
                       'jminus':ExchangeBC(),
                       'kminus':WallNoSlipFixedTBC(300.0)})
# Injector plane
vol10 = TFIVolume(p000=Vector3(X1,Y0,Z2), p100=Vector3(X2,Y0,Z2),
                  p110=Vector3(X2,Y0,Z1), p010=Vector3(X1,Y0,Z1),
                  p001=Vector3(X1,Y1,Z2), p101=Vector3(X2,Y1,Z2),
                  p111=Vector3(X2,Y1,Z1), p011=Vector3(X1,Y1,Z1))
grd10 = StructuredGrid(pvolume=vol10, niv=21, njv=51, nkv=71,
                       cf_list=[None,cf00_j,cf00_k])
blk10 = FluidBlock(i=1, j=0, grid=grd10, initialState=inflow,
                  bcs={'iminus':ExchangeBC(),
                       'iplus':ExchangeBC(),
                       'jplus':ExchangeBC(),
                       'kminus':WallNoSlipFixedTBC(300.0)})
vol11 = TFIVolume(p000=Vector3(X1,Y0,Z1), p100=Vector3(X2,Y0,Z1),
                  p110=Vector3(X2,Y0,Z0), p010=Vector3(X1,Y0,Z0),
                  p001=Vector3(X1,Y1,Z1), p101=Vector3(X2,Y1,Z1),
                  p111=Vector3(X2,Y1,Z0), p011=Vector3(X1,Y1,Z0))
grd11 = StructuredGrid(pvolume=vol11, niv=21, njv=11, nkv=71,
                       cf_list=[None,None,cf00_k])
blk11 = FluidBlock(i=1, j=1, grid=grd11, initialState=inflow,
                  bcs={'iminus':ExchangeBC(),
                       'iplus':ExchangeBC(),
                       'jminus':ExchangeBC(),
                       'kminus':InflowBC(injector)})
# Downstream of injector
vol20 = TFIVolume(p000=Vector3(X2,Y0,Z2), p100=Vector3(X3,Y0,Z2),
                  p110=Vector3(X3,Y0,Z1), p010=Vector3(X2,Y0,Z1),
                  p001=Vector3(X2,Y1,Z2), p101=Vector3(X3,Y1,Z2),
                  p111=Vector3(X3,Y1,Z1), p011=Vector3(X2,Y1,Z1))
cf20_i = RobertsFunction(True, False, 1.20)
grd20 = StructuredGrid(pvolume=vol20, niv=161, njv=51, nkv=71,
                       cf_list=[cf20_i,cf00_j,cf00_k])
blk20 = FluidBlock(i=2, j=0, grid=grd20, initialState=inflow,
                  bcs={'iminus':ExchangeBC(),
                       'iplus':OutflowBC(),
                       'jplus':ExchangeBC(),
                       'kminus':WallNoSlipFixedTBC(300.0)})
vol21 = TFIVolume(p000=Vector3(X2,Y0,Z1), p100=Vector3(X3,Y0,Z1),
                  p110=Vector3(X3,Y0,Z0), p010=Vector3(X2,Y0,Z0),
                  p001=Vector3(X2,Y1,Z1), p101=Vector3(X3,Y1,Z1),
                  p111=Vector3(X3,Y1,Z0), p011=Vector3(X2,Y1,Z0))
grd21 = StructuredGrid(pvolume=vol21, niv=161, njv=11, nkv=71,
                       cf_list=[cf20_i,None,cf00_k])
blk21 = FluidBlock(i=2, j=1, grid=grd21, initialState=inflow,
                  bcs={'iminus':ExchangeBC(),
                       'iplus':OutflowBC(),
                       'jminus':ExchangeBC(),
                       'kminus':WallNoSlipFixedTBC(300.0)})

