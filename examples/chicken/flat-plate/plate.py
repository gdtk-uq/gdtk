# plate.py
#
# This script creates and runs a simulation of a worked
# example (Example 5-11 on page 155) in Joseph Schetz's
# book "Boundary Layer Analysis".
#
# The example involves a laminar Mach 4.0 flow over a
# 1.0 m long flat plate. The boundary layer will grow on
# the 'kminus' boundary. Simulation results can be compared
# with those from the CLBL code from Schetz's book.
#
# Peter Jacobs & Wilson Chan,
# 2010-01-18 Python version for eilmer3
# 2015-06-07 Lua version for Eilmer4
# 2022-10-03 Python version for Chicken.
#
config.title = "Supersonic flow over a flat plate."
config.viscous = True
config.max_time = 2.5e-3
config.max_step = 300000
config.dt_init = 1.0e-9
add_cfl_value(0.0, 0.5)
add_dt_plot(0.0, 0.5e-3)
#
inflow = FlowState(p=1.013e3/10, T=300.0, velx=1390.0) # Lower pressure to increase BL thickness.
#
L = 1.1; H = 0.4*L; W = 0.1
vol0 = TFIVolume(p000=Vector3(0,0,W),      p100=Vector3(L,0,W),
                 p110=Vector3(L,0,0),      p010=Vector3(0,0,0),
                 p001=Vector3(0,0.25*H,W), p101=Vector3(L,H,W),
                 p111=Vector3(L,H,0),      p011=Vector3(0,0.25*H,0))
from gdtk.geom.cluster import RobertsFunction
cf_i = RobertsFunction(True, False, 1.05)
cf_k = RobertsFunction(True, False, 1.05)
grd0 = StructuredGrid(pvolume=vol0, niv=61, njv=3, nkv=81, cf_list=[cf_i, None, cf_k])
b0 = FluidBlock(i=0, grid=grd0, initialState=inflow,
                bcs={'iminus':InflowBC(inflow),'iplus':OutflowBC(),
                     'kminus':WallNoSlipFixedTBC(),'kplus':InflowBC(inflow)})
