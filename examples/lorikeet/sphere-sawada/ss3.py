# ss3.lua
#
# Sphere in equilibrium air modelling Case 3 from
#    K. Sawada & E. Dendou (2001)
#    Validation of hypersonic chemical equilibrium flow calculations
#    using ballistic-range data.
#    Shock Waves (2001) Vol. 11, pp 43--51
#
# Experimental shock stand-off distance is 2.59mm
# Sawada & Dendou CFD value:               2.56mm
#
# Authors: PAJ and RJG
# Versions:
#    Tcl version: 22-Jan-2004, derived from rbody.
#    Python version: ss3.py, 04-Apr-2005, 10-Aug-2006, 27-Nov-2006
#                    12-Nov-2008 by RJG for use in Elmer3
#    Lua version: 22-Aug-2016 by PJ for use in Eilmer4
#                 2019-05-26 shape grid with Billig correlation.
#    Back to a Python version for Lorikeet, 2022-01-02

config.title = "Sphere in hypersonic flow of air in chemical equilibrium."
print(config.title)
config.axisymmetric = True

init_gas_model('cea-lut-air.lua')
# Free-stream flow definition
p_inf = 20.0e3   # Pa
T_inf = 296.0    # degrees K
vx_inf = 4.68e3  # flow speed, m/s
inflow = FlowState(p=p_inf, T=T_inf, velx=vx_inf)
initial = FlowState(p=0.3*p_inf, T=T_inf)

print("Building grid.")
from gdtk.billig_patch import make_patch
from gdtk.geom.cluster import RobertsFunction
R = 31.8e-3  # radius of sphere, in metres
M_inf = vx_inf/inflow.gas.a
print("M_inf=", M_inf)
bp = make_patch(Minf=M_inf, R=R, scale=0.8)
cf_circum = RobertsFunction(end0=True, end1=False, beta=1.1)
grid = StructuredGrid(psurf=bp['patch'], niv=61, njv=61,
		      cf_list=[None, cf_circum, None, cf_circum])
blk0 = FluidBlock(grid=grid, initialState=initial,
		  bcs={'west':InflowBC(inflow), 'north':OutflowBC()})
# We have left east and south as (default) slip-walls

# Set a few more config options
config.max_time = 15.0*R/vx_inf # allow time to settle at nose
config.max_step = 10000000
config.dt_init = 1.0e-9
config.plot_dt = config.max_time/8
