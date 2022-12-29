# odw.py
# Oblique detonation wave with a simple reacting gas.
# Modelled on Rowan's example for Eilmer3, way back in the year 2006.
# Peter J. 2022-12-29 Adapted from the Chicken example.
#
config.title = "Oblique detonation wave with Powers-Aslam gas model."
print(config.title)
#
init_gas_model("powers-aslam-gas-model.lua")
inflow = FlowState(p=86.1e3, T=300.0, velx=964.302, massf=[1.0, 0.0])
initial = FlowState(p=28.7e3, T=300.0, velx=0.0, massf=[1.0, 0.0])
config.reacting = True
#
# Geometry
from gdtk.geom.path import Line, FnPath
from gdtk.geom.surface import CoonsPatch
#
xmin = -0.25; xmax = 1.75
ymin = 0.0; ymax = 2.0
#
from oblique_detonation import *
from math import pi
od =  ObliqueDetonation(pi/4.0, 300.0, 3.0, 1.0)
wall = FnPath(od.create_wall_function(0.0, xmax))
#
# We are going to build two patches in the x,y plane.
a0 = Vector3(xmin, 0.0); a1 = Vector3(xmin, ymax)
b0 = Vector3(0.0, 0.0);  b1 = Vector3(0.0, ymax)
c0 = wall(1.0); c1 = Vector3(xmax, ymax)
#
patch0 = CoonsPatch(p00=a0, p10=b0, p11=b1, p01=a1)
patch1 = CoonsPatch(west=Line(b0,b1), east=Line(c0,c1),
                    north=Line(b1,c1), south=wall)
#
# Mesh the patches, with particular discretisation.
#
factor = 2
nn = factor*40
nnx = nn; nny = nn # overall numbers of cells
fraction0 = (0-xmin)/(xmax-xmin) # fraction of domain upstream of wedge
nnx0 = int(fraction0*nnx) # share for blk0
nnx1 = nnx - nnx0 # share for blk1
grd0 = StructuredGrid(psurf=patch0, niv=nnx0+1, njv=nny+1)
grd1 = StructuredGrid(psurf=patch1, niv=nnx1+1, njv=nny+1)
#
blk0 = FluidBlock(i=0, grid=grd0, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow),'iplus':ExchangeBC()})
blk1 = FluidBlock(i=1, grid=grd1, initialState=initial,
                  bcs={'iminus':ExchangeBC(),'iplus':OutflowBC()})
#
config.flux_calc = FluxCalc.ausmdv_plus_hanel
config.max_time = 2.0e-2
config.max_step = 300000
config.plot_dt = config.max_time/40
