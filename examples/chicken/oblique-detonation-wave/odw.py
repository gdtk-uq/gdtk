# odw.py
# Oblique detonation wave with a simple reacting gas.
# Modelled on Rowan's example for Eilmer3, way back in the year 2006.
# Peter J. 2022-11-15
#
config.title = "Oblique detonation wave."
config.dt_init = 1.0e-6
config.max_time = 2.0e-2
config.max_step = 300000
add_dt_plot(0.0, config.max_time/40)
#
inflow = FlowState(p=86.1e3, T=300.0, velx=964.302)
initial = FlowState(p=28.7e3, T=300.0, velx=0.0)
#
# Geometry
from gdtk.geom.path import Line, FnPath
from gdtk.geom.surface import CoonsPatch
from gdtk.geom.volume import SweptSurfaceVolume
#
xmin = -0.25; xmax = 1.75
ymin = 0.0; ymax = 2.0
zmin = 0.0; zmax = 0.1
#
from oblique_detonation import *
from math import pi
od =  ObliqueDetonation(pi/4.0, 300.0, 3.0, 1.0)
wall = FnPath(od.create_wall_function(0.0, xmax))
#
# We are going to build two patches in the x,y plane and then extrude
# these patches into volumes.
a0 = Vector3(xmin, 0.0); a1 = Vector3(xmin, ymax)
b0 = Vector3(0.0, 0.0);  b1 = Vector3(0.0, ymax)
c0 = wall(1.0); c1 = Vector3(xmax, ymax)
#
patch0 = CoonsPatch(p00=a0, p10=b0, p11=b1, p01=a1)
patch1 = CoonsPatch(west=Line(b0,b1), east=Line(c0,c1),
                    north=Line(b1,c1), south=wall)
#
dz = Vector3(0.0, 0.0, zmax)
vol0 = SweptSurfaceVolume(face0123=patch0, edge04=Line(a0, a0+dz))
vol1 = SweptSurfaceVolume(face0123=patch1, edge04=Line(b0, b0+dz))
#
nn = 64
nnx = nn; nny = nn # overall numbers of cells
nnx0 = int(0.125*nnx) # share for blk0
nnx1 = nnx - int(0.125*nnx) # share for blk1
grd0 = StructuredGrid(pvolume=vol0, niv=nnx0+1, njv=nny+1, nkv=3)
grd1 = StructuredGrid(pvolume=vol1, niv=nnx1+1, njv=nny+1, nkv=3)
#
blk0 = FluidBlock(i=0, grid=grd0, initialState=inflow,
                  bcs={'iminus':InflowBC(inflow),'iplus':ExchangeBC()})
blk1 = FluidBlock(i=1, grid=grd1, initialState=inflow,
                  bcs={'iminus':ExchangeBC(),'iplus':OutflowBC()})
