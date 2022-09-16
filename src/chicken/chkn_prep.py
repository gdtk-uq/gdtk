# chkn_prep.py
# Prepare the files for a chicken run.
# Run with:
# $ python3 chkn_prep.py
#
# PJ 2022-09-16
# We'll start with just the build of a 3D grid for a simple box.

from gdtk.geom.vector3 import Vector3
from gdtk.geom.volume import TFIVolume
from gdtk.geom.sgrid import StructuredGrid

Z = 0.0
L = 1.0
H = 0.1
vol = TFIVolume(p000=Vector3(Z,Z,Z), p100=Vector3(L,Z,Z),
                p110=Vector3(L,H,Z), p010=Vector3(Z,H,Z),
                p001=Vector3(Z,Z,H), p101=Vector3(L,Z,H),
                p111=Vector3(L,H,H), p011=Vector3(Z,H,H))
print("vol=", vol)
c = vol(0.5,0.5,0.5)
print("c=", c)
grd = StructuredGrid(pvolume=vol, niv=21, njv=5, nkv=5)
grd.write_to_vtk_file("chkn_grid.vtk")
