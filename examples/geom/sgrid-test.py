# sgrid-test.py
#
# PJ, 2020-07-05

import math
from eilmer.geom.vector3 import Vector3, approxEqualVectors
from eilmer.geom.path import *
from eilmer.geom.surface import *
from eilmer.geom.sgrid import StructuredGrid
from eilmer.geom.cluster import RobertsFunction

print("Begin test of ParametricSurface classes...")
p00 = Vector3([0.0, 0.1, 3.0])
p10 = Vector3(1.0, 0.1, 3.0)
p11 = Vector3(1.0, 1.1, 3.0)
p01 = Vector3(0.0, 1.1, 3.0)
my_patch = CoonsPatch(p00=p00, p10=p10, p11=p11, p01=p01)
c = my_patch(0.5, 0.5)
assert approxEqualVectors(c, Vector3(0.5, 0.6, 3.0)), "CoonsPatch middle"
c = my_patch(0.1, 0.1)
assert approxEqualVectors(c, Vector3(0.1, 0.2, 3.0)), "CoonsPatch near corner"
print(f"my_patch={my_patch}, c={c}")

cf_y = RobertsFunction(True, False, 1.1)
g = StructuredGrid(psurf=my_patch, niv=11, njv=5,
                   cf_list=[None, cf_y, None, cf_y])
print("g=", g)
g.write_to_vtk_file("my_grid.vtk")

print("Done.")
