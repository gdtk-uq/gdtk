# surface-test.py
#
# PJ, 2020-07-05

import math
from gdtk.geom.vector3 import Vector3, approxEqualVectors
from gdtk.geom.path import *
from gdtk.geom.surface import *

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

print("Done.")
