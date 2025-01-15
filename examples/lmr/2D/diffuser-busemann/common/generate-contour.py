#!/usr/bin/env python
"""
Calls the busemann.py module to generate a contour file.

Usage:

  > ./generate-contour.py npts

  where 'npts' is an integer number of points used to
  represent the contour.

.. Author: RJG
.. Date: 2024-07-06
"""

from gdtk.busemann import BusemannDiffuser
from math import asin
import sys

if len(sys.argv) != 2:
    print("Incorrect number of arguments; exactly 1 expected.")
    print("Usage: ./generate-contour.py npts")
    print("")
    print("   where 'npts' is an integer")
    sys.exit(1)

npts = int(sys.argv[1])

M2 = 3.0
k = 1.4
theta_23 = asin(k/M2)
bd = BusemannDiffuser(M2, theta_23)

r = 1.0
dtheta = 0.001
bd.generate_contour(r, dtheta)
bd.write_contour('bd-contour.dat', npts)

props = bd.properties()

with open('bd-props.txt', "w") as f:
    for k, v in zip(props._fields, props):
        f.write(f"{k} = {v:.6f}\n")

  
  

