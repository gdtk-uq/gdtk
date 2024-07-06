#!/usr/bin/python3
"""
Calls the busemann.py module to generate a contour file.

Usage:

  > ./generate-contour.py

.. Author: RJG
.. Date: 2024-07-06
"""

from gdtk.busemann import BusemannDiffuser
from math import asin

M2 = 3.0
k = 1.4
theta_23 = asin(k/M2)
bd = BusemannDiffuser(M2, theta_23)

r = 1.0
dtheta = 0.001
bd.generate_contour(r, dtheta)
npts = 200
bd.write_contour('bd-contour.dat', npts)

props = bd.properties()

with open('bd-props.txt', "w") as f:
    for k, v in zip(props._fields, props):
        f.write(f"{k} = {v:.6f}\n")

  
  

