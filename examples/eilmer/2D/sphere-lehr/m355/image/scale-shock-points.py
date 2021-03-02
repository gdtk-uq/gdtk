#! /usr/bin/env python3
# scale-shock-points.py
# Adjust the coordinates of the digitized shock points
# so that origin is at the centre of the sphere.
# PJ 2020-03-10

import numpy as np
import math

lines = open('shock-points.csv', 'r').read().split('\n')
del lines[0]
x = []; y=[]
for line in lines:
    items = line.split()
    if len(items) < 2: break
    x.append(float(items[0]))
    y.append(float(items[1]))
x = np.array(x); y=np.array(y)
print('x=', x)
print('y=', y)

# In digitized units, we have from fitting the circle.
x0= 1.15390879046
y0= 0.0562897207593
R0= 0.7792234770096957

# In simulation, we set the nose radius.
Rn = 7.5e-3 # metres

# Shift to the body-centre coordinates.
xd = x - x0; yd = y - y0
r = np.sqrt(xd*xd + yd*yd)*Rn/R0
theta = np.arctan2(yd, xd)
# theta -= 1.0/27 # to rotate the shock a little

# Coordinates in simulation space.
xs = r*np.cos(theta); ys = r*np.sin(theta)
# Write in VTK classic format.
N = len(xs)
with open('shock-points.vtk', 'w') as f:
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Points on shock\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS %d float\n' % N)
    for i in range(N):
        f.write('%g %g 0.0\n' % (xs[i], ys[i]))
    f.write('VERTICES %d %d\n' % (N, 2*N))
    for i in range(N):
        f.write('%d %d\n' % (1, i))

