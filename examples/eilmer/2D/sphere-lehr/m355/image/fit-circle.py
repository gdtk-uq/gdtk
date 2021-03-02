#! /usr/bin/env python3
# fit-circle.py
# Least-squares fitting of a circle to Lehr body points, 
# as described in PJ's workbook 29-Oct-2018 (page 9)
# PJ 2020-03-10

import numpy as np
import math

lines = open('body-points.csv', 'r').read().split('\n')
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

N = len(x)
Sx = sum(x); Sy = sum(y)
Sx2 = sum(x*x); Sy2 = sum(y*y)
Sxy = sum(x*y);

A = np.array([[N,  Sx,  Sy ],
              [Sx, Sx2, Sxy],
              [Sy, Sxy, Sy2]])
print('A=', A)
b = np.array([-sum(x*x + y*y),
              -sum(x*(x*x + y*y)),
              -sum(y*(x*x + y*y))])
print('b=', b)

alpha = np.linalg.solve(A, b)
print('alpha=', alpha)
x0 = -alpha[1]/2; y0 = -alpha[2]/2
R0 = math.sqrt((x0*x0 + y0*y0) - alpha[0])
print('x0=', x0, 'y0=', y0, 'R0=', R0)

err2 = (x-x0)*(x-x0) + (y-y0)*(y-y0) - R0*R0
print('err2=', err2)
