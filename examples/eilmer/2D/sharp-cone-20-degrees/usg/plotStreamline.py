#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Example script for plotting streamline data
# author: Kyle Damm, June 2016

import matplotlib.pyplot as plt
from numpy import *
# expected solution
p1 = 95.84e3                # Pa
p2 = 153e3                  # Pa
x_shock = 0.297794758       # m
xa = linspace(0, 1.0, 10000)
pa = []
for i in xa:
    if i < x_shock:
        pa.append(p1)
    elif i >= x_shock:
        pa.append(p2)

# initialise arrays
x = []
p = []

# read streamline.dat
with open('streamline.dat', 'r') as f:
    data = f.readlines()
    count = 0
    for line in data:
        dat = line.split()
        if (count>0 and dat[0] == '#'):
            print line
        if (count>1):
            x.append(float(dat[0]))
            p.append(float(dat[8]))
        count += 1

# plot streamline data
plt.plot(x, p, 'k+')
plt.plot(xa, pa, 'k')
plt.title('Pressure along a streamline')
plt.xlabel('x (m)')
plt.ylabel('p (Pa)')
plt.ylim(80000,)
plt.show()
