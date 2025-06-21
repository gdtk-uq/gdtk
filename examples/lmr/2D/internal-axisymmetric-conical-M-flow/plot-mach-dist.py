# Author: RJG
# Date: 2025-06-21

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.size'] = 14

tm_soln = np.loadtxt('mach-distribution.dat', skiprows=1, unpack=True)
w0 = np.loadtxt('wall-0.dat', skiprows=1, usecols=(0, 13), unpack=True)
w1 = np.loadtxt('wall-1.dat', skiprows=1, usecols=(0, 13), unpack=True)
w2 = np.loadtxt('wall-2.dat', skiprows=1, usecols=(0, 13), unpack=True)
w3 = np.loadtxt('wall-3.dat', skiprows=1, usecols=(0, 13), unpack=True)

cm = 1/2.54
fig, ax = plt.subplots(figsize=(20*cm, 14*cm))

# filter our pre-shock values by creating masks
M_inf = 2.6
m0 = w0[1] < M_inf
m1 = w1[1] < M_inf
m2 = w2[1] < M_inf
m3 = w3[1] < M_inf

ax.plot(tm_soln[0], tm_soln[1], "-", linewidth=2, label='Taylor-Maccoll solution for M-flow')
ax.plot(w0[0][m0], w0[1][m0], "s", markerfacecolor='none', markersize=4, label="wall-adjacent cells")
ax.plot(w1[0][m1], w1[1][m1], "o", markerfacecolor='none', markersize=4, label="+1 off wall")
ax.plot(w2[0][m2], w2[1][m2], "^", markerfacecolor='none', markersize=4, label="+2 off wall")
ax.plot(w3[0][m3], w3[1][m3], "D", markerfacecolor='none', markersize=4, label="+3 off wall")
ax.set_xlabel('x')
ax.set_ylabel('Mach number')
ax.set_ylim(2.1, 2.6)
ax.legend(fontsize='small', loc='lower right')

plt.savefig("mach-distribution.pdf")
