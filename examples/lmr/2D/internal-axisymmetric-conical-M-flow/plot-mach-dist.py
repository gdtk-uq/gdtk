import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'monospace'

tm_soln = np.loadtxt('mach-distribution.dat', skiprows=1, unpack=True)
w0 = np.loadtxt('wall-0.dat', skiprows=1, usecols=(0, 13), unpack=True)
w1 = np.loadtxt('wall-1.dat', skiprows=1, usecols=(0, 13), unpack=True)
w2 = np.loadtxt('wall-2.dat', skiprows=1, usecols=(0, 13), unpack=True)
w3 = np.loadtxt('wall-3.dat', skiprows=1, usecols=(0, 13), unpack=True)

fig, ax = plt.subplots()

ax.plot(tm_soln[0], tm_soln[1], "-", linewidth=2, label='Taylor-Maccoll solution for M-flow')
ax.plot(w0[0], w0[1], "s", markerfacecolor='none', markersize=4, label="wall-adjacent cells")
ax.plot(w1[0], w1[1], "o", markerfacecolor='none', markersize=4, label="+1 off wall")
ax.plot(w2[0], w2[1], "^", markerfacecolor='none', markersize=4, label="+2 off wall")
ax.plot(w3[0], w3[1], "D", markerfacecolor='none', markersize=4, label="+3 off wall")
ax.set_xlabel('x')
ax.set_ylabel('Mach number')
ax.legend()

plt.savefig("mach-distribution.pdf")
