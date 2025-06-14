#! /usr/bin/env python3
#
# estimate_shock_angle.py
# An example of custom postprocessing for an Eilmer5 simulation.
# Peter J. 2025-06-15
#
from gdtk.lmr import LmrConfig, SimInfo
import matplotlib.pyplot as plt
import scipy.stats.mstats as mstats
import math

lmrcfg = LmrConfig()
sim = SimInfo(lmrcfg)
snap_last = sim.read_snapshot(sim.snapshots[-1])

# Decide a trigger-level for detecting the shock.
xs, ys, ps = snap_last.get_slice(var='p', j=10)
p_min = ps[0] # Free-stream pressure
p_max = max(ps) # Post-shock pressure
p_level = p_min + 0.3*(p_max - p_min)

def locate_shock_along_strip(xs, ys, ps, plevel):
    x = xs[0]; y = ys[0]
    for i in range(1, len(ps)):
        if ps[i] >= p_level:
            frac = (p_level - ps[i-1])/(ps[i] - ps[i-1])
            x = xs[i-1]*(1.0-frac) + xs[i]*frac
            y = ys[i-1]*(1.0-frac) + ys[i]*frac
            break
    return x, y

x_loc = []; y_loc = []
for j in range(sim.grids[0].njc):
    xs, ys, ps = snap_last.get_slice(var='p', j=j)
    if ps[-1] > p_level:
        x, y = locate_shock_along_strip(xs, ys, ps, p_level)
        x_loc.append(x); y_loc.append(y)

# print("x_loc, y_loc=", x_loc, y_loc)
result = mstats.linregress(x_loc, y_loc)
print(f"shock_angle= {math.degrees(math.atan(result.slope))} degrees")

fig, ax = plt.subplots()
ax.set_title("Shock location")
ax.set_xlabel('x, m'); ax.set_ylabel('y, m')
ax.plot(x_loc, y_loc, '+')
plt.show()
