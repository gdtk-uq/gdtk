# plot_surface_conditions.py
# Assemble the ramp loads by picking up the loads files.
# Peter J. 2024-03-13, 2025-08-27 simplify the reading of the dataframe
# RJG, 2024-03-19
#   + change to pyplot usage
#   + add option to save figure
#
# Usage:
#
# For plotting to screen:
# $ python[3] plot_surface_conditions.py
#
# For saving figure to disk:
# $ python[3] plot_surface_conditions.py --save-figure
#
from gdtk.lmr import LmrConfig, SimInfo
import numpy as np
import matplotlib.pyplot as plt
import sys

lmrcfg = LmrConfig()
sim = SimInfo(lmrcfg)
# Pick up the final loads files as a Pandas DataFrame.
df = sim.read_loads(indx=sim.loads_indices[-1], group="loads")
df = df.sort_values(by=['pos.x'])

xp = np.loadtxt('mohammadian-figure-12-p_p_inf.data')
xq = np.loadtxt('mohammadian-figure-13-heat-flux.data')

ig, (ax0, ax1) = plt.subplots(2,1)
ax0.plot(df['pos.x']*1000, df['q_total']/1000, '-b', label='Eilmer')
ax0.plot(xq[:,0]*25.4, xq[:,1]*11.5, '+k', label='Mohammadian (1972)')
ax0.set_xlabel(''); ax0.set_ylabel('q, kW/m^2')
ax0.set_xlim(0,250); ax0.set_ylim(0,120)
ax0.legend(loc='upper right')
ax0.set_title('Surface conditions on convex ramp')
ax1.plot(df['pos.x']*1000, df['p']/1000, '-b', label='Eilmer')
ax1.plot(xp[:,0]*25.4, xp[:,1]*66.43/1000, '+k', label='Mohammadian (1972)')
ax1.set_xlabel('x, mm'); ax1.set_ylabel('p, kPa')
ax1.set_xlim(0,250); ax1.set_ylim(0,2.5)
ax1.legend(loc='upper right')
ax1.set_title('')
if len(sys.argv) == 2 and sys.argv[1] == "--save-figure":
    fig.savefig("convex-ramp-surf-conds.png", transparent=True, dpi=300)
else:
    plt.show()
