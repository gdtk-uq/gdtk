# plot_surface_conditions.py
# Assemble the ramp loads by picking data out of the loads files.
# Peter J. 2024-03-13
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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
#
# We know that the ramp surface consists of the southern boundaries
# of the lower blocks.  These will be the even-number blocks.
blk_ids = list(range(0,27,2))
#
# For the final snapshot, assemble the data for southern boundaries
# into a single data frame.
meta_df = pd.read_table('./lmrsim/loads/loads-metadata', sep='\s+')
snaps = meta_df['loads_index']
final_snap = snaps[snaps.index[-1]]
blk_dfs = []
for i in blk_ids:
    f = './lmrsim/loads/%04d/blk-%04d-bndry-2-loads.dat' % (final_snap, i)
    blk_dfs.append(pd.read_table(f, sep='\s+'))
ramp_df = pd.concat(blk_dfs)

xp = np.loadtxt('mohammadian-figure-12-p_p_inf.data')
xq = np.loadtxt('mohammadian-figure-13-heat-flux.data')

plt.rcParams.update({'font.size': 14,
                     'font.family': 'monospace',
                     'font.monospace': ['Courier New']})
fig, (ax0, ax1) = plt.subplots(2,1)
ax0.plot(ramp_df['pos.x']*1000, ramp_df['q_total']/1000, '-b', label='Eilmer')
ax0.plot(xq[:,0]*25.4, xq[:,1]*11.5, '+k', label='Mohammadian (1972)')
ax0.set_xlabel(''); ax0.set_ylabel('q, kW/m^2')
ax0.set_xlim(0,250); ax0.set_ylim(0,120)
ax0.legend(loc='upper right')
ax0.set_title('Surface conditions on convex ramp')
ax1.plot(ramp_df['pos.x']*1000, ramp_df['p']/1000, '-b', label='Eilmer')
ax1.plot(xp[:,0]*25.4, xp[:,1]*66.43/1000, '+k', label='Mohammadian (1972)')
ax1.set_xlabel('x, mm'); ax1.set_ylabel('p, kPa')
ax1.set_xlim(0,250); ax1.set_ylim(0,2.5)
ax1.legend(loc='upper right')
ax1.set_title('')
if len(sys.argv) == 2 and sys.argv[1] == "--save-figure":
    fig.savefig("convex-ramp-surf-conds.png", transparent=True, dpi=300)
else:
    plt.show()
