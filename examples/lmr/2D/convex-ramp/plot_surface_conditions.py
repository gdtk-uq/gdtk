# plot_surface_conditions.py
# Assemble the ramp loads by picking data out of the loads files.
# Peter J. 2024-03-13
#
import pandas as pd
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

import numpy as np
xp = np.loadtxt('mohammadian-figure-12-p_p_inf.data')
xq = np.loadtxt('mohammadian-figure-13-heat-flux.data')

import pylab as pl
pl.subplot(2,1,1)
pl.plot(ramp_df['pos.x']*1000, ramp_df['q_total']/1000, '-b', label='Eilmer')
pl.plot(xq[:,0]*25.4, xq[:,1]*11.5, '+k', label='Mohammadian (1972)')
pl.xlabel(''); pl.ylabel('q, kW/m^2')
pl.xlim(0,250); pl.ylim(0,120)
pl.legend(loc='upper right')
pl.title('Surface conditions on convex ramp')
pl.subplot(2,1,2)
pl.plot(ramp_df['pos.x']*1000, ramp_df['p']/1000, '-b', label='Eilmer')
pl.plot(xp[:,0]*25.4, xp[:,1]*66.43/1000, '+k', label='Mohammadian (1972)')
pl.xlabel('x, mm'); pl.ylabel('p, kPa')
pl.xlim(0,240); pl.ylim(0,2.5)
pl.legend(loc='upper right')
pl.title('')
pl.show()

