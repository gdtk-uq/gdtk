# plot_surface_conditions.py
# Assemble the body loads by picking data out of the loads files.
#
import pandas as pd
#
# We know that the ramp surface consists of the eastern boundaries
# of blocks 2 and 3 plus the southern boundary of block 4.
blkbndy_ids = [(2,1), (3,1), (4,2)]
#
# For the final snapshot, assemble the data for southern boundaries
# into a single data frame.
meta_df = pd.read_table('./lmrsim/loads/loads-metadata', sep='\s+')
snaps = meta_df['loads_index']
final_snap = snaps[snaps.index[-1]]
blk_dfs = []
for b_b in blkbndy_ids:
    # print("b_b=", b_b)
    f = './lmrsim/loads/%04d/blk-%04d-bndry-%d-loads.dat' % (final_snap, b_b[0], b_b[1])
    blk_dfs.append(pd.read_table(f, sep='\s+'))
body_df = pd.concat(blk_dfs)
# print(body_df)
# print(body_df.columns)
drag = sum(body_df['p']*body_df['area']*body_df['n.x']*body_df['outsign'])
print('drag=', drag)

import numpy as np

import pylab as pl
# pl.subplot(1,1,1)
pl.plot(body_df['pos.x']*1000, body_df['p']/1000, '-b')
pl.xlabel('x'); pl.ylabel('p, kPa')
pl.xlim(0,17.0); pl.ylim(0,70)
# pl.legend(loc='upper right')
pl.title('Surface pressure')
pl.show()

