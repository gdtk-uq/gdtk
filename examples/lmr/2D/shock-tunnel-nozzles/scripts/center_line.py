"""
Python code for plotting with PJ's new gdtk.lmr interface

@author: Nick Gibbons
"""

from gdtk.lmr import LmrConfig, SimInfo
import matplotlib.pyplot as plt
from sys import argv
from os import getcwd, chdir
from numpy import concatenate

plt.rcParams.update({'font.size': 12})
plt.rcParams['svg.fonttype'] = 'none'

if len(argv) == 1:
    directories = ['.']
else:
    directories = argv[1:]
datas = []
for directory in directories:
    savedir = getcwd()
    print("directory: ", directory)
    print("savedir: ", savedir)
    chdir(directory)

    lmrcfg = LmrConfig()
    sim = SimInfo(lmrcfg)
    ss = sim.read_snapshot(sim.snapshots[-1])
    data = {key:list() for key in ss.fields[0].keys()}
    keys = list(ss.fields[0].keys())
    outsidex = list()
    outsidey = list()

    # I've built the blocks weird, with the expanding section first and the throat last
    blks = [len(ss.fields)-1] + [n for n in range(len(ss.fields)-1)]
    for n in blks:
        outsidex.append(ss.fields[n]['pos.x'][:,-1])
        outsidey.append(ss.fields[n]['pos.y'][:,-1])
        for key,value in ss.fields[n].items():
            data[key].append(value[:,0].copy())
    data = {key:concatenate(val) for key,val in data.items()}
    outsidex = concatenate(outsidex)
    outsidey = concatenate(outsidey)
    datas.append(data)
    chdir(savedir)

fig = plt.figure(figsize=(9,7))
ax0,ax1,ax2 = fig.subplots(3,1, gridspec_kw={'height_ratios': [1, 1, 0.5]})

fig.suptitle("Nozzle Centreline Conditions")

linestyles = ['-', '--', ':', '-.']
for i,d in enumerate(datas):
    ax0.plot(d['pos.x']*1000.0, d['T'], color='red', linestyle=linestyles[i], label=directories[i])
    ax1.semilogy(d['pos.x']*1000.0, d['massf-O'], color='mediumblue', linestyle=linestyles[i], label=directories[i])
    ax2.plot(outsidex*1000.0, outsidey*1000.0, color='black', linestyle=linestyles[i], linewidth=2.0)
    ax2.plot(data['pos.x']*1000.0, data['pos.y']*1000.0, color='black', linestyle='--', linewidth=1.0)

ax2.set_xlabel('x (mm)')

ax0.set_ylabel('T (K)')
ax1.set_ylabel('Y[O]')
ax2.set_ylabel('y (mm)')

ax0.set_xticklabels([])
ax1.set_xticklabels([])


ax0.legend(framealpha=1.0)
ax1.legend(framealpha=1.0)
ax0.grid()
ax1.grid()
fig.tight_layout()
plt.savefig('five_six_centre_line.svg')
plt.show()
