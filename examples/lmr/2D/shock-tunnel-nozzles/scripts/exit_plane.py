"""
Python code for plotting with PJ's new gdtk.lmr interface

@author: Nick Gibbons
"""

from gdtk.lmr import LmrConfig, SimInfo
import matplotlib.pyplot as plt
from sys import argv
from os import getcwd, chdir

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
    print(list(ss.fields[0].keys()))

    xs = ss.fields[14]['pos.x'][-1,:]*1000.0
    print("x pos; ", xs[0])
    data = {key:val[-1,:].copy() for key,val in ss.fields[14].items()}
    data['pitot_p'] = 0.93872*data['rho']*data['vel.x']**2
    datas.append(data)
    chdir(savedir)

fig = plt.figure(figsize=(13,5))
ax0,ax1,ax2,ax3 = fig.subplots(1,4)

fig.suptitle("Nozzle Outflow Conditions")

linestyles = ['-', '--', ':', '-.']
for i,d in enumerate(datas):
    ax0.plot(d['T'], d['pos.y']*1000.0, color='red', linestyle=linestyles[i], label=directories[i])
    ax1.plot(d['p'], d['pos.y']*1000.0, color='blue', linestyle=linestyles[i])
    ax2.plot(d['vel.x'], d['pos.y']*1000.0, color='black', linestyle=linestyles[i])
    ax3.plot(d['pitot_p']/19.33e6, d['pos.y']*1000.0, color='darkviolet', linestyle=linestyles[i])

ax0.set_ylabel('y (mm)')
ax0.set_xlabel('T (K)')
ax1.set_xlabel('p (Pa)')
ax2.set_xlabel('vel.x (m/s)')
ax3.set_xlabel('pp/ps')

ax0.set_title('Temperature')
ax1.set_title('Pressure')
ax2.set_title('Velocity')
ax3.set_title('Pitot/Stag Pressure')

ax0.legend(framealpha=1.0)
ax0.grid()
ax1.grid()
ax2.grid()
ax3.grid()
fig.tight_layout()
plt.savefig('exit_plane.svg')
plt.show()
