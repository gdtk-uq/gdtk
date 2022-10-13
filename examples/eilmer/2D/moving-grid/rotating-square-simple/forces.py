#! /bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sim_data = pd.read_csv('square_data.csv', delimiter='\ ', engine="python")

t = sim_data['time(s)'].to_numpy()
theta = sim_data['thetab'].to_numpy()
omega = sim_data['omegab'].to_numpy()
Fx    = sim_data['Fxb'].to_numpy()
Fy    = sim_data['Fyb'].to_numpy()
Mz    = sim_data['Mzbg'].to_numpy()

def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

fig1,ax1 = plt.subplots()
ax1.set_xlabel('$t$ (ms)')
ax1.set_ylabel('Force (N/m)')
ax1.grid()
ax1.set_xticks([0,0.5,1.0,1.5,2.0])
ax1.set_xlim([-0.1, 2.1])
ax1.plot(t*1000,Fx,linewidth=0.7,color='black',label='Fx')
ax1.plot(t*1000,Fy,'--',linewidth=0.7,color='black',label='Fy')
ax2 = ax1.twiny()
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels([int(180* theta[find_nearest(t, val/1000)] / np.pi) for val in ax1.get_xticks()])
ax2.set_xlabel(r'$\theta ^\circ$ ')
ax1.legend(loc='center')
#fig1.savefig('./forces.svg')


fig2,ax3 = plt.subplots()
ax3.set_xlabel('$t$ (ms)')
ax3.set_xticks([0,0.5,1.0,1.5,2.0])
ax3.set_ylabel('Moment (Nm/m)')
ax3.set_xlim([-0.1, 2.1])
ax3.grid()
ax3.plot(t*1000,Mz,linewidth=0.7,color='black')
ax4 = ax3.twiny()
ax4.set_xticks(ax1.get_xticks())
ax4.set_xbound(ax1.get_xbound())
ax4.set_xticklabels([int(180* theta[find_nearest(t, val/1000)] / np.pi) for val in ax1.get_xticks()])
ax4.set_xlabel(r'$\theta ^\circ$ ')
#fig2.savefig('./moment.svg')

plt.show()
