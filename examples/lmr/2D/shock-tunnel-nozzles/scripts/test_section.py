"""
 References:
   "Experimental validation of the T4 Mach 7.0 nozzle"
   W. Y. K. Chan, M. K. Smart and P. A. Jacobs
   School of Mechanical & Mining Engineering
   Research Report Research Report Number 2014/14
   https://espace.library.uq.edu.au/view/UQ:378568

@author: Nick Gibbons
"""

from numpy import argmax, array
from gdtk.lmr import LmrConfig, SimInfo
import matplotlib.pyplot as plt
from sys import argv
from os import getcwd, chdir

plt.rcParams.update({'font.size': 12})
plt.rcParams['svg.fonttype'] = 'none'

pstag = 19.33e6
target_x = 1.0 + 141e-3 # mm

pitot_data = [
    [0.002883341823739175,0.010600291179164265],
    [0.01229750382068263,0.010499358348017558],
    [0.027030056036678556,0.01081227295123723],
    [0.03295975547631176,0.01057021565630837],
    [0.035038206826286306,0.010504074002934058],
    [0.056984207845134985,0.010858733565497214],
    [0.0629750382068263,0.010572943721136097],
    [0.0651146204788589,0.010670719791553573],
    [0.06688741721854306,0.010768596076597376],
    [0.06688741721854306,0.010604661650367869],
    [0.0720224146714213,0.010821837880571998],
    [0.0720224146714213,0.01065790345434249],
    [0.0863881813550688,0.010435399146505427],
    [0.08638818135506877,0.010686765266724008],
    [0.0870606214977076,0.010686581539909078],
    [0.09054508405501782,0.01058726884522129],
    [0.09048395313295973,0.010412422159680867],
    [0.09292919001528276,0.010586617450150179],
    [0.09500764136525723,0.010706268146502645],
    [0.09500764136525729,0.010902989457978053],
    [0.11077941925624046,0.009882286786423139],
    [0.11401935812531838,0.008438778606404265],
    [0.11401935812531838,0.008045335983453448],
    [0.11701477330616405,0.008164736143240101],
    [0.12300560366785535,0.006709547391775157],
    [0.12508405501782988,0.006249963115450029],
    [0.12502292409577181,0.005605171074718352],
    [0.13737137035150282,0.0026291195866703365],
    [0.1373713703515028,0.002268463848965419],
    [0.14000000000000004,0.001688510671465828],
    [0.14,0.001448073512995883],
]
exp_y = array([p[0] for p in pitot_data])
pp_on_pe = array([p[1] for p in pitot_data])

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

    xs = ss.fields[8]['pos.x'][:,0]
    idx = argmax(xs>target_x)
    print("x pos {} @ idx {}".format(xs[idx], idx))
    data8 = {key:val[idx, :].copy() for key,val in ss.fields[8].items()}

    xs = ss.fields[11]['pos.x'][:,0]
    idx = argmax(xs>target_x)
    print("x pos {} @ idx {}".format(xs[idx], idx))
    data11 = {key:val[idx, :].copy() for key,val in ss.fields[11].items()}

    data ={}
    for key in data8.keys():
        data[key] = array(list(data8[key]) + list(data11[key]))
    data['pitot_p'] = 0.93872*data['rho']*data['vel.x']**2

    datas.append(data)
    chdir(savedir)

fig = plt.figure(figsize=(7,4))
ax0 = fig.subplots(1,1)

fig.suptitle("Pitot Pressure/Stagnation Pressure @ x=141mm")

linestyles = ['-', '--', ':', '-.']
for i,d in enumerate(datas):
    ax0.plot(d['pos.y'], d['pitot_p']/pstag,  color='black', linestyle=linestyles[i], label=directories[i], linewidth=2.0)
    ax0.errorbar(exp_y, pp_on_pe, yerr=pp_on_pe*0.1,  color='black', marker='o', capsize=1, elinewidth=1, markeredgewidth=1, linestyle="None")

ax0.set_ylabel('p/pp')
ax0.set_xlabel('y (mm)')

ax0.set_xlim(0.0, 0.16)
ax0.set_ylim(0.0, 0.014)
ax0.legend(framealpha=1.0)
ax0.grid()
fig.tight_layout()
plt.show()
