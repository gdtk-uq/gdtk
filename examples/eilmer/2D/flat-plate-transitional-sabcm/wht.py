"""
Plot the heat transfer to the plate.

@author: Nick Gibbons
"""

from glob import glob
import json
from numpy import array, concatenate, argsort, stack
import matplotlib.pyplot as plt

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber

def read_loads_file(filename):
    header = None
    body = []
    with open(filename) as fp:
        for line in fp:
            line = line.strip().split()
            if line[0].startswith('#'):
                header = line[1:]
            else:
                body.append(list(map(number,line)))
            
    assert header!=None
    header = [i.split(':')[-1] if ':' in i else i for i in header]
    cols = [array(col) for col in zip(*body)]
    data = dict(zip(header,cols))
    return data

def collate_datafiles(datafiles, sortkey='pos.x'):
    keys = set()
    for f in datafiles: keys = keys.union(set(f.keys())) 

    data = {}
    for k in keys:
        values = [f[k] for f in datafiles]
        data[k] = concatenate(values)

    if not sortkey is None:
        idxs = argsort(data[sortkey])
        for k,v in data.items():
            data[k] = v[idxs].copy()

    return data

def read_config_file(filename):
    with open(filename) as fp:
        config = json.load(fp)
    return config

plt.rcParams['svg.fonttype'] = 'none'

jobname = glob('./config/*.config')[0].split('/')[-1].split('.config')[0]
metadata = read_loads_file('./loads/{}-loads.times'.format(jobname))

index = metadata['loads_index'][-1]
time_string = str(index).zfill(4)
datafilenames = sorted(glob('./loads/t{}/*.wall.dat'.format(time_string)))
datafiles = [read_loads_file(filename) for filename in datafilenames]
data = collate_datafiles(datafiles, sortkey='pos.x')

config = read_config_file('config/{}.config'.format(jobname))
Tu = config['freestream_turbulent_intensity']
niv = config['fluid_block_array_0']['niv']-1
njv = config['fluid_block_array_0']['njv']-1

fig = plt.figure(figsize=(7,5))
ax = fig.subplots(1,1) 

# First entry in He and Morgan, table 1
Reu = 4.99e6
Ret = 1.28e6
Ret_uncertainty = 0.15e6
xt_upper = (Ret+Ret_uncertainty/2)/Reu
xt_lower = (Ret-Ret_uncertainty/2)/Reu

ax.plot(data['pos.x'], -1.0*data['q_total']/1e4, 'k-', linewidth=2, label='Mach=6.52, Re=4.99e6/m')
ax.plot([xt_lower, xt_lower], [0,30], 'r--', linewidth=1, label="Experimental Transition Location")
ax.plot([xt_upper, xt_upper], [0,30], 'r--', linewidth=1)

ax.set_title("Heat Transfer Rate, Turbulent Intensity={:.1f}%, ({:d}x{:d})".format(Tu*100, niv, njv))
ax.set_xlabel('x (mm)')
ax.set_ylabel('q (W/cm2)')
ax.set_ylim(0,100)
ax.grid()
ax.legend(framealpha=1.0)
plt.tight_layout()
plt.savefig('wht_{}x{}_{:.1f}.svg'.format(niv,njv,Tu*100))
#plt.close()
plt.show()
