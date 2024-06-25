"""
Compute wall drag from a set of loads files

@author: Nick Gibbons
"""

from numpy import array, argsort, concatenate, sqrt
from sys import argv
from glob import glob

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
        lines = [line.strip().split() for line in fp]

    header = lines[0]
    body = [list(map(number,line)) for line in lines[1:]]
            
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

    idxs = argsort(data[sortkey])
    for k,v in data.items():
        data[k] = v[idxs].copy()

    return data


def wall_data(groupname='wall', time=-1):
    times = read_loads_file('lmrsim/loads/loads-metadata')
    finalindex = times['loads_index'][time]
    finalindex_string = str(finalindex).zfill(4)

    datafilepattern = 'lmrsim/loads/{ftime}/blk-*-{gname}.dat'.format(ftime=finalindex_string, gname=groupname)
    datafilenames = glob(datafilepattern)
    datafiles = [read_loads_file(filename) for filename in datafilenames]

    data = collate_datafiles(datafiles)
    return data

def read_aoa(filename):
    with open(filename) as fp:
        contents = fp.read()
    aoa = float(contents.strip().split('=')[-1])
    return aoa

def read_com(filename):
    with open(filename) as fp:
        contents = fp.read()
    com = [float(i) for i in contents.strip().split(',')]
    return com

if __name__=='__main__':
    data = wall_data()
    aoa = read_aoa('angle_of_attack.lua')
    com = read_com('centre_of_mass.lua')

    p = data['p']
    nx = data['n.x']
    ny = data['n.y']
    area = data['area']
    outsign = data['outsign']
    xforce = (p*area*nx*outsign).sum()
    yforce = (p*area*ny*outsign).sum()
    taux = (data['tau_wall.x']*data['area']*data['outsign']).sum()
    tauy = (data['tau_wall.y']*data['area']*data['outsign']).sum()
    # Compute pitching moment 
    Fx = (p*nx + data['tau_wall.x'])*area*outsign
    Fy = (p*ny + data['tau_wall.y'])*area*outsign
    dx = data['pos.x'] - com[0]
    dy = data['pos.y'] - com[1]
    M = (dx*Fy - dy*Fx).sum()

    print("Alpha: ", aoa)
    print("Inviscid x-Force: {} (N)".format(xforce))
    print("Inviscid y-Force: {} (N)".format(yforce))
    print("Viscous x-Force: {} (N)".format(taux))
    print("Viscous y-Force: {} (N)".format(tauy))
    print("Pitching Moment: {} (N.m)".format(M))

