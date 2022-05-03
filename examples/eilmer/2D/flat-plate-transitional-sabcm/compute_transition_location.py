"""
Compute wall drag from a set of loads files

@author: Nick Gibbons
"""

from numpy import array, argsort, concatenate, sqrt, diff, argmin, argmax
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

    idxs = argsort(data[sortkey])
    for k,v in data.items():
        data[k] = v[idxs].copy()

    return data


def wall_data(jobname, groupname='wall', time=-1):
    times = read_loads_file('loads/{}-loads.times'.format(jobname))
    finalindex = times['loads_index'][time]
    finalindex_string = str(finalindex).zfill(4)

    datafilepattern = 'loads/t{ftime}/*.t{ftime}.{gname}.dat'.format(ftime=finalindex_string, gname=groupname)
    datafilenames = glob(datafilepattern)
    datafiles = [read_loads_file(filename) for filename in datafilenames]

    data = collate_datafiles(datafiles)
    return data

def transition_location(data):
    # Heat transfer out of the domain is defined as negative, so we flip it first off
    q = data['q_total']*-1.0
    diffs = diff(q)
    idx_t = argmax(diffs)
    x_t = data['pos.x'][idx_t]
    return x_t

if __name__=='__main__':
    jobname = argv[1]
    data = wall_data(jobname)
    x_t = transition_location(data)
    print("Transition Location: {} (m)".format(x_t))

