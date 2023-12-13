"""
Abridged code for reading legacy eilmer files directly from flow

@author: Nick Gibbons
"""
from numpy import array, concatenate
from glob import glob
import json
import gzip

def read_config_file(filename):
    with open(filename) as fp:
        config = json.load(fp)
    return config

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber


def get_fluid_block_array(directory_name, tidx=-1):
    job_name = glob(directory_name+'/config/*.config')[0].split('/')[-1].split('.config')[0]
    config = read_config_file(directory_name+'/config/{}.config'.format(job_name))
    times = sorted(glob(directory_name+'/flow/t*'))
    if len(times)==0: raise Exception("No flow solutions found!")
    timedir = times[tidx]
    print("Picked up flow solution at time: ", timedir)

    flowfilenames = sorted(glob('{}/*.gz'.format(timedir)))

    blocks = [to_2D_structured_block(*read_flow_block(name)) for name in flowfilenames]

    # Block indexes are ordered y first fsr. 
    nib = len(blocks)
    njb = 1
    f = concatenate_fluid_block_array(blocks, nib=nib, njb=njb)
    return f

def read_flow_block(filename):
    with gzip.open(filename, mode='rt') as fp:
        lines = [line.strip() for line in fp]

    lineiter = iter(lines)

    metadata = {}
    metadata['flow_type'] = next(lineiter)
    metadata['label'] = next(lineiter)
    metadata['sim_time'] = float(next(lineiter).split(':')[-1])
    metadata['variables'] = int(next(lineiter).split(':')[-1])
    metadata['var_names'] = [i.strip('"') for i in next(lineiter).split()]
    metadata['dimensions'] = int(next(lineiter).split(':')[-1])
    metadata['nicell'] = int(next(lineiter).split(':')[-1])
    metadata['njcell'] = int(next(lineiter).split(':')[-1])
    metadata['nkcell'] = int(next(lineiter).split(':')[-1])

    values = [list(map(number, line.split())) for line in lineiter]
    cols = [array(col) for col in list(zip(*values))]
    assert len(cols)==len(metadata['var_names'])
    data = dict(zip(metadata['var_names'],cols))
    return metadata, data 

def to_2D_structured_block(metadata, data):
    nx = metadata['nicell']
    ny = metadata['njcell']
    newdata = {}
    for k,v in data.items():
        newdata[k] = v.reshape(ny,nx).copy()
    return newdata

def concatenate_dictionaries(dicts, axis):
     keys = dicts[0].keys()
     return {k:concatenate([d[k] for d in dicts], axis) for k in keys} 

def concatenate_fluid_block_array(blocks, nib, njb):
    hblocks = []
    for j in range(njb):
        hblock = []
        for i in range(nib):
            idx = i*njb + j
            hblock.append(blocks[idx])
        hblock_dict = concatenate_dictionaries(hblock, axis=1)
        hblocks.append(hblock_dict)
    blockarray = concatenate_dictionaries(hblocks, axis=0)
    return blockarray
