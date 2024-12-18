"""
Abridged code for reading legacy eilmer files directly from flow

@author: Nick Gibbons
"""
from numpy import array, concatenate
from glob import glob
import yaml
import json
import gzip

def read_config_file(filename):
    with open(filename) as fp:
        config = json.load(fp)
    return config

def read_flow_block_gz(filename, ncells, variables):
    with gzip.open(filename, mode='rt') as fp:
        rawdata = array([float(line.strip()) for line in fp])

    data = {}
    for i,k in enumerate(variables):
        data[k] = rawdata[i*ncells:(i+1)*ncells].copy()

    return data
    

def get_fluid_block_array(dir, tidx=-1):
    config = read_config_file(dir+'/lmrsim/config')

    with open(dir+'/lmrsim/snapshots/snapshot-times-metadata') as fp:
       sstfmd = fp.read()
    times = yaml.safe_load(sstfmd)
    times_idx = [k for k in times.keys()][tidx]

    with open(dir+'/lmrsim/snapshots/fluid.metadata') as fp:
       fluidmd_yaml = fp.read()
    fluidmd = yaml.safe_load(fluidmd_yaml)

    blks = []
    fba = config['fluid_block_array_0']
    for j in range(fba['njb']):
        for i in range(fba['nib']):
            id = fba['idarray'][i][j]
            ncells = fba['njcs'][j]*fba['nics'][i]

            flowfilename = dir+'/lmrsim/snapshots/{}/fluid-{}.gz'.format(times_idx,str(id).zfill(4))
            blkdata = read_flow_block_gz(flowfilename, ncells, fluidmd['variables'])
            sblkdata = to_2D_structured_block(blkdata, fba['nics'][i], fba['njcs'][j])
            blks.append(sblkdata) 

    f = concatenate_fluid_block_array(blks, nib=fba['nib'], njb=fba['njb'])
    return f

def to_2D_structured_block(data, nx, ny):
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

if __name__=='__main__':
    blks = get_fluid_block_array('.')
    print("blks[T].shape", blks['T'].shape)
