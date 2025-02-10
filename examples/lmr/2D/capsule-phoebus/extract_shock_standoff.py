"""
Python script for shock shape based grid tailoring for blunt bodies, lmr version

@author: Nick Gibbons
"""

from numpy import array, concatenate, argmin, amax, arange, mean, arctan2, linspace, pi, cos, sin, degrees, sqrt, inf, zeros, diff
from glob import glob
#from pylab import plot,show
import json
import gzip
import yaml

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

def read_grid_block(filename):
    with gzip.open(filename, mode='rt') as fp:
        lines = [line.strip() for line in fp]

    lineiter = iter(lines)

    grid = {}
    grid['gtype'] = next(lineiter)
    grid['label'] = next(lineiter).split(':')[-1]
    grid['dimensions'] = int(next(lineiter).split(':')[-1])
    grid['niv'] = int(next(lineiter).split(':')[-1])
    grid['njv'] = int(next(lineiter).split(':')[-1])
    grid['nkv'] = int(next(lineiter).split(':')[-1])

    npoints = grid['niv']*grid['njv']*grid['nkv']
    points = []
    for i in range(npoints):
        line = next(lineiter)
        point_string = line.split()
        point_numbers = list(map(number, point_string))
        points.append(point_numbers)

    grid['points'] = array(points)
    return grid

def to_2D_structured_grid(grid):
    nx = grid['niv']
    ny = grid['njv']
    newgrid = {k:v for k,v in grid.items() if k!='points'}
    newgrid['points'] = grid['points'].reshape(ny,nx,3)
    return newgrid

def plot_block(g):
    plot(g[0,:,0],  g[0,:,1],  'b-')
    plot(g[-1,:,0], g[-1,:,1], 'g-')
    plot(g[:,0,0],  g[:,0,1],  'r-')
    plot(g[:,-1,0], g[:,-1,1], 'k-')

def plot_grid(g):
    plot(g['points'][:,0,0],  g['points'][:,0,1],  'r-')
    plot(g['points'][:,-1,0], g['points'][:,-1,1], 'k-')
    plot(g['points'][0,:,0],  g['points'][0,:,1],  'b-')
    plot(g['points'][-1,:,0], g['points'][-1,:,1], 'g-')

def outline_grid(g):
    for i in range(len(g['points'])):
        plot(g['points'][i,:,0],  g['points'][i,:,1],  'k.')

def plot_block_points(g, color, marker):
    plot(g[0,:,0],  g[0,:,1],  color=color, linestyle='none', marker=marker)
    plot(g[-1,:,0], g[-1,:,1], color=color, linestyle='none', marker=marker)
    plot(g[:,0,0],  g[:,0,1],  color=color, linestyle='none', marker=marker)
    plot(g[:,-1,0], g[:,-1,1], color=color, linestyle='none', marker=marker)

def write_shock_shape_to_eilmer_file(x,y, filename):
    with open(filename, 'w') as fp:
        for xi,yi in zip(x,y):
            fp.write('{:.7e} {:.7e} {:.7e}\n'.format(xi, yi, 0.0))
    return

def concatenate_grids(grids, axis):
     assert len(grids)>1
     trimgrids = []

     colors = ['red','blue','green','black','magenta','darkgrey','lightgrey','yellow']
     for i,g in enumerate(grids):
         if i==0:
             start=0
         else:
             start=1

         if axis==0:
             #print("    ", g[0,0,0:2], g[-1,0,0:2], g.shape)
             trimgrids.append(g[start:,:,:].copy())
             #plot_block(g[start:,:,:].copy() )
         if axis==1:
             #print("    ", g[0,0,0:2], g[0,-1,0:2], g.shape)
             trimgrids.append(g[:,start:,:].copy())
             #plot_block(g[:,start:,:].copy())

     return concatenate(trimgrids, axis)

def concatenate_fluid_block_array_grids(grids, nib, njb):
    hgrids = []
    for j in range(njb):
        hgrid = []
        for i in range(nib):
            idx = i*njb + j
            hgrid.append(grids[idx]['points'])
        hgridcat = concatenate_grids(hgrid, axis=1)
        hgrids.append(hgridcat)
    gridarray = concatenate_grids(hgrids, axis=0)
    return gridarray

def get_fluid_block_array_grid():
    """

    """
    snapshot_dirs = sorted(glob('lmrsim/snapshots/0*'))
    snapshot_idxs = [i.split('/')[-1] for i in snapshot_dirs]
    snapshot_idx = snapshot_idxs[-1]

    gridfilenames = sorted(glob('lmrsim/snapshots/{}/grid-*.gz'.format(snapshot_idx)))
    grids = [to_2D_structured_grid(read_grid_block(name)) for name in gridfilenames]

    config = read_config_file('lmrsim/config')
    nib = config['fluid_block_array_0']['nib']
    njb = config['fluid_block_array_0']['njb']
    g = concatenate_fluid_block_array_grids(grids, nib=nib, njb=njb)
    return g

def get_shock_standoff():
    g = get_fluid_block_array_grid()
    ss = -1.0*g[0,0,0]
    return ss

if __name__=='__main__':
    ss = get_shock_standoff() 
    print("Shock Standoff: ", ss*1000.0, "mm")
