"""
Python script for shock shape based grid tailoring for blunt bodies

Note: 

@author: Nick Gibbons
"""

from numpy import array, concatenate, argmin, amax, arange, mean, arctan2, linspace, pi, cos, sin, degrees, sqrt, inf, zeros, diff
from glob import glob
import gzip
from pylab import plot,show
from numpy.polynomial.polynomial import Polynomial
from libe4post import read_config_file, number
from sys import argv
from configparser import ConfigParser

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

    values = [list(map(number, line.split())) for line in lineiter]
    grid['points'] = array(values)
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

def get_fluid_block_array_grid(directory_name, tidx=-1):
    job_name = glob(directory_name+'/config/*.config')[0].split('/')[-1].split('.config')[0]
    config = read_config_file(directory_name+'/config/{}.config'.format(job_name))
    times = sorted(glob(directory_name+'/grid/t*'))
    if len(times)==0: raise Exception("No grids found!")
    timedir = times[tidx]
    print("Picked up grid directory at time: ", timedir)

    gridfilenames = sorted(glob('{}/*.gz'.format(timedir)))

    grids = [to_2D_structured_grid(read_grid_block(name)) for name in gridfilenames]
    #for g in grids: plot_grid(g)

    # Block indexes are ordered y first fsr. 
    nib = config['fluid_block_array_0']['nib']
    njb = config['fluid_block_array_0']['njb']
    g = concatenate_fluid_block_array_grids(grids, nib=nib, njb=njb)
    return g

def write_gpdat_file(filename, points):
    """ Write collection of points to a GridPro tube file """
    print("in write_dat_file, fname: {} points.shape: {}".format(filename, points.shape))
    with open(filename,'w') as fp:
        fp.write('{}\n'.format(len(points)))

        for p in points:
            fp.write("{} {} {}\n".format(*p))
    return

def displace_shock_slightly(points, displacement):
    xs,ys = points[:,0], points[:,1]
    dxdy = zeros(xs.shape)
    dxdy[1:] = diff(xs)/diff(ys)
    dxdy[0] = dxdy[1]

    Tx = dxdy.copy()
    Ty = 0.0*Tx + 1.0
    Tmag = sqrt(Tx**2 + Ty**2)
    Tx/=Tmag
    Ty/=Tmag

    N = xs.size
    Nx = -1.0*Ty.copy()
    Ny =  1.0*Tx.copy()
    d0 = displacement
    d1 = 5*displacement
    a = d1 - d0
    d = d0 + a*(ys/ys.max())**2

    outpoints = zeros(points.shape)
    outpoints[:,0] = xs + d*Nx
    outpoints[:,1] = ys + d*Ny

    return outpoints
   
def stretch_last_vertex_to_x_position(points, xout):
    """ No need to add more points, just stretch the last cell a little"""
    xb,yb = points[:,0], points[:,1]
    Tx = xb[-1] - xb[-2]
    Ty = yb[-1] - yb[-2]
    d = sqrt(Tx**2 + Ty**2) # d is the spacing of the last two points
    Tx /= d
    Ty /= d

    yout = yb[-1] + (xout - xb[-1])/Tx*Ty
    D = sqrt((yout-yb[-1])**2 + (xout-xb[-1])**2) # D is the length of boundary we need to interpolate over

    newpoints = points.copy()
    newpoints[-1,0] = xout
    newpoints[-1,1] = yout
    return newpoints


if __name__=='__main__':
    displacement = float(argv[1])

    g = get_fluid_block_array_grid('.')
    points = g[:,0,:].copy()
    dpoints = displace_shock_slightly(points, displacement)
    epoints = stretch_last_vertex_to_x_position(dpoints, points[-1,0])
    write_shock_shape_to_eilmer_file(epoints[:,0], epoints[:,1], 'inflated_shock_shape.dat')
    #write_gpdat_file('inflated_shock_shape2.gpdat', epoints)

    plot_block(g)
    plot(points[:,0], points[:,1], 'ro')
    plot(points[:,0], points[:,1], 'r-')
    plot(dpoints[:,0], dpoints[:,1], 'b-')
    plot(epoints[:,0], epoints[:,1], 'g.')
    show()
