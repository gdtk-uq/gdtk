"""
Convert a p2d file to eilmer grid format

@author: Nick Gibbons
"""

from numpy import array, zeros, absolute, argmax, argmin
from pylab import plot,show
from sys import argv

def read_p2d(filename):
    """ Read in a single block p2d file """
    with open(filename) as fp:
        lines = [line.strip().split() for line in fp]  
    
    nblocks = int(lines[0][0])
    nx, ny = int(lines[1][0]), int(lines[1][1])
    data = []
    for line in lines[2:]:
        data.extend(list(map(float, line)))
    data = array(data)
    print("data.shape", data.shape)
    print("nx: {} ny: {}".format(nx, ny))
    npoints = nx*ny
    print("npoints", npoints)
    print("data.size", data.size)
    assert data.size==npoints*2
    
    x = data[:npoints]
    y = data[npoints:]
    x = x.reshape((ny,nx))
    y = y.reshape((ny,nx))
    gridpoints = zeros((ny,nx,2))
    gridpoints[:,:,0] = x
    gridpoints[:,:,1] = y
    return gridpoints 


def split_blocks_at_x(grid, xsplit):
    """ Divide a grid into two blocks at the nearest point to x """
    startnx = argmin(absolute(grid[0,:,0]-xsplit))
    print("startnx", startnx)
    block0 = grid[:,:startnx+1,:].copy()
    block1 = grid[:,startnx:,:].copy()
    return block0, block1

def write_eilmer3_sblock(block, filename):
    """ Write a block of points to an old eilmer3 format file """
    # Strangely, eilmer4 can only read in eilmer3 format during prep
    with open(filename,'w') as fp:
        ny,nx,dim = block.shape
        fp.write('{} {} {}\n'.format(nx, ny, 1))
        for j in range(ny):
            for i in range(nx):
                fp.write('{:1.18e} '.format(block[j,i,0]))
                fp.write('{:1.18e} '.format(block[j,i,1]))
                fp.write('{:1.18e}\n'.format(0.0))
    return


#filename = 'flatplate_clust2.p2dfmt'
#filename = 'flatplate_clust2_3levelsdown_69x49.p2dfmt'
filename = argv[1]
grid = read_p2d(filename)

xsplit = 0.0
block0, block1 = split_blocks_at_x(grid, xsplit)
#plot(block0[:,:,0], block0[:,:,1], 'r+')
#plot(block1[:,:,0], block1[:,:,1], 'b.')
#show()

write_eilmer3_sblock(block0, 'larc.b0000.txt')
write_eilmer3_sblock(block1, 'larc.b0001.txt')
print("Done")

