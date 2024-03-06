"""
Integration test script for multi-temperature sphere-nonaka case

@author: Nick Gibbons
"""

from math import sqrt
#from pylab import plot,show

def interpolate(x,y,xtarget):
    assert len(x)==len(y)
    n = len(x)

    for l,u in zip(range(0,n-1), range(1,n)):
        xl = x[l]
        xu = x[u]

        if (xl<xtarget) and (xtarget<xu):
            yl = y[l]
            yu = y[u]
            ytarget = (yu-yl)/(xu-xl)*(xtarget-xl) + yl
            return ytarget
    else:
        raise Exception("No matching pair found for interpolation of x: {}".format(xtarget))

def read_dat_file(filename):
    with open(filename) as fp:
        data = fp.read()

    header = None
    body = []

    for line in data.splitlines():
        if line.startswith('#'):
            header = line.strip('#').strip().split()
        else:
            body.append(list(map(float,line.strip().split())))

    cols = [list(i) for i in zip(*body)]
    return header, cols


_, explower = read_dat_file("expt-data/fig10-lower.g3data")
exptheta, expstandoff = explower

cfd = dict(zip(*read_dat_file("shock-shape.dat")))
cfdstandoff = [interpolate(cfd['theta'], cfd['d/R'], theta) for theta in exptheta]

rms = sqrt(1.0/len(cfdstandoff)*sum((exps-cfds)**2 for exps,cfds in zip(expstandoff, cfdstandoff)))
print("Shock Standoff RMS: ", rms)

#plot(cfd['theta'], cfd['d/R'])
#plot(exptheta, expstandoff, 'ro')
#plot(exptheta, cfdstandoff, 'bo')
#show()
