"""
Plot residual files from an e4-nk-shared simulation

@author: Nick Gibbons
"""

from numpy import log10, array
from pylab import semilogy, show, legend, xlabel

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber

def read_diagnostics_file(filename):
    with open(filename) as fp:
        lines = [line.strip() for line in fp]

    keys = []
    values = []

    for line in lines:
        if line.startswith('#'):
            keys.append(line)
        else:
            values.append(line)

    keys = [k.split(':')[-1].strip() for k in keys]
    values = [array(list(map(number, line.split()))) for line in values]
    cols = list(zip(*values))
    assert len(cols)==len(keys)
    data = dict(zip(keys,cols))
    return data

if __name__=='__main__':
    data = read_diagnostics_file('e4-nk.diagnostics.dat')
    semilogy(data['step'], data['global-residual-rel'], label='global-residual-rel')
    semilogy(data['step'], data['mass-rel']  , label = 'mass-rel')
    semilogy(data['step'], data['y-mom-rel'] , label = 'y-mom-rel')
    semilogy(data['step'], data['x-mom-rel'] , label = 'x-mom-rel')
    semilogy(data['step'], data['energy-rel'], label = 'energy-rel')
    semilogy(data['step'], data['nuhat-rel'] , label = 'nuhat-rel')
    xlabel('step')
    legend()
    show()
