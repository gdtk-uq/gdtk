"""
Compute wall drag from a set of loads files

@author: Nick Gibbons
"""

from numpy import array, argsort, concatenate, mean, argmin, log, diff, log10
from sys import argv
from glob import glob
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})
plt.rcParams['svg.fonttype'] = 'none'

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


def wall_data(dirname='lmrsim', groupname='wall', time=-1):
    times = read_loads_file('{}/loads/loads-metadata'.format(dirname))
    finalindex = times['loads_index'][time]
    finalindex_string = str(finalindex).zfill(4)

    datafilepattern = '{dname}/loads/{ftime}/blk-*-{gname}.dat'.format(dname=dirname, ftime=finalindex_string, gname=groupname)
    datafilenames = glob(datafilepattern)
    datafiles = [read_loads_file(filename) for filename in datafilenames]

    data = collate_datafiles(datafiles)
    data['q_total'] *= -1.0*data['outsign']
    return data

def linear_fit(x, y):
    sx = x.sum()
    sy = y.sum()
    sx2= (x*x).sum()
    sxy= (x*y).sum()
    M = y.size
    m = (sxy - sx*sy/M)/(sx2 - sx*sx/M)
    c = (sy - sx*m)/M
    return m, c

def cut_out_segment(data, startx, endx):
    x = data['pos.x']
    s = argmin(x<startx)
    e = argmin(x<endx)
    xx = x[s:e].copy()
    cutout = {k:v[s:e].copy() for k,v in data.items()}

    return cutout

def get_orders_of_convergence(names, spacings, whts):
    """ We sort the incoming arrays and use the tightest spacing as the truth value """

    sortidxs = argsort(spacings)
    spacings = spacings[sortidxs].copy() 
    whts = whts[sortidxs].copy() 

    names = [names[i] for i in sortidxs][1:]
    errors = (whts - whts[0])[1:]
    dxs = spacings[1:]
    n = len(dxs)

    labels = []
    ps = []
    for i in range(0,n-1):
        label = "{}->{}".format(names[i+1], names[i])
        p = log(errors[i]/errors[i+1])/log(dxs[i]/dxs[i+1])
        ps.append(p)
        labels.append(label)

    return labels, ps

def get_orders_of_convergence2(names, spacings, whts):
    """ We sort the incoming arrays and use the tightest spacing as the truth value """


    n = len(whts)

    labels = []
    ps = []
    for i in range(0,n-2):
        truth = whts[i]
        error1 = whts[i+1] - truth
        error2 = whts[i+2] - truth
        label = "{}->{}".format(names[i+2], names[i+1])
        p = log(error1/error2)/log(spacings[i+1]/spacings[i+2])
        ps.append(p)
        labels.append(label)

    return labels, ps

if __name__=='__main__':
    directory_names = argv[1:]

    datas = []

    for directory_name in directory_names:
        print("Reading name: ", directory_name)
        data = wall_data(directory_name)
        datas.append(data)

    fig = plt.figure(figsize=(11,4.5))
    ax,ax2 = fig.subplots(1,2) 
    colours = ['blue','red','darkgreen','magenta', 'goldenrod', 'teal', 'olive']

    whts = []
    spacings = []
    cutouts = []
    for data in datas:
        cutout = cut_out_segment(data, 0.4-0.02, 0.4+0.02)
        whts.append(mean(cutout['q_total']))
        spacings.append(mean(cutout['cellWidthNormalToSurface']))
        cutouts.append(cutout)
    spacings = array(spacings)
    whts = array(whts)

    sortidxs = argsort(spacings)
    spacings = spacings[sortidxs].copy() 
    whts = whts[sortidxs].copy() 
    names = [directory_names[i] for i in sortidxs]

    labels, ps = get_orders_of_convergence2(names, spacings, whts)
    for l,p in zip(labels, ps):
        print("Convergence Order {}={}".format(l,p))

    for i in range(1, len(names)):
        change = (whts[0] - whts[i])/whts[0]*100
        pair = "{}->{}".format(names[i], names[0])
        print("Change in WHT {}: {:3.3f}%".format(pair, change))

    for j, dirs in enumerate(directory_names):
        data = datas[j]
        colour = colours[j]
        cutout = cutouts[j]
        
        ax.plot(data['pos.x']*1000,  data['q_total']/1e4,    color=colour, linewidth=1.5, linestyle='-', label=dirs)
        ax2.plot(data['pos.x']*1000,  data['q_total']/1e4,    color=colour, linewidth=2.0, linestyle='-', label=dirs)
        #ax.plot(cutout['pos.x']*1000,-1.0*cutout['q_total']/1e4,  color=colour, linewidth=1.0, linestyle='-', label=dirs)

    zoom_xmin = 400-40; zoom_xmax = 400+40
    zoom_ymin = 11; zoom_ymax = 16
    ax.plot([zoom_xmin, zoom_xmax, zoom_xmax, zoom_xmin, zoom_xmin],
            [zoom_ymin, zoom_ymin, zoom_ymax, zoom_ymax, zoom_ymin],
            color='black', linewidth=2.0, linestyle='--')

    ax.legend(framealpha=1.0)
    ax.set_ylabel('Heat Transfer (W/cm2)')
    ax.set_xlabel('X Position (mm)')
    ax.set_ylim((-4, 144))
    ax2.set_xlabel('X Position (mm)')
    ax2.set_xlim((zoom_xmin, zoom_xmax))
    ax2.set_ylim((zoom_ymin,zoom_ymax))
    ax.grid()
    plt.tight_layout()
    #plt.savefig('wht.svg')
    #plt.close()
    plt.show()

