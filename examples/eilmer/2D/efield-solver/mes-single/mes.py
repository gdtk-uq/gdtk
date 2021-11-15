"""
Pure python Method of Exact Solutions checker for electric field problem.
 - FIXME: This is a very crude way of reading the new flow format...

@author: Nick
"""

from numpy import array, cos, sin, exp, sqrt
from zipfile import ZipFile

def read_new_flow_file(filename, variablename):
    variablefile = '{}.dat'.format(variablename)
    with ZipFile(filename) as zf:
        with zf.open(variablefile) as fp:
            values = fp.read()

    data = array(list(map(float,values.strip().split())))
    return data

def test_field(x,y):
    return exp(x)*sin(y)

efieldfile = 'CellData/efield/t0001/elec.efield.b0000.t0001.zip'
flowfile = 'CellData/field/t0001/elec.field.b0000.t0001.zip'

electric_potential = read_new_flow_file(efieldfile, 'electric_potential')
x = read_new_flow_file(flowfile, 'pos.x')
y = read_new_flow_file(flowfile, 'pos.y')
target_potential = test_field(x,y)
n = x.size

difference = (target_potential - electric_potential)
difference_squared = difference**2
rms = sqrt(difference_squared.sum()/float(n))

print("RMS Error: ", rms)

