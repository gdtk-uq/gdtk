"""
Pure python Method of Exact Solutions checker for electric field problem.

@author: Nick
"""

import pytest
import subprocess
from numpy import array, sin, exp, sqrt
from yaml import safe_load
import os
import gzip

def read_flow(blkid, ssid):
    with open("lmrsim/snapshots/fluid.metadata") as fp:
        md = {key:val for key,val in safe_load(fp).items()}
    
    blkstr = str(blkid).zfill(4)
    ssstr  = str(ssid).zfill(4)

    with gzip.open("lmrsim/snapshots/{}/fluid-{}.gz".format(ssstr,blkstr)) as fp: 
        flow = fp.read()
    flow = array([float(val) for val in flow.splitlines()])

    nvar = len(md['variables'])
    ncells = flow.size//nvar
    flow = flow.reshape((nvar,ncells))
    data = {key:var.copy() for key,var in zip(md['variables'], flow)}
    del flow
    return data

def target_field(x,y):
    return exp(x)*sin(y)

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_prep():
    cmds = [
        "lmr prep-gas -i ideal-air.inp -o ideal-air.lua",
        "lmr prep-grid --job=grid.lua",
        "lmr prep-flow --job=elec.lua",
    ]

    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_run():
    cmd = "lmr-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_exact_solution():
    data = read_flow(blkid=0, ssid=1)

    electric_potential = data['phi']
    x = data['pos.x']
    y = data['pos.y']
    target_potential = target_field(x,y)
    n = x.size

    difference = (target_potential - electric_potential)
    difference_squared = difference**2
    rms = sqrt(difference_squared.sum()/float(n))
    target_rms = 0.0002598072
    assert abs(rms - target_rms) < 1e-6, "Computed field does not match target exact solution"

def test_cleanup():
    cmd = "rm -rf lmrsim ideal-air.lua __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

if __name__=="__main__":
    data = read_flow(blkid=0, ssid=1)

    electric_potential = data['phi']
    x = data['pos.x']
    y = data['pos.y']
    target_potential = target_field(x,y)
    n = x.size

    difference = (target_potential - electric_potential)
    difference_squared = difference**2
    rms = sqrt(difference_squared.sum()/float(n))
    print("RMS Error: ", rms)

