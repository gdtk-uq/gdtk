# Author: Nick Gibbons
# Date: 2025-09-15
#
# Integration test for laminar flat plate, for comparison to bloxide

import pytest
import subprocess
import sys
from re import sub
from gdtk.lmr import LmrConfig, SimInfo
from numpy import argmin, absolute, array, sqrt
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

POSX_TARGET = 280e-3 # how far along the plate to check the boundary layer

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber

def read_bloxide_file(filename):
    header = None
    body = []
    with open(filename) as fp:
        for line in fp:
            line = line.strip().split()
            if line[0].startswith('#'):
                header = line[1:]
            else:
                body.append(list(map(number,line)))
            
    assert header!=None
    header = [i.split(':')[-1] if ':' in i else i for i in header]
    cols = [array(col) for col in zip(*body)]
    data = dict(zip(header,cols))
    return data

def test_prep():

    cmds = [
            "lmr prep-gas -i ideal-air.lua -o ideal-air.gas",
            "lmr prep-grid --job=job.lua",
            "lmr prep-flow --job=job.lua",
    ]

    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_run():
    cmd = "mpirun -np 8 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    reason = ""
    steps = 0
    t = 0.0
    lines = proc.stdout.split("\n")
    for line in lines:
        if line.find("STOP-REASON") != -1:
            reason = ' '.join(line.split()[1:]).strip()
    assert reason.startswith("relative-global-residual-target"), \
      "Failed to stop for the expected reason:" + reason

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def read_bloxide_stdout(filename, searchstring):
    with open(filename) as fp:
        lines = [line for line in fp]

    for line in lines:
        if line.startswith(searchstring):
            target = float(line.split(':')[1].split()[0])
            return target
    else:
        raise Exception("Couldn't find line that starts with '{}' in file {}".format(searchstring, filename))


def test_loads(verbose=False):
    lmrcfg = LmrConfig()
    sim = SimInfo(lmrcfg)
    # Pick up the final loads files as a Pandas DataFrame.
    df = sim.read_loads(indx=sim.loads_indices[-1], group="wall")
    df = df.sort_values(by=['pos.x'])
    loads = {key:val.to_numpy() for key,val in df.items()}

    POSX_TARGET = 280e-3
    idxx = argmin(absolute(POSX_TARGET - loads['pos.x']))
    outsign = loads['outsign'][idxx]
    tauw = loads['tau_wall.x'][idxx]*outsign
    qw   = loads['q_total'][idxx]*outsign
    posx = loads['pos.x'][idxx]

    bloxideqw = read_bloxide_stdout('bloxide/bloxide.txt', 'Heat Transfer')
    # bloxide uses a different sign convention to Eilmer, and the units are reported
    # in W/cm2, so we do a little conversion here.
    bloxideqw *= -1e4
    bloxidetauw = read_bloxide_stdout('bloxide/bloxide.txt', 'Skin Friction')

    if verbose:
        print("qw", qw, "bloxide qw: ", bloxideqw, "({})".format(abs(qw - bloxideqw)/qw))
        print("tauw", tauw, "bloxide tau: ", bloxidetauw, "({})".format(abs(tauw - bloxidetauw)/tauw))

    assert abs(qw - bloxideqw)/qw < 0.1, "CFD heat transfer does not match self-similar solution"
    assert abs(tauw - bloxidetauw)/tauw < 0.01, "CFD skin friction does not match self-similar solution"

def test_profile(makeplot=False):
    lmrcfg = LmrConfig()
    sim = SimInfo(lmrcfg)

    flow = sim.read_snapshot(sim.snapshots[-1])
    data = flow.fields[3]
    ii = argmin(absolute(POSX_TARGET - data['pos.x'][:,0]))
    print(data['pos.x'].shape, ii, "x_slice=", data['pos.x'][ii,0])
    cfdbl = {key:val[ii,:].copy() for key,val in data.items()}

    x = cfdbl['pos.x']*1000.0
    y = cfdbl['pos.y']*1000.0
    velx = cfdbl['vel.x']
    T = cfdbl['T']

    bloxide = read_bloxide_file('bloxide/bloxide.dat')
    bloxide['y']*=1000.0

    inside_idxs = y<bloxide['y'].max()
    y_inside = y[inside_idxs]
    bl_interp_objs = {key:interp1d(bloxide['y'], val) for key,val in bloxide.items()}
    bl_interped = {key:interp(y_inside) for key,interp in bl_interp_objs.items()}
    bl_interped['y'] = y_inside

    n = bl_interped['vel'].size
    vrms = sqrt(((bl_interped['vel']-velx[inside_idxs])**2).sum()/n)
    Trms = sqrt(((bl_interped['T']-T[inside_idxs])**2).sum()/n)

    if makeplot:
        print("vrms", vrms, "Trms", Trms)

    assert vrms < 1.2, "CFD velocity profile does not match self-similar solution"
    assert Trms < 1.0, "CFD temperature profile does not match self-similar solution"

    if makeplot:
        fig = plt.figure(figsize=(9,4))
        ax = fig.subplots(1,2)

        fig.suptitle("Boundary Layer at x={:4.1f}mm".format(x[0]))
        ax[0].set_title("Velocity")
        ax[0].plot(velx, y, color='black', marker='o', markersize=4.0, linewidth=2.0, label="CFD")
        #ax[0].plot(bl_interped['vel'], bl_interped['y'], color='green', marker='o', markersize=4.0, linestyle="None")
        ax[0].plot(bloxide['vel'], bloxide['y'], 'b--', linewidth=2.0, label="bloxide")
        ax[0].set_ylabel('y (mm)')
        ax[0].set_xlabel('vel.x (m/s)')
        ax[0].grid()
        ax[0].legend(framealpha=1.0)

        ax[1].set_title("Temperature")
        ax[1].plot(T, y, color='black', marker='o', markersize=4.0, linewidth=2.0, label="CFD")
        #ax[1].plot(bl_interped['T'], bl_interped['y'], marker='o', markersize=4.0, linestyle="None", color="green")
        ax[1].plot(bloxide['T'], bloxide['y'], 'r--', linewidth=2.0, label="bloxide")
        ax[1].set_xlabel('T (K)')
        ax[1].grid()
        ax[1].legend(framealpha=1.0)
        ax[1].yaxis.set_ticklabels([])

        plt.tight_layout()
        plt.savefig('bl.svg')
        plt.show()

def test_cleanup():
    cmd = "rm -rf lmrsim ideal-air.gas __pycache__ bl.svg"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

if __name__=='__main__':
    test_prep()
    test_run()
    test_loads(verbose=True)
    test_profile(makeplot=True)

