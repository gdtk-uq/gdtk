# Author: Nick Gibbons
# Date: 2024-08-19
#
# Integration test for the "Ghostblood" turbulent flat plates, structured version

import pytest
import subprocess
import sys
from re import sub

sys.path.append(__file__.replace('test_splate.py', '.'))
from wall_heat_transfer import *

params = [
    {'niv':128, 'njv':48, 'a':8e-5},   # 0
    {'niv':128, 'njv':56, 'a':4e-5},   # 1
    {'niv':128, 'njv':64, 'a':2e-5},   # 2
    {'niv':128, 'njv':72, 'a':1e-5},   # 3
    {'niv':128, 'njv':86, 'a':5e-6},   # 4
    {'niv':128, 'njv':92, 'a':2.5e-6}, # 5
    {'niv':128, 'njv':108, 'a':1.25e-6}, # 6 This one went bad
]

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def prep(stage):
    niv = params[stage]['niv'] 
    njv = params[stage]['njv'] 
    a   = params[stage]['a'] 

    with open('job.lua') as fp:
        jobscript = fp.read()

    newscript = sub("niv = [0-9+]+;", "niv = {}+1;".format(niv), jobscript)
    newscript = sub("njv = [0-9+]+;", "njv = {}+1;".format(njv), newscript)
    newscript = sub("a = [0-9+-e]+;" ,"a = {:e};".format(a), newscript)

    with open('test.lua', 'w') as fp:
        fp.write(newscript)

    cmds = ["mkdir -p lmrsim",
            "lmr prep-gas -i ideal-air.lua -o ideal-air.gas",
            "lmr prep-grid --job=test.lua",
            "lmr prep-flow --job=test.lua",
    ]

    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def run_steady(stage):
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

    
    cmd = 'lmr extract-line -l "0.4,0.0,0.0,0.4,0.013,0.0,200" -f -o=lmrsim/line.txt'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

    a   = params[stage]['a'] 
    fch = "{}".format(int(a/1e-6)).zfill(2)
    cmd = "mv lmrsim {}um".format(fch)
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd


def test_stage_zero():
    prep(2)
    run_steady(2)

def test_stage_one():
    prep(3)
    run_steady(3)

def test_stage_two():
    prep(4)
    run_steady(4)


def test_wht_convergence_order():
    
    #directory_names = ['10um','05um','02um']
    directory_names = ['20um', '10um','05um']

    datas = []

    for directory_name in directory_names:
        data = wall_data(directory_name)
        datas.append(data)

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

    labels, ps = get_orders_of_convergence2(directory_names, spacings, whts)
    #for l,p in zip(labels, ps): print("Convergence {}={}".format(l,p))

    #target_p = 1.7967148081854336 # from the 10um, 05um, 02um triplet
    target_p = 1.915883199296974 # from the 20um, 10um, 05um triplet
    assert target_p - ps[0] < 1e-2, "Computed incorrect order of convergence for heat transfer"
    return

def test_cleanup():
    cmd = "rm -rf ./lmrsim 20um 10um 05um test.lua ideal-air.gas __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

if __name__=='__main__':
    # So we can run this as a normal script
    #for i in range(2,5): #len(params)):
    for i in range(len(params)):
        print("Running stage ", i)
        prep(i)
        run_steady(i)
