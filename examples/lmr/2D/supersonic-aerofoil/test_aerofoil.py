# Author: Nick Gibbons
# Date: 2024-08-04
#
# Integration test for the BD aerofoil at 3 degrees AoA

import pytest
import subprocess
import sys
#sys.path.append('.') # So that we can import the below line when running in pytest
#sys.path.append('2D/supersonic-aerofoil') # So that we can import the below line when running in pytest
sys.path.append(__file__.replace('test_aerofoil.py', '.')) # So that we can import the below line when running in pytest
from loads import *

# Solution generated on 21st of August 2024
# KAD, 2025-08-21
# Values for macos added.
# Compiler: LLVM D compiler v 1.41.0 with LLVM 20.1.6
# CPU: Apple M1 Pro
if (sys.platform == 'linux'):
    REF_LOADS={"Alpha":  3.0,
               "Inviscid x-Force": 15.479138121629932,
               "Inviscid y-Force": 262.59714106272236,
               "Viscous x-Force": 30.459896093918754,
               "Viscous y-Force": 0.2528720799466463,
               "Pitching Moment": -4.584229142740343,}
else:
    REF_LOADS={"Alpha":  3.0,
               "Inviscid x-Force": 15.479127952676537,
               "Inviscid y-Force": 262.5972502859998,
               "Viscous x-Force": 30.45982960589959,
               "Viscous y-Force": 0.2528709172772733,
               "Pitching Moment": -4.584223054459798,}
    
# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_gen_grid():
    cmd = "lmr custom-script --job=gengrid.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_partition_grid():
    cmds = ["ugrid_partition grid.su2 mapped-cells 4 2",
            "mkdir -p su2grid",
            "mv block_0_grid.su2 block_1_grid.su2 block_2_grid.su2 block_3_grid.su2 su2grid/",
    ]

    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_prep():
    cmds = ["mkdir -p lmrsim",
            "cp mapped-cells lmrsim/mapped-cells",
            "lmr prep-gas -i gm-air.inp -o lmrsim/gm-air.lua",
            "lmr prep-grid --job=af-grid.lua",
            "lmr prep-flow --job=af.lua",
    ]


    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_run_steady():
    cmd = "mpirun -np 4 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    reason = ""
    steps = 0
    t = 0.0
    lines = proc.stdout.split("\n")
    for line in lines:
        if line.find("STOP-REASON") != -1:
            reason = ' '.join(line.split()[1:]).strip()
        if line.find("FINAL-STEP") != -1:
            steps = int(line.split()[1])
    assert reason.startswith("relative-global-residual-target"), \
      "Failed to stop for the expected reason:" + reason
    assert abs(steps-102) < 4, "Failed to take correct number of steps:" + str(steps)

def test_solution():
    data = wall_data()
    aoa = read_aoa('angle_of_attack.lua')
    com = read_com('centre_of_mass.lua')
    xforce, yforce, taux, tauy, M = compute_totals(data, com)

    alphadiff = abs(aoa - REF_LOADS["Alpha"])
    xforceiff = abs(xforce - REF_LOADS["Inviscid x-Force"])
    yforcediff = abs(yforce - REF_LOADS["Inviscid y-Force"])
    tauxdiff = abs(taux - REF_LOADS["Viscous x-Force"])
    tauydiff = abs(tauy - REF_LOADS["Viscous y-Force"])
    Mdiff = abs(M - REF_LOADS["Pitching Moment"])

    #print("alphadiff", alphadiff )
    #print("xforceiff", xforceiff )
    #print("yforcediff",yforcediff)
    #print("tauxdiff",  tauxdiff  )
    #print("tauydiff",  tauydiff  )
    #print("Mdiff",     Mdiff )

    assert alphadiff < 1e-6, "Computed incorrect angle of attack"
    assert xforceiff < 1e-6, "Computed incorrect Inviscid X-Force"
    assert yforcediff < 1e-6, "Computed incorrect Inviscid Y-Force"
    assert tauxdiff < 1e-6, "Computed incorrect Viscous X-Force"
    assert tauydiff < 1e-6, "Computed incorrect Viscous Y-Force"
    assert Mdiff < 1e-6, "Computed incorrect Pitching Moment"

def test_cleanup():
    cmd = "rm -rf ./lmrsim mapped-cells su2grid grid.su2 "
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
