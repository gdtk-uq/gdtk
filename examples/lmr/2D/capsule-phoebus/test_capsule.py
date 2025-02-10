# Author: Nick Gibbons
# Date: 2025-02-10
#
# Integration test for the steady-state shock fitting

import pytest
import subprocess
import sys
from re import sub

sys.path.append(__file__.replace('test_capsule.py', '.'))
from extract_shock_standoff import *


# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_prep():

    cmds = ["mkdir -p lmrsim",
            "lmr prep-gas -i ideal-air.inp -o ideal-air.gas",
            "lmr prep-grids --job=grid.lua",
            "lmr prep-flow --job=sim.lua",
    ]

    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_run_mpi():
    cmd = "mpirun --use-hwthread-cpus -np 8 lmrZ-mpi-run"
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
    assert abs(steps-978) < 5, "Failed to take correct number of steps:" + str(steps)

def test_shock_standoff():
    ss = get_shock_standoff()
    ss_in_mm = 1000.0*ss
    
    target_ss_in_mm = 2.978677601030077
    assert abs(target_ss_in_mm - ss_in_mm)  < 1e-2, "Computed incorrect shock standoff"
    return

def test_cleanup():
    cmd = "rm -rf ./lmrsim ideal-air.gas __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
