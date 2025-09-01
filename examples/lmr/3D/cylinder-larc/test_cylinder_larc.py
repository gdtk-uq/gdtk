# Author: Kyle A. Damm
# Date: 2024-09-13
#
# Integration test for the LARC cylinder simulation.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
import shutil
import numpy as np
from math import *
import importlib

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_partition_grid():
    cmds = [
	"ugrid_partition cylinder.su2 mapped_cells 8 3 true",
	"mkdir su2_grid",
	"mv block_0_cylinder.su2 block_1_cylinder.su2 block_2_cylinder.su2 block_3_cylinder.su2 block_4_cylinder.su2 block_5_cylinder.su2 block_6_cylinder.su2 block_7_cylinder.su2 su2_grid/",
	"mv mapped_cells su2_grid/"
    ]
    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_gas_and_kinetics():
    cmds = [
        "prep-gas gm-air-5sp-2T.inp air-5sp-2T.gas",
	"prep-chem air-5sp-2T.gas rm-park-5sp-2T.inp air-5sp-2T.chem",
	"prep-kinetics air-5sp-2T.gas air-5sp-2T.chem eem-park-5sp-2T.inp air-5sp-2T.exch"
    ]
    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_grid():
    cmd = "lmr prep-grid"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_init():
    cmd = "lmr prep-sim"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run():
    cmd = "mpirun --use-hwthread-cpus -np 8 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_cfl_check = 0.01
    expected_reason_for_stop = "relative-global-residual-target"
    # KAD, 2025-09-03
    # Values for macos added.
    # Compiler: LLVM D compiler v 1.41.0 with LLVM 20.1.6
    # CPU: Apple M1 Pro
    if (sys.platform == 'linux'):
        expected_number_steps = 316
        expected_final_cfl = 1.000e+06
    else:
        expected_number_steps = 310
        expected_final_cfl = 1.000e+06
    reason = ""
    steps = 0
    cfl = 0.0
    lines = proc.stdout.split("\n")
    for line in lines:
        if line.find("STOP-REASON") != -1:
            reason = line.split()[1]
        if line.find("FINAL-STEP") != -1:
            steps = int(line.split()[1])
        if line.find("FINAL-CFL") != -1:
            cfl = float(line.split()[1])
    assert reason == expected_reason_for_stop, "Failed to stop for the expected reason."
    assert abs(steps-expected_number_steps) < 5, "Failed to take correct number of steps."
    assert abs(cfl - expected_final_cfl)/expected_final_cfl < tolerance_on_cfl_check, \
        "Failed to arrive at expected CFL value on final step."

def test_cleanup():
    cmd = "rm -rf ./lmrsim su2_grid air-5sp-2T.gas air-5sp-2T.chem air-5sp-2T.exch __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd


if __name__=='__main__':
    test_check_jacobian()
