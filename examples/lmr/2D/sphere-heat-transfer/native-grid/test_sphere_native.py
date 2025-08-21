# Author: Kyle A. Damm
# Date: 2024-08-23
#
# Integration test for heat transfer on a sphere discretized with an unstructured grid.
#
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
import shutil
from math import *
import numpy as np

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_prep_gas():
    cmd = "prep-gas gm-ideal-air.inp ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_grid():
    cmd = "lmr prep-grid"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_convert_grid():
    cmd = "lmr structured2unstructured"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_init():
    cmd = "lmr prep-sim"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run():
    cmd = "mpirun --use-hwthread-cpus -np 1 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_cfl_check = 0.01
    expected_reason_for_stop = "relative-global-residual-target"
    # KAD, 2025-08-20
    # Values for macos added.
    # Compiler: LLVM D compiler v 1.41.0 with LLVM 20.1.6
    # CPU: Apple M1 Pro
    if (sys.platform == 'linux'):
        expected_number_steps = 260
        expected_final_cfl = 1.660e+05
    else:
        expected_number_steps = 264
        expected_final_cfl = 2.319e+05
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

def read_file(filename):
    container = []
    with open(filename, 'r') as f:
        data = f.readlines()
        idx = data[0].split().index("q_total")
        print(idx)
        print(data[0])
        for line in data[1:]:
            dat = line.split()
            container.append(float(dat[idx]))
    return container

def test_solution():

    # file names
    if (sys.platform == 'linux'):
        ref_solution_file = "reference_solution_linux.dat"
        solution_file = "lmrsim/loads/0259/blk-0000-bndry-1-wall.dat"
    else:
        ref_solution_file = "reference_solution_macos.dat"
        solution_file = "lmrsim/loads/0263/blk-0000-bndry-1-wall.dat"

    # read in data
    ref_data = read_file(ref_solution_file)
    sol_data = read_file(solution_file)

    # check data has consistent dimension
    assert len(ref_data) == len(sol_data), "Reference solution and computed solution have inconsistent dimension."
    
    # compute error
    error = 0.0
    for i in range(0, len(sol_data)):
        error += (ref_data[i]-sol_data[i])**2
    error = sqrt(error)

    assert error < 1e-6, "Computed incorrect heat flux along sphere surface."

def test_cleanup():
    cmd = "rm -rf ./lmrsim ideal-air.gas __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

