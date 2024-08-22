# Author: Kyle A. Damm
# Date: 2024-08-21
#
# Integration test for supersonic ramp on unstructured grids.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
import shutil
from math import *
import numpy as np
import importlib
# the following lines enable us to import the local python module from any higher-level directory
# without having module name clash issues with any of the other automated tests
path_to_module = os.path.relpath(__file__).replace('test_supersonic_ramp.py', 'plot_over_line')
plot_over_line = importlib.import_module(path_to_module.replace('/', '.'))

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_partition_grid():
    cmd = "ugrid_partition ramp.su2 mapped_cells 1 2 true"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_gas():
    cmd = "prep-gas gm-ideal-air.inp ideal-air.gas"
    proc = subprocess.run(cmd.split())
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
    cmd = "mpirun --use-hwthread-cpus -np 1 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_cfl_check = 0.01
    expected_reason_for_stop = "relative-global-residual-target"
    expected_number_steps = 97
    expected_final_cfl = 6.057e+04
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

def test_snapshot():
    cmd = "lmr snapshot2vtk --final"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_gradient():
    cmd = "lmr gradient2vtk --final"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_limiter():
    cmd = "lmr limiter2vtk --final"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_solution():

    # file names
    ref_solution_file = "reference_solution.dat" # generated on 2024-08-21 using commit a97145f696a7f10db47a0ff3cc2d508f16fec7e7
    solution_file     = "solution_over_line.dat"

    # call post process from plot_over_line.py
    plot_over_line.post_process(solution_file)

    # read in data
    ref_data = plot_over_line.read_file(ref_solution_file)
    sol_data = plot_over_line.read_file(solution_file)

    # check data has consistent dimension
    assert len(ref_data['rho']) == len(sol_data['rho']), "Reference solution and computed solution have inconsistent dimension."
    
    # compute error
    error = 0.0
    for i in range(0, len(sol_data['rho'])):
        error += (ref_data['rho'][i]-sol_data['rho'][i])**2
    error = sqrt(error)

    assert error < 1e-6, "Computed incorrect density profile over line."

def test_cleanup():
    cmd = "rm -rf ./lmrsim mapped_cells block_0_ramp.su2 solution_over_line.dat ideal-air.gas __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd


