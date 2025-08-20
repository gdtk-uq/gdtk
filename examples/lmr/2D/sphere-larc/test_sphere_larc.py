# Author: Kyle A. Damm
# Date: 2024-08-21
#
# Integration test for the LARC sphere simulation.
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
	"ugrid_partition sphere.su2 mapped_cells 1 2",
	"mkdir su2_grid",
	"mv block_0_sphere.su2 su2_grid/",
	"mv mapped_cells su2_grid/"
    ]
    for cmd in cmds:
        proc = subprocess.run(cmd.split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_gas_and_kinetics():
    cmds = [
        "prep-gas gm-air-7sp-2T.inp air-7sp-2T.gas",
	"prep-chem air-7sp-2T.gas rm-park-7sp-2T.inp air-7sp-2T.chem",
	"prep-kinetics air-7sp-2T.gas air-7sp-2T.chem eem-park-7sp-2T.inp air-7sp-2T.exch"
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
    cmd = "mpirun --use-hwthread-cpus -np 1 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_cfl_check = 0.01
    expected_reason_for_stop = "relative-global-residual-target"
    expected_number_steps = 309
    expected_final_cfl = 5.641e+05
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

def test_check_jacobian():
    cmd = "lmrZ-check-jacobian --read-frozen-limiter-values=true --output=norms.dat"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_norm_check = 1.0e-14
    norm0 = 0.0
    norm1 = 0.0
    lines = proc.stdout.split("\n")
    for line in lines:
        if line.find("c0, 2-norm") != -1:
            norm0 = float(line.split()[2])
        if line.find("c1, 2-norm") != -1:
            norm1 = float(line.split()[2])
    assert abs((norm0 - norm1)/norm0) < tolerance_on_norm_check, \
        "Failed to compute consistent (Jacobian) matrix-vector multiplication effect."

def test_cleanup():
    cmd = "rm -rf ./lmrsim su2_grid air-7sp-2T.gas air-7sp-2T.chem air-7sp-2T.exch norms.dat __pycache__"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd


if __name__=='__main__':
    test_check_jacobian()
