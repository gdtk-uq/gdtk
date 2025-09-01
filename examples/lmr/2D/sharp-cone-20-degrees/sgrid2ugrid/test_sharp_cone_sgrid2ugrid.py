# Author: Rowan J. Gollan & Peter J.
# Date: 2024-02-24
#
# Integration test for 20-degree sharp cone on structured grids.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
import shlex

pytestmark = pytest.mark.short

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_prep_gas():
    cmd = "lmr prep-gas -i ideal-air.lua -o ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_grid_steady():
    cmd = "lmr prep-grid --job=../sg-minimal/grid.lua"
    proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_sgrid2ugrid():
    cmd = "lmr structured2unstructured"
    proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_init_steady():
    cmd = "lmr prep-sim --job=steady.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run_steady():
    cmd = "lmr run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_cfl_check = 0.01
    expected_reason_for_stop = "relative-global-residual-target"
    # RJG, 2024-02-25
    # This case behaves subtly differently on linux and macos
    # The macos version takes one extra step to convergence,
    # but because of how we grow the CFL with a power law
    # that extra step makes quite a difference in expected CFL.
    #
    # So we specialise the expected values based on OS
    # KAD, 2025-08-20
    # Values for macos updated.
    # Compiler: LLVM D compiler v 1.41.0 with LLVM 20.1.6
    # CPU: Apple M1 Pro
    if (sys.platform == 'linux'):
        expected_number_steps = 46
        expected_final_cfl = 1.066e+04
    else:
        expected_number_steps = 45
        expected_final_cfl = 6.095e+03
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

def test_snapshot_steady():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_cleanup_steady():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
