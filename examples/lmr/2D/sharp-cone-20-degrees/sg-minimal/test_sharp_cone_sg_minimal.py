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

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_prep_gas():
    cmd = "lmr prep-gas -i ideal-air.lua -o ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_grid_steady():
    cmd = "lmr prep-grid --job=grid.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_init_steady():
    cmd = "lmr prep-sim --job=steady.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run_steady():
    cmd = "lmr run --max-cpus=1"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_reason_for_stop = "relative-global-residual-target"
    expected_number_steps = 60
    if (sys.platform == 'linux'):
        expected_final_cfl = 6.226e+05
    else:
        expected_final_cfl = 5.787e+05
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
    assert abs(steps-expected_number_steps) < 9, "Failed to take correct number of steps."
    assert abs(cfl - expected_final_cfl)/expected_final_cfl < 0.01, \
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

def test_prep_grid_transient():
    cmd = "lmr prep-grid --job=grid.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_init_transient():
    cmd = "lmr prep-sim --job=transient.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run_transient():
    cmd = "lmr run"
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
        if line.find("FINAL-TIME") != -1:
            t = float(line.split()[1])
    assert reason.startswith("Reached target simulation time"), \
      "Failed to stop for the expected reason."
    assert abs(steps-833) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.005)/0.005 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_transient_restart():
    cmd = "lmr run -s 2"
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
        if line.find("FINAL-TIME") != -1:
            t = float(line.split()[1])
    assert reason.startswith("Reached target simulation time"), \
      "Failed to stop for the expected reason."
    assert abs(steps-333) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.005)/0.005 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot_transient():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_cleanup_transient():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

