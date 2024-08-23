# Author: Peter J. anf Rowan G.
# Date: 2024-08-23
#
# Integration test for the moving-grid example of a piston with fixed velocity.
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

def test_prep_grid():
    cmd = "lmr prep-grid --job=grid.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_sim():
    cmd = "lmr prep-sim --job=transient.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run():
    cmd = "lmr run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    reason = ""
    steps = 0
    t = 0.0
    lines = proc.stdout.split("\n")
    print(lines)
    for line in lines:
        if line.find("STOP-REASON") != -1:
            reason = ' '.join(line.split()[1:]).strip()
        if line.find("FINAL-STEP") != -1:
            steps = int(line.split()[1])
        if line.find("FINAL-TIME") != -1:
            t = float(line.split()[1])
    assert reason.startswith("maximum-time"), \
      "Failed to stop for the expected reason."
    assert abs(steps-841) < 15, "Failed to take correct number of steps."
    assert abs(t - 0.040)/0.040 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_piston_motion():
    f = open('piston.data', 'r')
    lines = f.readlines()
    f.close()
    items = lines[-1].strip().split()
    t = float(items[0])
    x = float(items[1])
    v = float(items[2])
    assert abs(t - 0.040) < 0.0001, "Failed to reach final time."
    assert abs(x - 6.550) < 0.1, "Failed to reach correct position."
    assert abs(v - 276.6) < 1.0, "Failed to reach correct velocity."

def test_energy_balance():
    cmd = 'lmr custom-script --job=balanceCheck.lua'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    text = proc.stdout
    lines = proc.stdout.split("\n")
    while not lines[0].strip().startswith("Energy-error ="): lines = lines[1:]
    energy_error = float(lines[0].split()[2])
    assert abs(energy_error) < 500.0, "Failed to get small energy error."

def test_cleanup_transient():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas piston.data userPadSave.lua"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

