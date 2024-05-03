# Author: Rowan J. Gollan & Peter J.
# Date: 2024-03-15
#
# Another integration test for shock-fitting on structured grids.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_prep_gas():
    cmd = "prep-gas ideal-air.lua ideal-air.gas"
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
    assert reason.startswith("maximum-time"), \
      "Failed to stop for the expected reason."
    assert abs(steps-509) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.0006)/0.0006 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_shock_location():
    cmd = 'lmr custom-script --job=shock-position.lua'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    x = 0.0
    lines = proc.stdout.split("\n")
    for line in lines:
        if line.find("location=") != -1:
            x = float(line.split()[1])
    assert abs(x - 0.1668)/0.1668 < 0.01, \
      "Failed to see expected shock location."

def test_cleanup():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

