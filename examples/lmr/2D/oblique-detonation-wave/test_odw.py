# Author: Rowan J. Gollan & Peter J.
# Date: 2024-03-02
#
# Integration test for the oblique-detonation wave simulation,
# exercising finite-rate chemistry and custom-script processing.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


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
    assert abs(steps-2113) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.0200)/0.0200 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_custom_script():
    cmd = 'lmr custom-script --job=estimate_shock_angle.lua'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    lines = proc.stdout.split("\n")
    shock_angle_deg = 0.0
    average_deviation_metres = 1.0
    for line in lines:
        if line.find("shock_angle_deg") != -1:
            shock_angle_deg = float(line.split()[1])
        if line.find("average_deviation_metres") != -1:
            average_deviation_metres = float(line.split()[1])
    assert abs(shock_angle_deg - 45.139) < 0.1, "Failed to get correct shock angle."
    assert abs(average_deviation_metres) < 0.005, "Shock is not straight enough."

def test_cleanup():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
