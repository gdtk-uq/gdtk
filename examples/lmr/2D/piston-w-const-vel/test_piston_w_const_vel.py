# Author: Peter J. anf Rowan G.
# Date: 2024-08-21
#
# Integration test for the moving-grid example of a piston with fixed velocity.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
# import pyyaml

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
    assert abs(steps-506) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.0006)/0.0006 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_probe():
    cmd = 'lmr probe-flow --location=0.28,0,0 --names=pos.x,T,vel.x,p'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    text = proc.stdout
    lines = proc.stdout.split("\n")
    # Pull apart the YAML output manually because
    # we don't want to depend upon the PyYAML package.
    # Look for first line with cell position.
    while not lines[0].strip().startswith("pos.x"): lines = lines[1:]
    posx = float(lines[0].split()[-1])
    assert abs(posx - 0.283) < 0.01, "Failed to report correct probe position."
    T = float(lines[1].split()[-1])
    assert abs(T - 398.5) < 1.0, "Failed to see correct temperature."
    velx = float(lines[2].split()[-1])
    assert abs(velx - 293.0) < 1.0, "Failed to see correct velocity."
    velx = float(lines[3].split()[-1])
    assert abs(velx - 30258.0) < 100.0, "Failed to see correct pressure."

def test_cleanup_transient():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

