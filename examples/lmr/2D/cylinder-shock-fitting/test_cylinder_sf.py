# Author: Rowan J. Gollan & Peter J.
# Date: 2024-03-02
#
# Integration test for shock-fitting on structured grids.
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
    assert reason.startswith("Reached target simulation time"), \
      "Failed to stop for the expected reason."
    assert abs(steps-4067) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.016463)/0.016463 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_slice():
    cmd = 'lmr slice-flow --slice-list="0,:,0,0;2,:,0,0;4,:,0,0"'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    reason = ""
    steps = 0
    t = 0.0
    lines = proc.stdout.split("\n")
    # Discard comment lines and look at first data line for cell position.
    while lines[0].startswith("#"):
        lines = lines[1:]
    xpos = float(lines[0].split()[0])
    assert abs(xpos+1.41994) < 0.01, "Failed to move to correct shock position."
    # Discard tailing empty lines and pick stagnation pressure from last line.
    while len(lines[-1].strip()) == 0: lines = lines[0:-1]
    pstag = float(lines[-1].split()[3])
    # We make an estimate of the expected pressure from the free-stream conditions.
    rho_inf = 100.0e3/(287.1*300)
    V_inf = 2430.0
    pestimate = 0.926*rho_inf*V_inf*V_inf
    assert abs(pstag - pestimate)/pestimate < 0.01, "failed to see correct stagnation pressure."

def test_cleanup():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

