# Author: Rowan J. Gollan & Peter J.
# Date: 2024-03-03
#
# Integration test for the 3D simulation of the sod shock tube (unstructured grid),
# exercising some 3D volume building and initial flow setting.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import yaml

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_prep_gas():
    cmd = "lmr prep-gas -i ideal-air.lua -o ideal-air.gas"
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
    assert abs(steps-75) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.0006)/0.0006 < 0.01, \
      "Failed to arrive at expected time on final step."
    print("reason=", reason)

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_probe_post_shock_region():
    expected_result = {'rho':0.2647, 'p':30.2e3, 'T':398.0, 'vel.x':293.0}
    cmd = 'lmr probe-flow --names=rho,p,T,vel.x --location=0.78,0.025,0.025'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    probe_result = yaml.safe_load(proc.stdout)
    for key in expected_result.keys():
        val = float(probe_result['pointdata'][0][key])
        v = expected_result[key]
        assert abs(val - v)/(abs(v)+1.0) < 0.01, "Failed to see correct "+key
    return

def test_probe_expanded_driver_region():
    expected_result = {'rho':0.4271, 'p':30.2e3, 'T':247.0, 'vel.x':293.0}
    cmd = 'lmr probe-flow --names=rho,p,T,vel.x --location=0.6,0.025,0.025'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    probe_result = yaml.safe_load(proc.stdout)
    for key in expected_result.keys():
        val = float(probe_result['pointdata'][0][key])
        v = expected_result[key]
        assert abs(val - v)/(abs(v)+1.0) < 0.01, "Failed to see correct "+key
    return

def test_cleanup():
    cmd = "rm -rf ./lmrsim ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
