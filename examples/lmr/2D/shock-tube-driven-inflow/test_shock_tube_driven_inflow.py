# Author: Peter J. and Rowan Gollan
# Date: 2024-03-03
#
# Integration test for the simulation of a shock tube with driven inflow,
# exercising InFlowBC_StaticProfile.
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

def test_prep_static_profile():
    cmd = "python3 prepare-static-profile.py"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_grid():
    cmd = "lmr prep-grid --job=grid.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_sim():
    cmd = "lmr prep-sim --job=transient-A.lua"
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
    assert abs(steps-435) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.0005)/0.0005 < 0.01, \
      "Failed to arrive at expected time on final step."
    print("reason=", reason)

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_probe_post_shock_region():
    # Expected values come from the analytic solution.
    expected_result = {'rho':0.0417124, 'p':7152.19, 'T':597.22, 'vel.x':587.33}
    cmd = 'lmr probe-flow --names=rho,p,T,vel.x --location=0.90,0.025,0.0'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    probe_result = yaml.safe_load(proc.stdout)
    for key in expected_result.keys():
        val = float(probe_result['pointdata'][0][key])
        v = expected_result[key]
        assert abs(val - v)/(abs(v)+1.0) < 0.01, "Failed to see correct "+key
    return

def test_probe_before_shock_region():
    expected_result = {'rho':0.0124931, 'p':1.0e3, 'T':278.8, 'vel.x':0.0}
    cmd = 'lmr probe-flow --names=rho,p,T,vel.x --location=0.95,0.025,0.0'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    probe_result = yaml.safe_load(proc.stdout)
    for key in expected_result.keys():
        val = float(probe_result['pointdata'][0][key])
        v = expected_result[key]
        assert abs(val - v)/(abs(v)+1.0) < 0.01, "Failed to see correct "+key
    return

def test_cleanup():
    cmd = "rm -rf ./lmrsim ideal-air.gas static-profile.data"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
