# Author: Rowan J. Gollan
# Date: 2023-07-05
#
# Integration test for the 15 deg wedge case.

import pytest
import subprocess
import re
import os

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

expected_reason_for_stop = "relative-global-residual-target"
expected_number_steps = 33
expected_final_cfl = 2.072e+04
tolerance_on_cfl_check = 0.01
expected_number_steps_on_restart = 33
expected_final_cfl_on_restart = 2.072e+04
expected_restart_step = 26

def expected_output(proc, expected_n_steps, expected_final_cfl, check_start_step=False):
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
        if line.find("RESTART-STEP") != -1:
            restart_step = float(line.split()[1])
    assert reason == expected_reason_for_stop, "Failed to stop for the expected reason."
    assert steps == expected_n_steps, "Failed to take correct number of steps."
    assert abs(cfl - expected_final_cfl)/expected_final_cfl < tolerance_on_cfl_check, \
        "Failed to arrive at expected CFL value on final step."
    if (check_start_step):
        assert  restart_step == expected_restart_step, "Failed to restart at the expected step."

def test_prep():
    make_targets = ['gas', 'grid', 'init']
    for tgt in make_targets:
        cmd = "make " + tgt
        proc = subprocess.run(cmd.split())
        assert proc.returncode == 0, "Failed during: " + cmd


def test_run():
    cmd = "lmr run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_output(proc, expected_number_steps, expected_final_cfl)

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')


def test_restart():
    cmd = "lmr run -s 1"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_output(proc, expected_number_steps_on_restart, expected_final_cfl_on_restart, True)

def test_cleanup():
    cmd = "make clean"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

