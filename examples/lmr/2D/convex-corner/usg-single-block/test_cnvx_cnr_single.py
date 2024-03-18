# Author: Rowan J. Gollan
# Date: 2023-05-31
#
# Integration test for convex corner case on single block.

import pytest
import subprocess
import re
import os

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

expected_reason_for_stop = "relative-global-residual-target"
expected_number_steps = 28
expected_final_cfl = 6.611e+03
expected_restart_step = 21

def expected_output(proc, check_start_step=False):
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
    assert steps == expected_number_steps, "Failed to take correct number of steps."
    assert abs(cfl - expected_final_cfl)/expected_final_cfl < 0.005, \
        "Failed to arrive at expected CFL value on final step."
    if (check_start_step):
        assert  restart_step == expected_restart_step, "Failed to restart at the expected step."

def test_prep():
    make_targets = ['prep-gas', 'grid', 'init']
    for tgt in make_targets:
        cmd = "make " + tgt
        proc = subprocess.run(cmd.split())
        assert proc.returncode == 0, "Failed during: " + cmd


def test_run():
    cmd = "lmr run --max-cpus=1"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_output(proc)

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')


def test_restart():
    cmd = "lmr run -s 1"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_output(proc, True)

def test_cleanup():
    cmd = "make clean"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

