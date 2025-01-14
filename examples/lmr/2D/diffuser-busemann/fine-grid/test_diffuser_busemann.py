# Author: Rowan J. Gollan
# Date: 2024-08-06
#
# Integration test for diffuser Busemann.
# Note: This fine grid depends on a completed coarse-grid solution.
#       This complicates the testing a little, but this is precisely
#       the workflow sequence we want to exercise with this test:
#       computing a coarse-grid solution and using it to warm-start
#       a simulation on a finer grid

import pytest
import subprocess
import re
import os
import shlex
import sys

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

expected_reason_for_stop = "relative-global-residual-target"
expected_final_cfl = 1.000e+03
# RJG, 2025-01-14
# The expected number of steps differs on linux and macos due to
# floating-point precision subtleties. This affects growth of
# CFL and the test for convergence based on relative residuals.
#
# Values for MacOs
# ----------------
# Compiler: LLVM D Compiler v1.40.0 with LLVM 19.1.6
# CPU: Apple M2 Pro
if (sys.platform == 'linux'):
    expected_number_steps = 39
else:
    expected_number_steps = 46

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
    assert abs(steps - expected_number_steps)<2, "Failed to take correct number of steps."
    assert reason == expected_reason_for_stop, "Failed to stop for the expected reason."
    assert abs(cfl - expected_final_cfl)/expected_final_cfl < 0.005, \
        "Failed to arrive at expected CFL value on final step."
    if (check_start_step):
        assert  restart_step == expected_restart_step, "Failed to restart at the expected step."

def test_prepare_coarse_grid_solution():
    cg_dir = "../coarse-grid"
    make_targets = ['prep', 'run']
    for tgt in make_targets:
        cmd = f"cd {cg_dir}; make {tgt}"
        proc = subprocess.run(cmd, shell=True)
        assert proc.returncode == 0, "Failed during: " + cmd


def test_prep():
    make_targets = ['gas', 'grid', 'prep-sim']
    for tgt in make_targets:
        cmd = "make " + tgt
        proc = subprocess.run(shlex.split(cmd))
        assert proc.returncode == 0, "Failed during: " + cmd

def test_run():
    cmd = "mpirun -np 4 lmrZ-mpi-run"
    proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_output(proc)


def test_cleanup():
    # first clean coarse-grid area
    cmd = "cd ../coarse-grid; make clean"
    proc = subprocess.run(cmd, shell=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    # then clean in this directory
    cmd = "make clean"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

