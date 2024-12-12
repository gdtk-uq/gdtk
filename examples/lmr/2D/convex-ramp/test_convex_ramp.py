# Author: Kyle A. Damm & Rowan J. Gollan & Peter J.
# Date: 2024-02-24
#
# Integration test for convex ramp on structured grids.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
import shutil

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

gas_model_file = 'air-5sp-gas-model.lua'

def test_prep_gas():
    repo = os.getenv('DGD_REPO')
    locn = 'examples/kinetics/air-chemistry-2T/'
    shutil.copyfile(os.path.join(repo, locn, gas_model_file), gas_model_file)
    cmd = "lmr prep-gas -i air-5sp-gas-model.lua -o air-5sp-2T.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

def test_prep_grid_steady():
    cmd = "lmr prep-grid --job=grid.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_init_steady():
    cmd = "lmr prep-sim --job=steady.lua"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd

def test_run_steady():
    cmd = "mpirun --use-hwthread-cpus -np 8 lmrZ-mpi-run"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    tolerance_on_cfl_check = 0.01
    expected_reason_for_stop = "relative-global-residual-target"
    # KAD & RJG, 2024-03-22
    #            2024-09-05
    # This case behaves subtly differently on linux and macos
    # in terms of how the CFL grows.
    #
    # It appears that the macos version needed to drop
    # the CFL on the way to convergence.
    #
    # So we specialise the expected values based on OS
    if (sys.platform == 'linux'):
        expected_number_steps = 246
        #expected_final_cfl = 2.808e+04
        expected_final_cfl = 13700.0 
    else:
        expected_number_steps = 246
        expected_final_cfl = 1.420e+04
    reason = ""
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
    assert reason == expected_reason_for_stop, "Failed to stop for the expected reason."
    assert abs(steps-expected_number_steps) < 5, "Failed to take correct number of steps."
    assert abs(cfl - expected_final_cfl)/expected_final_cfl < tolerance_on_cfl_check, \
        "Failed to arrive at expected CFL value on final step."

def test_snapshot_steady():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_cleanup_steady():
    os.remove(gas_model_file)
    os.remove('air-5sp-2T.gas')
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd


