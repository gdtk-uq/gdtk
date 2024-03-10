# Author: Rowan J. Gollan & Peter J.
# Date: 2024-03-08
#
# Integration test for shock-fitting and 2T chemistry on structured grids.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import shutil

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

gas_model_files = ['air-5sp-gas-model.lua',
                   'GuptaEtAl-air-reactions-2T.lua',
                   'air-energy-exchange.lua']

def test_prep_gas():
    repo = os.getenv('DGD_REPO')
    locn = 'examples/kinetics/air-chemistry-2T'
    for f in gas_model_files:
        shutil.copyfile(os.path.join(repo, locn, f), f)
    f = open('GuptaEtAl-air-reactions-2T.lua', 'r')
    lines = f.readlines()
    f.close()
    f = open('GuptaEtAl-air-reactions-2T.lua', 'w')
    for line in lines:
        if line.startswith('SUBSET_SELECTION'):
            changedline = re.sub('11', '5', line)
            f.write(changedline)
        else:
            f.write(line)
    f.close()
    cmd = "lmr prep-gas -i air-5sp-gas-model.lua -o air-5sp-2T.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "lmr prep-reactions -g air-5sp-2T.gas -i GuptaEtAl-air-reactions-2T.lua -o air-5sp-6r-2T.chem"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "lmr prep-energy-exchange -g air-5sp-2T.gas -r air-5sp-6r-2T.chem -i air-energy-exchange.lua -o air-VT.exch"
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
    assert abs(steps-4745) < 5, "Failed to take correct number of steps."
    assert abs(t - 20.0e-6)/20.0e-6 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_shock_standoff():
    cmd = 'lmr custom-script --job=shock-shape.lua'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = 'python3 compute-error.py'
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    rms = 100.0  # something obviously large
    lines = proc.stdout.split("\n")
    for line in lines:
        if line.find("Shock Standoff RMS") != -1:
            t = float(line.strip().split()[-1])
    assert abs(t - 0.01205)/0.01205 < 0.01, \
      "Failed to compute expected shock shape."

def test_cleanup():
    for f in gas_model_files: os.remove(f)
    os.remove('air-5sp-2T.gas')
    os.remove('air-5sp-6r-2T.chem')
    os.remove('air-VT.exch')
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
