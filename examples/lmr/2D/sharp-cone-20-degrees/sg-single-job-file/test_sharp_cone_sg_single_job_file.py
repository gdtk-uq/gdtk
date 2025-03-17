# Author: Rowan J. Gollan & Peter J.
# Date: 2024-02-24
#
# Integration test for 20-degree sharp cone on structured grids.
# Note that the tests are not independent and must be run in order of appearance.

import pytest
import subprocess
import re
import os
import sys
import shlex

pytestmark = pytest.mark.short

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
    assert abs(steps-833) < 5, "Failed to take correct number of steps."
    assert abs(t - 0.005)/0.005 < 0.01, \
      "Failed to arrive at expected time on final step."

def test_history():
    files = ["lmrsim/hist/hc-00-blk-0001-cell-9.dat.0",
             "lmrsim/hist/hc-01-blk-0001-cell-20.dat.0"]
    for fname in files:
        cmd = f"wc -l {fname}"
        proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
        lines = proc.stdout.split("\n")
        n_lines = int(lines[0].split()[0])
        expected_n_lines = 51
        assert n_lines == expected_n_lines, f"Failed to find expected number of lines in file: {fname}"

def test_snapshot():
    cmd = "lmr snapshot2vtk --all"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    assert os.path.exists('lmrsim/vtk')

def test_cleanup():
    cmd = "rm -r ./lmrsim"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd
    cmd = "rm ideal-air.gas"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

