# Author: Rowan J. Gollan
# Date: 2024-08-21
#
# Spot check on verification of k-omega model via Manufactured Solutions.
#
# This test compares on errors in a few key flow field quantities.
# The reference values for these errors were generated with revision:
#
#     76e87cbe9877a366ba2e9f1a0f64bae8dc120f89
#
# Reference output generated with:
#
# $ lmr compute-norms --norms='rho,vel.x,tq-tke,tq-omega' --reference-solution=ref-soln.lua --final
#

REF_OUTPUT_FROM_COMPUTE_NORMS = """
---
snapshot: 0002
field_type: fluid
error-norms-computed: yes
reference-solution-file: ref-soln.lua
rho:
   L1: 3.832147915973724110e-03
   L2: 4.075837778488668184e-03
   Linf: 9.155588503992784233e-03
   Linf-pos: {x: 9.375000e-01, y: 9.375000e-01, z: 9.375000e-01 }
vel.x:
   L1: 2.669499024611020621e-01
   L2: 3.167181137212501385e-01
   Linf: 1.136080229751740944e+00
   Linf-pos: {x: 6.250000e-02, y: 9.375000e-01, z: 6.875000e-01 }
tq-tke:
   L1: 1.100557346792062408e+00
   L2: 1.394327229755474074e+00
   Linf: 3.747313355211872477e+00
   Linf-pos: {x: 9.375000e-01, y: 5.625000e-01, z: 6.875000e-01 }
tq-omega:
   L1: 3.242030390844224397e-01
   L2: 4.551411514184353080e-01
   Linf: 1.939404442581690091e+00
   Linf-pos: {x: 9.375000e-01, y: 6.875000e-01, z: 9.375000e-01 }
...
"""

import pytest
import subprocess
import re
import os
import shlex
import yaml

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

expected_reason_for_stop = "relative-global-residual-target"
expected_number_steps = 66
expected_final_cfl = 8.813e+05
tolerance_on_cfl_check = 0.01
tolerance_on_norms = 1.0e-6

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

def expected_norms(proc):
    test_output = yaml.safe_load(proc.stdout)
    reference = yaml.safe_load(REF_OUTPUT_FROM_COMPUTE_NORMS)
    fields = ['rho', 'vel.x', 'tq-tke', 'tq-omega']
    norms = ['L1', 'L2', 'Linf']
    for var in fields:
        for norm in norms:
            t_val = float(test_output[var][norm])
            r_val = float(reference[var][norm])
            assert abs(t_val - r_val)/r_val < tolerance_on_norms, \
                f"Failure when comparing {norm} norm for flow field variable: {var}"

def test_prep():
    make_targets = ['prep-area', 'prep-grid', 'prep-sim']
    for tgt in make_targets:
        cmd = "make " + tgt
        proc = subprocess.run(shlex.split(cmd))
        assert proc.returncode == 0, "Failed during: " + cmd


def test_run():
    cmd = "lmr run"
    proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_output(proc, expected_number_steps, expected_final_cfl)

def test_norms():
    cmd = "lmr compute-norms --norms='rho,vel.x,tq-tke,tq-omega' --reference-solution=ref-soln.lua --final"
    proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed during: " + cmd
    expected_norms(proc)

def test_cleanup():
    cmd = "make clean"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

