# Author: Nick Gibbons and Lachlan Whyborn, with help from Daryl Bond
# Date: 2024-12-18
#
# Integration test for codename "Mercy" aka the Steepening Wave Problem from:
# "Fluid Mechanics", L. D. Landau and E. M. Lifshitz
# Chapter 10, Heinemann, 1987
#

import pytest
import subprocess
import sys
import shutil
import os
from re import sub

sys.path.append(__file__.replace('test_swp.py', '.'))
from post_process import *

target_orders = {
    "p_L0":3.92,
    "p_L2":3.96,
  "rho_L0":3.73,
  "rho_L2":3.57,
"vel.x_L0":3.70,
"vel.x_L2":3.81,
}
Ns = [32, 64, 128]

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def test_prep():
    with open('blank/swp.lua') as fp:
        simfile = fp.read()
    with open('blank/swp-grid.lua') as fp:
        gridfile = fp.read()

    if os.path.isdir('asf'):
        shutil.rmtree('asf')

    if not os.path.isdir('asf'):
        os.makedirs('asf')

    for N in Ns:
        dir = 'asf/{}'.format(str(N).zfill(3))
        shutil.copytree('blank', dir)

        newsimfile  = sub("N=[0-9+]+", "N={};".format(N), simfile)
        newgridfile = sub("N=[0-9+]+", "N={};".format(N), gridfile)
        with open(dir+'/swp.lua', 'w') as fp:
            fp.write(newsimfile)
        with open(dir+'/swp-grid.lua', 'w') as fp:
            fp.write(newgridfile)


        os.chdir(dir)
        cmds = [
            "lmr prep-gas -i ideal-air.lua -o ideal-air.gas",
            "lmr prep-grid --job=swp-grid.lua",
            "lmr prep-flow --job=swp.lua",
        ]

        for cmd in cmds:
            proc = subprocess.run(cmd.split(), capture_output=True, text=True)
            assert proc.returncode == 0, "Failed during: " + cmd
        os.chdir('../..')
    return

def test_run():
    for N in Ns:
        dir = 'asf/{}'.format(str(N).zfill(3))
        os.chdir(dir)
        proc = subprocess.run("lmr run".split(), capture_output=True, text=True)
        assert proc.returncode == 0, "Failed during: " + cmd
        os.chdir('../..')
    return


def test_order_of_accuracy():
    #dirs = ['asf/032', 'asf/064', 'asf/128']
    dirs = ['asf/{}'.format(str(N).zfill(3)) for N in Ns]
    variables, L0s, L2s = order_of_convergence(dirs)
    for var in variables:
        m = L0s[var+'m']
        #print(" {:>5s} L0={:6.6f}".format(var, m))
    
        m = L2s[var+'m']
        #print(" {:>5s} L2={:6.6f}".format(var, m))

    assert abs(target_orders['p_L0'] - L0s['pm']) < 1e-2, "Computed incorrect order of convergence for p_L0"
    assert abs(target_orders['p_L2'] - L2s['pm']) < 1e-2, "Computed incorrect order of convergence for p_L2"

    assert abs(target_orders['rho_L0'] - L0s['rhom']) < 1e-2, "Computed incorrect order of convergence for rho L0"
    assert abs(target_orders['rho_L2'] - L2s['rhom']) < 1e-2, "Computed incorrect order of convergence for rho L2"

    assert abs(target_orders['vel.x_L0'] - L0s['vel.xm']) < 1e-2, "Computed incorrect order of convergence for vel L0"
    assert abs(target_orders['vel.x_L2'] - L2s['vel.xm']) < 1e-2, "Computed incorrect order of convergence for vel L2"

def test_cleanup():
    cmd = "rm -rf asf"
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0, "Failed during: " + cmd

#test_prep()
#test_run()
#test_order_of_accuracy()

