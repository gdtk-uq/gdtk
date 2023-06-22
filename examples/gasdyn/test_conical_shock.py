# test_conical_shock.py
#
# To run on Linux Mint 21.1:
# $ PYTHONDONTWRITEBYTECODE=1 pytest-3
#
# PJ, 2019-12-01
#     2023-06-12 convert to pytest
#
from gdtk.gas import GasModel, GasState, GasFlow
import math
import os
import subprocess
import shutil
import pytest

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_0_conical_shock():
    cmd = 'prep-gas ideal-air.inp ideal-air-gas-model.lua'
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0
    #
    m1 = 1.5
    print("Conical-shock test for m1=%g" % m1)
    #
    gmodel = GasModel('ideal-air-gas-model.lua')
    state1 = GasState(gmodel)
    state1.p = 100.0e3 # Pa
    state1.T = 300.0 # K ideal air, not high T
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    print("state1: %s" % state1)
    v1 = m1*state1.a
    print("v1=%g" % v1)
    #
    beta = math.radians(45.0)
    print("  given beta(degrees)=%g" % math.degrees(beta))
    state_c = GasState(gmodel)
    flow = GasFlow(gmodel)
    theta_c, v_c = flow.theta_cone(state1, v1, beta, state_c)
    print("  theta_c=%g degrees" % math.degrees(theta_c))
    print("  v_c=%g" % (v_c))
    print("  state_c: %s" % state_c)
    assert math.isclose(math.degrees(theta_c), 15.0, abs_tol=0.1), "cone angle"
    #
    print("Conical shock angle from deflection.")
    beta2 = flow.beta_cone(state1, v1, theta_c)
    print("  beta2(degrees)=%g" % math.degrees(beta2))
    assert math.isclose(beta, beta2, rel_tol=1.0e-2, abs_tol=1.0e-5), "shock wave angle"
    return

if __name__ == '__main__':
    test_0_conical_shock()
