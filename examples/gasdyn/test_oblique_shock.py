# test_oblique_shock.py
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


def test_0_oblique_shock():
    repo = os.getenv('DGD_REPO')
    gas_model_file = 'cea-air13species-gas-model.lua'
    shutil.copyfile(os.path.join(repo, 'src/gas/sample-data', gas_model_file), gas_model_file)
    #
    m1 = 1.5
    print("Oblique-shock demo for m1=%g" % m1)
    #
    gmodel = GasModel(gas_model_file)
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
    state2 = GasState(gmodel)
    flow = GasFlow(gmodel)
    theta, v2 = flow.theta_oblique(state1, v1, beta, state2)
    print("  theta=%g degrees" % math.degrees(theta))
    print("  v2=%g" % (v2))
    print("  state2: %s" % state2)
    assert math.isclose(math.degrees(theta), 2.78667, abs_tol=0.1), "wedge angle"
    #
    print("Oblique shock angle from deflection.")
    beta2 = flow.beta_oblique(state1, v1, theta)
    print("  beta2(degrees)=%g" % math.degrees(beta2))
    assert math.isclose(beta, beta2, rel_tol=1.0e-3), "shock-wave angle"
    return

if __name__ == '__main__':
    test_0_oblique_shock()
