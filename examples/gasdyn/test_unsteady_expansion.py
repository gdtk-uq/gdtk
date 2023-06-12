# test_unsteady_expansion.py
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


def test_0_unsteady_expansion():
    repo = os.getenv('DGD_REPO')
    gas_model_file = 'cea-air13species-gas-model.lua'
    shutil.copyfile(os.path.join(repo, 'src/gas/sample-data', gas_model_file), gas_model_file)
    #
    gmodel = GasModel(gas_model_file)
    state1 = GasState(gmodel)
    state1.p = 100.0e3 # Pa
    state1.T = 320.0 # K  ideal air, not high T
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    print("  state1: %s" % state1)
    v1 = 0.0
    jplus = v1 + 2*state1.a/(1.4-1)
    print("  v1=%g jplus=%g" % (v1,jplus))
    #
    print("Finite wave process along a cplus characteristic, stepping in pressure.")
    state2 = GasState(gmodel)
    flow = GasFlow(gmodel)
    v2 = flow.finite_wave_dp(state1, v1, "cplus", 60.0e3, state2, 500)
    print("  v2=%g" % v2)
    print("  state2: %s" % state2)
    print("  ideal v2=%g" % (jplus - 2*state2.a/(1.4-1)))
    assert math.isclose(v2, 126.2, rel_tol=1.0e-3), "velocity after finite_wave_dp"
    assert math.isclose(state2.p, 60.0e3, rel_tol=1.0e-3), "pressure after finite_wave_dp"
    assert math.isclose(state2.T, 276.5, rel_tol=1.0e-3), "temperature after finite_wave_dp"
    #
    print("Finite wave process along a cplus characteristic, stepping in velocity.")
    v2 = flow.finite_wave_dv(state1, v1, "cplus", 125.0, state2)
    print("  v2=%g" % v2)
    print("  state2: %s" % state2)
    print("  ideal v2=%g" % (jplus - 2*state2.a/(1.4-1)))
    assert math.isclose(v2, 125.0, rel_tol=1.0e-3), "velocity after finite_wave_dv"
    assert math.isclose(state2.p, 60.23e3, rel_tol=0.5e-3), "pressure after finite_wave_dv"
    assert math.isclose(state2.T, 276.9, rel_tol=1.0e-3), "temperature after finite_wave_dv"
    return

if __name__ == '__main__':
    test_0_unsteady_expansion()
