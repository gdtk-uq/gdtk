# test_normal_shock.py
#
# To run on Linux Mint 21.1:
# $ PYTHONDONTWRITEBYTECODE=1 pytest-3
#
# PJ, 2019-11-28
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


def test_0_normal_shock():
    repo = os.getenv('DGD_REPO')
    gas_model_file = 'cea-air13species-gas-model.lua'
    # gas_model_file = 'air-5sp-eq.lua'
    shutil.copyfile(os.path.join(repo, 'src/gas/sample-data', gas_model_file), gas_model_file)
    #
    gmodel = GasModel(gas_model_file)
    state1 = GasState(gmodel)
    state1.p = 125.0e3 # Pa
    state1.T = 300.0 # K
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    print("state1: %s" % state1)
    #
    print("normal shock, given shock speed")
    vs = 2414.0
    print("vs=%g" % vs)
    state2 = GasState(gmodel)
    flow = GasFlow(gmodel)
    v2, vg = flow.normal_shock(state1, vs, state2)
    print("v2=%g vg=%g" % (v2, vg))
    print("state2: %s" % state2)
    assert math.isclose(v2, 361.9, rel_tol=1.0e-3), "v2 number after shock"
    assert math.isclose(vg, 2052.1, rel_tol=1.0e-3), "vg number after shock"
    assert math.isclose(state2.p, 7.314e6, rel_tol=1.0e-3), "p2 number after shock"
    assert math.isclose(state2.T, 2630.0, rel_tol=1.0e-3), "T2 number after shock"
    #
    print("normal shock, given pressure ratio")
    p2p1 = 58.516
    print("p2p1=%g" % p2p1)
    vs, v2, vg = flow.normal_shock_p2p1(state1, p2p1, state2)
    print("vs=%g v2=%g vg=%g" % (vs, v2, vg))
    print("state2: %s" % state2)
    assert math.isclose(vs, 2414.0, rel_tol=1.0e-3), "vs number after p2p1 shock"
    assert math.isclose(v2, 361.9, rel_tol=1.0e-3), "v2 number after p2p1 shock"
    assert math.isclose(vg, 2052.1, rel_tol=1.0e-3), "vg number after p2p1 shock"
    #
    print("reflected shock")
    state5 = GasState(gmodel)
    vr_b = flow.reflected_shock(state2, vg, state5)
    print("vr_b=%g" % vr_b)
    print("state5: %s" % state5)
    assert math.isclose(vr_b, 573.9, rel_tol=1.0e-3), "vr_b number after reflected shock"
    assert math.isclose(state5.p, 59.47e6, rel_tol=1.0e-3), "p5 number after reflected shock"
    assert math.isclose(state5.T, 4551.8, rel_tol=1.0e-3), "T5 number after reflected shock"
    return

if __name__ == '__main__':
    test_0_normal_shock()
