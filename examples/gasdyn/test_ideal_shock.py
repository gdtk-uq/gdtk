# test_ideal_shock.py
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


def test_0_ideal_shock():
    cmd = 'prep-gas ideal-air.inp ideal-air-gas-model.lua'
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0
    #
    gmodel = GasModel('ideal-air-gas-model.lua')
    state1 = GasState(gmodel)
    state1.p = 125.0e3 # Pa
    state1.T = 300.0 # K
    state1.update_thermo_from_pT()
    state1.update_sound_speed()
    print("state1: %s" % state1)
    print("normal shock (in ideal gas), given shock speed")
    vs = 2414.0
    state2 = GasState(gmodel)
    flow = GasFlow(gmodel)
    v2, vg = flow.ideal_shock(state1, vs, state2)
    print("v2=%g vg=%g" % (v2, vg))
    print("state2: %s" % state2)
    assert math.isclose(state2.p, 7.0268e6, rel_tol=1.0e-3), "post-shock pressure"
    assert math.isclose(state2.T, 3101.5, rel_tol=1.0e-3), "post-shock temperature"
    assert math.isclose(v2, 443.96, rel_tol=1.0e-3), "post-shock velocity"
    assert math.isclose(vg, 1970.04, rel_tol=1.0e-3), "post-shock velocity (lab frame)"
    return

if __name__ == '__main__':
    test_0_ideal_shock()
