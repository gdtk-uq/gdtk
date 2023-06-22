# test_fvreactor.py
# A simple fixed-volume reactor test for PyGasState.
# PJ & RJG 2019-11-25
#          2023-06-11 convert to pytest
#
# To run on Linux Mint 21.1:
# $ PYTHONDONTWRITEBYTECODE=1 pytest-3

from gdtk.gas import GasModel, PyGasState as GasState, ThermochemicalReactor
from math import isclose
import os, subprocess
import pytest

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_0_finite_volume_reactor():
    cmd = 'prep-gas nitrogen-2sp.inp nitrogen-2sp.lua'
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0
    #
    cmd = 'prep-chem nitrogen-2sp.lua nitrogen-2sp-2r-Keq.lua chem.lua'
    proc = subprocess.run(cmd.split())
    assert proc.returncode == 0
    #
    gm = GasModel("nitrogen-2sp.lua")
    reactor = ThermochemicalReactor(gm, "chem.lua")
    #
    gs = GasState(gm)
    gs.p = 1.0e5 # Pa
    gs.T = 4000.0 # degree K
    gs.molef = {'N2':2/3, 'N':1/3}
    gs.update_thermo_from_pT()
    #
    tFinal = 300.0e-6 # s
    t = 0.0
    dt = 1.0e-6
    dtSuggest = 1.0e-11
    while t <= tFinal:
        dtSuggest = reactor.update_state(gs, dt, dtSuggest)
        t = t + dt
        # dt = dtSuggest # uncomment this to get quicker stepping
        gs.update_thermo_from_rhou()
    #
    assert isclose(t, tFinal, rel_tol=1.0e-6)
    assert isclose(gs.T, 6177.3, rel_tol=1.0e-4)
    assert isclose(gs.p, 145.518e3, rel_tol=1.0e-4)
    assert isclose(gs.massf[0], 8.692812671003e-01, rel_tol=1.0e-4)
    assert isclose(gs.massf[1], 1.307187328997e-01, rel_tol=1.0e-4)
    assert isclose(gs.conc[0], 2.178123112237e+00, rel_tol=1.0e-4)
    assert isclose(gs.conc[1], 6.550733441684e-01, rel_tol=1.0e-4)
    #
    return
