# test_fpreactor.py
# A simple fixed-pressure reactor test for PyGasState.
# PJ & RJG 2019-11-25
#          2023-06-11 convert to pytest
#          2025-08-16 adapted from fixed-volume variant
#
# To run on Linux Mint 21.1:
# $ PYTHONDONTWRITEBYTECODE=1 pytest-3

from gdtk.gas import GasModel, PyGasState as GasState, ThermochemicalReactor
from math import isclose, pow
import os, subprocess
import pytest

# This is used to change to local directory so that subprocess runs nicely.
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)


def test_0_fixed_pressure_reactor():
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
    p_fixed = 1.0e5 # Pa
    gs.p = p_fixed
    gs.T = 4000.0 # degree K
    gs.molef = {'N2':2/3, 'N':1/3}
    gs.update_thermo_from_pT()
    #
    tFinal = 600.0e-6 # s
    t = 0.0
    dt = 1.0e-6
    dtSuggest = 1.0e-11
    while t <= tFinal:
        dtSuggest = reactor.update_state(gs, dt, dtSuggest)
        t = t + dt
        # dt = dtSuggest # uncomment this to get quicker stepping
        gs.update_thermo_from_rhou()
        # Allow the gas mixture to relax isentropically to the original pressure.
        Rgas = gs.R
        gmma = gs.gamma
        c = gs.p * pow(1.0/gs.rho, gmma)
        gs.rho = 1.0/pow(c/p_fixed, 1.0/gmma)
        gs.T = p_fixed / (Rgas*gs.rho)
        gs.update_thermo_from_rhoT()
        assert abs((p_fixed - gs.p)/(p_fixed + 1.0)) < 1.0e-6, \
            f"Pressures mismatch p_fixed={p_fixed}, gs.p={gs.p}"
    #
    assert isclose(t, tFinal, rel_tol=1.0e-6)
    assert isclose(gs.T, 6001.86, rel_tol=1.0e-4)
    assert isclose(gs.p, 100.0e3, rel_tol=1.0e-4)
    assert isclose(gs.rho, 5.01465e-02, rel_tol=1.0e-4)
    assert isclose(gs.massf[0], 8.805569519277e-01, rel_tol=1.0e-4)
    assert isclose(gs.massf[1], 1.194430480723e-01, rel_tol=1.0e-4)
    assert isclose(gs.conc[0], 1.576277248560e+00, rel_tol=1.0e-4)
    assert isclose(gs.conc[1], 4.276278979184e-01, rel_tol=1.0e-4)
    #
    return
