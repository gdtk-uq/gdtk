# fpreactor.py
# A simple fixed-pressure reactor.
# PJ & RJG 2025-08-16 adapted from fixed-volume reactor.
#
# To prepare:
#   $ prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
#   $ prep-chem nitrogen-2sp.lua nitrogen-2sp-2r-Keq.lua chem.lua
#
# To run:
#   $ python3 fpreactor.py

from gdtk.gas import GasModel, GasState, ThermochemicalReactor
import math

gm = GasModel("nitrogen-2sp.lua")
reactor = ThermochemicalReactor(gm, "chem.lua")

gs = GasState(gm)
p_fixed = 1.0e5 # Pa
gs.p = p_fixed
gs.T = 4000.0 # degree K
gs.molef = {'N2':2/3, 'N':1/3}
gs.update_thermo_from_pT()

tFinal = 600.0e-6 # s
t = 0.0
dt = 1.0e-6
dtSuggest = 1.0e-11
print("# Start integration")
f = open("fpreactor.data", 'w')
f.write('# 1:t(s)  2:T(K)  3:p(Pa)  4:massf_N2  5:massf_N  6:conc_N2  7:conc_N 8:rho(kg/m**3)\n')
f.write("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e %20.12e\n" %
        (t, gs.T, gs.p, gs.massf[0], gs.massf[1], gs.conc[0], gs.conc[1], gs.rho))
while t <= tFinal:
    # Allow the reactions to occur for a fixed-volume.
    # The basic reactor update changes only the mass fractions of the species.
    dtSuggest = reactor.update_state(gs, dt, dtSuggest)
    t = t + dt
    # dt = dtSuggest # uncomment this to get quicker stepping
    gs.update_thermo_from_rhou()
    # Allow the gas mixture to relax isentropically to the original pressure.
    Rgas = gs.R
    gmma = gs.gamma
    c = gs.p * math.pow(1.0/gs.rho, gmma)
    gs.rho = 1.0/math.pow(c/p_fixed, 1.0/gmma)
    gs.T = p_fixed / (Rgas*gs.rho)
    gs.update_thermo_from_rhoT()
    assert abs((p_fixed - gs.p)/(p_fixed + 1.0)) < 1.0e-6, \
        f"Pressures mismatch p_fixed={p_fixed}, gs.p={gs.p}"
    f.write("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e %20.12e\n" %
            (t, gs.T, gs.p, gs.massf[0], gs.massf[1], gs.conc[0], gs.conc[1], gs.rho))
f.close()
print("# Done.")
