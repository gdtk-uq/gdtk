# fvreactor.lua
# A simple fixed-volume reactor.
# PJ & RJG 2019-11-25
#
# To prepare:
#   $ prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
#   $ prep-chem nitrogen-2sp.lua nitrogen-2sp-2r.lua chem.lua
#
# To run:
#   $python3 fvreactor.py

from eilmer.gas import GasModel, GasState, ChemicalReactor

gm = GasModel("nitrogen-2sp.lua")
chem_update = ChemicalReactor("chem.lua", gm)

q = GasState(gm)
q.p = 1.0e5 # Pa
q.T = 4000.0 # degree K
q.molef = {'N2':2/3, 'N':1/3}
q.update_thermo_from_pT()

tFinal = 200.0e-6 # s
t = 0.0
dt = 1.0e-6
dtSuggest = 1.0e-11
print("# Start integration")
f = open("fvreactor.data", 'w')
f.write('# 1:t(s)  2:T(K)  3:p(Pa)  4:massf_N2  5:massf_N  6:conc_N2  7:conc_N\n')
f.write("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n" %
        (t, q.T, q.p, q.massf[0], q.massf[1], q.conc[0], q.conc[1]))
while t <= tFinal:
    dtSuggest = chem_update.update_state(q, dt, dtSuggest)
    t = t + dt
    # dt = dtSuggest # uncomment this to get quicker stepping
    q.update_thermo_from_rhou()
    f.write("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n" %
            (t, q.T, q.p, q.massf[0], q.massf[1], q.conc[0], q.conc[1]))
f.close()
print("# Done.")
