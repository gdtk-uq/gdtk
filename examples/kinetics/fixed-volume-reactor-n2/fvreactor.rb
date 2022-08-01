# fvreactor.rb
# A simple fixed-volume reactor.
# PJ & RJG 2019-11-27
#
# To prepare:
#   $ prep-gas nitrogen-2sp.inp nitrogen-2sp.lua
#   $ prep-chem nitrogen-2sp.lua nitrogen-2sp-2r.lua chem.lua
#
# To run:
#   $ ruby fvreactor.py
$LOAD_PATH << '~/dgdinst/lib'
require 'gdtk/gas'

gm = GasModel.new("nitrogen-2sp.lua")
reactor = ThermochemicalReactor.new(gm, "chem.lua")

gs = GasState.new(gm)
gs.p = 1.0e5 # Pa
gs.T = 4000.0 # degree K
gs.molef = {'N2'=>2.0/3, 'N'=>1.0/3}
gs.update_thermo_from_pT()

tFinal = 300.0e-6 # s
t = 0.0
dt = 1.0e-6
dtSuggest = 1.0e-11
puts "# Start integration"
f = open("fvreactor.data", 'w')
f.write("# 1:t(s)  2:T(K)  3:p(Pa)  4:massf_N2  5:massf_N  6:conc_N2  7:conc_N\n")
f.write("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n" %
        [t, gs.T, gs.p, gs.massf[0], gs.massf[1], gs.conc[0], gs.conc[1]])
while t <= tFinal do
  dtSuggest = reactor.update_state(gs, dt, dtSuggest)
  t = t + dt
  # dt = dtSuggest # uncomment this to get quicker stepping
  gs.update_thermo_from_rhou()
  f.write("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n" %
          [t, gs.T, gs.p, gs.massf[0], gs.massf[1], gs.conc[0], gs.conc[1]])
end
f.close()
puts "# Done."
