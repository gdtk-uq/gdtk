# unsteady-expansion-test.rb
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ ruby unsteady-expansion-test.rb
#
# PJ, 2019-12-01
#
$LOAD_PATH << '~/dgdinst/lib'
require 'gdtk/gas'

puts "Unsteady expansion."
gmodel = GasModel.new('cea-air13species-gas-model.lua')
state1 = GasState.new(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 320.0 # K  ideal air, not high T
state1.update_thermo_from_pT()
state1.update_sound_speed()
puts "  state1: %s" % [state1]
v1 = 0.0
jplus = v1 + 2*state1.a/(1.4-1)
puts "  v1=%g jplus=%g" % [v1,jplus]

puts "Finite wave process along a cplus characteristic, stepping in pressure."
state2 = GasState.new(gmodel)
flow = GasFlow.new(gmodel)
v2 = flow.finite_wave_dp(state1, v1, "cplus", 60.0e3, state2, 500)
puts "  v2=%g" % [v2]
puts "  state2: %s" % [state2]
puts "  ideal v2=%g" % [jplus - 2*state2.a/(1.4-1)]

puts "Finite wave process along a cplus characteristic, stepping in velocity."
v2 = flow.finite_wave_dv(state1, v1, "cplus", 125.0, state2)
puts "  v2=%g" % [v2]
puts "  state2: %s" % [state2]
puts "  ideal v2=%g" % [jplus - 2*state2.a/(1.4-1)]
