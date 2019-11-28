# normal-shock-test.rb
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ ruby normal-shock-test.rb
#
# PJ, 2019-11-28
# 
$LOAD_PATH << '~/dgdinst/lib'
require 'eilmer/gas'

gmodel = GasModel.new('cea-air13species-gas-model.lua')
state1 = GasState.new(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 300.0 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
puts "state1: %s" % state1

puts "normal shock, given shock speed"
vs = 2414.0
puts "vs=%g" % [vs]
state2 = GasState.new(gmodel)
flow = GasFlow.new(gmodel)
v2, vg = flow.normal_shock(state1, vs, state2)
puts "v2=%g vg=%g" % [v2, vg]
puts "state2: %s" % [state2]

puts "normal shock, given pressure ratio"
p2p1 = 58.516
puts "p2p1=%g" % [p2p1]
vs, v2, vg = flow.normal_shock_p2p1(state1, p2p1, state2)
puts "vs=%g v2=%g vg=%g" % [vs, v2, vg]
puts "state2: %s" % [state2]

puts "reflected shock"
state5 = GasState.new(gmodel)
vr_b = flow.reflected_shock(state2, vg, state5)
puts "vr_b=%g" % [vr_b]
puts "state5: %s" % [state5]
