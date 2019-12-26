# t4-reflected-shock-tunnel.rb
# Compute the state-to-state processes
# for a particular shot of the the T4 shock tunnel.
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ ruby t4-reflected-shock-tunnel.rb
#
# PJ, 2019-11-30
# 
$LOAD_PATH << '~/dgdinst/lib'
require 'eilmer/gas'

gmodel = GasModel.new('cea-air13species-gas-model.lua')
state1 = GasState.new(gmodel)
state1.p = 125.0e3 # Pa
state1.T = 300.0 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
puts "Initial test gas:"
puts "  state1: %s" % state1

puts "normal shock, given shock speed"
vs = 2414.0
puts "  vs=%g" % [vs]
state2 = GasState.new(gmodel)
flow = GasFlow.new(gmodel)
v2, vg = flow.normal_shock(state1, vs, state2)
puts "  v2=%g vg=%g" % [v2, vg]
puts "  state2: %s" % [state2]

puts "normal shock, given pressure ratio"
p2p1 = 58.516
puts "  p2p1=%g" % [p2p1]
vs, v2, vg = flow.normal_shock_p2p1(state1, p2p1, state2)
puts "  vs=%g v2=%g vg=%g" % [vs, v2, vg]
puts "  state2: %s" % [state2]

puts "reflected shock"
state5 = GasState.new(gmodel)
vr_b = flow.reflected_shock(state2, vg, state5)
puts "  vr_b=%g" % [vr_b]
puts "  state5: %s" % [state5]

puts "Expand from stagnation (with ratio of pressure to match observation)"
state5s = GasState.new(gmodel)
v5s = flow.expand_from_stagnation(state5, 34.37/59.47, state5s)
puts "  v5s=%g Mach=%g" % [v5s, v5s/state5s.a]
puts "  state5s: %s" % [state5s]
puts "  (h5s-h1)=%g" % [state5s.enthalpy - state1.enthalpy]

puts "Expand to throat condition (Mach 1.0001)"
state6 = GasState.new(gmodel)
v6 = flow.expand_to_mach(state5s, 1.0001, state6)
puts "  v6=%g Mach=%g" % [v6, v6/state6.a]
puts "  state6: %s" % [state6]
puts "  ceaSavedData=%s" % [state6.ceaSavedData]

puts "Something like Mach 4 nozzle"
state7 = GasState.new(gmodel)
v7 = flow.steady_flow_with_area_change(state6, v6, 27.0, state7)
puts "  v7=%g Mach=%g" % [v7, v7/state7.a]
puts "  state7: %s" % [state7]

puts "Total condition"
state8 = GasState.new(gmodel)
flow.total_condition(state7, v7, state8)
puts "  state8: %s" % [state8]

puts "Pitot condition"
state9 = GasState.new(gmodel)
flow.pitot_condition(state7, v7, state9)
puts "  state9: %s" % [state9]
