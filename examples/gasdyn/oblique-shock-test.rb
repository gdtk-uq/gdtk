# oblique-shock-test.rb
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ ruby oblique-shock-test.rb
#
# PJ, 2019-12-01
#
$LOAD_PATH << '~/dgdinst/lib'
require 'eilmer/gas'

m1 = 1.5
puts "Oblique-shock demo for m1=%g" % [m1]

gmodel = GasModel.new('cea-air13species-gas-model.lua')
state1 = GasState.new(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 300.0 # K ideal air, not high T
state1.update_thermo_from_pT()
state1.update_sound_speed()
puts "state1: %s" % state1
v1 = m1*state1.a
puts "v1=%g" % [v1]

beta = 45.0 * Math::PI/180.0
puts "  given beta(degrees)=%g" % [beta*180/Math::PI]
state2 = GasState.new(gmodel)
flow = GasFlow.new(gmodel)
theta, v2 = flow.theta_oblique(state1, v1, beta, state2)
puts "  theta=%g degrees" % [theta*180/Math::PI]
puts "  v2=%g" % [v2]
puts "  state2: %s" % [state2]

puts "Oblique shock angle from deflection."
beta2 = flow.beta_oblique(state1, v1, theta)
puts "  beta2(degrees)=%g" % [beta2*180/Math::PI]
