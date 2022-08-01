# conical-shock-test.rb
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ ruby conical-shock-test.rb
#
# PJ, 2019-12-01
#
$LOAD_PATH << '~/dgdinst/lib'
require 'gdtk/gas'

m1 = 1.5
puts "Conical-shock demo for m1=%g" % [m1]

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
state_c = GasState.new(gmodel)
flow = GasFlow.new(gmodel)
theta_c, v_c = flow.theta_cone(state1, v1, beta, state_c)
puts "  theta_c=%g degrees" % [theta_c*180/Math::PI]
puts "  v_c=%g" % [v_c]
puts "  state_c: %s" % [state_c]

puts "Conical shock angle from deflection."
beta2 = flow.beta_cone(state1, v1, theta_c)
puts "  beta2(degrees)=%g" % [beta2*180/Math::PI]
