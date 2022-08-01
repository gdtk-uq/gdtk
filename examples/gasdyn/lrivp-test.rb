# lrivp-test.rb
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ ruby lrivp-test.rb
#
# PJ, 2020-03-30 adpated from osher-riemann-test.rb
#
$LOAD_PATH << '~/dgdinst/lib'
require 'gdtk/gas'

# Set up something like Sod's shock tube.
gmodel = GasModel.new('ideal-air-gas-model.lua')
state1 = GasState.new(gmodel)
state1.p = 10.0e3 # Pa
state1.T = 278.8 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
puts "state1: %s" % state1
state4 = GasState.new(gmodel)
state4.p = 100.0e3 # Pa
state4.T = 348.4 # K
state4.update_thermo_from_pT()
state4.update_sound_speed()
puts "state4: %s" % state4

puts "Solve Lagrangian flavour of Riemann problem"
flow = GasFlow.new(gmodel)
pstar, wstar = flow.lrivp(state4, state1, 0.0, 0.0)
puts "pstar=%g wstar=%g" % [pstar, wstar]

puts "Solve piston-at-left problem (given constact-surface speed)"
pstar = flow.piston_at_left(state1, 0.0, 293.4)
puts "pstar=%g" % (pstar)

puts "Solve piston-at-right problem (given constact-surface speed)"
pstar = flow.piston_at_right(state4, 0.0, 293.4)
puts "pstar=%g" % (pstar)
