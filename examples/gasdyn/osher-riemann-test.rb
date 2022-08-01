# osher-riemann-test.rb
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ ruby osher-riemann-test.rb
#
# PJ, 2020-02-07
#     2021-04-28 flux calculation
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

# Intermediate states
state2 = GasState.new(gmodel)
state3 = GasState.new(gmodel)
state0 = GasState.new(gmodel)

puts "Solve Riemann problem"
flow = GasFlow.new(gmodel)
pstar, wstar, wL, wR, velX0 = \
  flow.osher_riemann(state4, state1, 0.0, 0.0, state3, state2, state0)
puts "pstar=%g wstar=%g wL=%g wR=%g velX0=%g" % [pstar, wstar, wL, wR, velX0]
puts "state2: %s" % [state2]
puts "state3: %s" % [state3]
puts "state0: %s" % [state0]

puts "Compute flux at x0 directly."
f_mass, f_x_momentum, f_energy = \
    flow.osher_flux(state4, state1, 0.0, 0.0)
puts "F_mass=%g F_x_momentum=%g F_energy=%g" % [f_mass, f_x_momentum, f_energy]
