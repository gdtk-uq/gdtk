# classic-shock-tube.rb
#
# Run with commands like:
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ ruby classic-shock-tube.rb
#
# PJ, 2020-02-10 adapted from the Lua and Python scripts

$LOAD_PATH << '~/dgdinst/lib'
require 'eilmer/gas'
require 'eilmer/zero_solvers'
include ZeroSolvers

puts "Compute the flow conditions expected in the Sod shock tube."
#
puts "shock-tube fill conditions with air driving air"
gm = GasModel.new('ideal-air-gas-model.lua')
flow = GasFlow.new(gm)
states = []
5.times do states << GasState.new(gm) end
states[0].p = 1.0e4; states[0].T = 278.8 # spare, to be used later
states[1].p = 1.0e4; states[1].T = 278.8 # driven tube, initial
states[2].p = 1.0e4; states[2].T = 278.8 # intermediate, post-shock
states[3].p = 1.0e4; states[3].T = 278.8 # intermediate, post-expansion
states[4].p = 1.0e5; states[4].T = 348.4 # driver tube, initial
states.each do |gs| gs.update_thermo_from_pT() end

puts "state 1: %s" % states[1]
#
# For the unsteady expansion of the driver gas, regulation of the amount
# of expansion is determined by the shock-processed test gas.
# Across the contact surface between these gases, the pressure and velocity
# have to match so we set up some trials of various pressures and check 
# that velocities match.
error_in_velocity = Proc.new do |p3p4|
  # Compute the velocity mismatch for a given pressure ratio across the expansion.
  # Across the expansion, we get a test-gas velocity, v3g.
  p3 = p3p4*states[4].p
  v3g = flow.finite_wave_dp(states[4], 0.0, 'cplus', p3, states[3])
  # Across the contact surface.
  p2 = p3
  puts "current guess for p3 and p2=%g" % p2
  v1s, v2, v2g = flow.normal_shock_p2p1(states[1], p2/states[1].p, states[2])
  (v3g - v2g)/v3g
end
#
p3p4 = ZeroSolvers.secant(error_in_velocity, 0.105, 0.11, 1.0e-3)
puts "From secant solve: p3/p4=%g" % p3p4
puts "Expanded driver gas:"
p3 = p3p4*states[4].p
v3g = flow.finite_wave_dp(states[4], 0.0, 'cplus', p3, states[3])
puts "v3g=%g" % v3g
puts "state 3: %s" % states[3]
puts "Shock-processed test gas:"
v1s, v2, v2g = flow.normal_shock_p2p1(states[1], p3/states[1].p, states[2])
puts "v1s=%g v2g=%g" % [v1s, v2g]
puts("state 2: %s" % states[2])
if (v2g - v3g).abs/v3g > 1.0e-3 then
  raise "mismatch in velocities"
end
#
# Make a record for plotting against the Eilmer3 simulation data.
# We reconstruct the expected data along a tube 0.0 <= x <= 1.0
# at t=100us, where the diaphragm is at x=0.5.
x_centre = 0.5 # metres
t = 600.0e-6 # seconds
f = open('analytic.data', 'w')
f.write("# 1:x(m)  2:rho(kg/m**3) 3:p(Pa) 4:T(K) 5:V(m/s)\n")
puts 'Left end'
x = 0.0
f.write("%g %g %g %g %g\n" % [x, states[4].rho, states[4].p, states[4].T, 0.0])
puts 'Upstream head of the unsteady expansion.'
x = x_centre - states[4].a * t
f.write("%g %g %g %g %g\n" % [x, states[4].rho, states[4].p, states[4].T, 0.0])
puts 'The unsteady expansion in n steps.'
n = 100
dp = (states[3].p - states[4].p) / n
states[0].copy_values(states[4])
s = gm.entropy(states[4])
v = 0.0
p = states[4].p
n.times do
  rhoa = states[0].rho * states[0].a
  dv = -dp / rhoa
  v = v + dv
  p = p + dp
  states[0].p = p
  states[0].update_thermo_from_ps(s)
  x = x_centre + t * (v - states[0].a)
  f.write("%g %g %g %g %g\n" % [x, states[0].rho, states[0].p, states[0].T, v])
end
puts 'Downstream tail of expansion.'
x = x_centre + t * (v3g - states[3].a)
f.write("%g %g %g %g %g\n" % [x, states[3].rho, states[3].p, states[3].T, v3g])
puts 'Contact surface.'
x = x_centre + t * v3g
f.write("%g %g %g %g %g\n" % [x, states[3].rho, states[3].p, states[3].T, v3g])
x = x_centre + t * v2g  # should not have moved
f.write("%g %g %g %g %g\n" % [x, states[2].rho, states[2].p, states[2].T, v2g])
puts 'Shock front'
x = x_centre + t * v1s  # should not have moved
f.write("%g %g %g %g %g\n" % [x, states[2].rho, states[2].p, states[2].T, v2g])
f.write("%g %g %g %g %g\n" % [x, states[1].rho, states[1].p, states[1].T, 0.0])
puts 'Right end'
x = 1.0
f.write("%g %g %g %g %g\n" % [x, states[1].rho, states[1].p, states[1].T, 0.0])
f.close
