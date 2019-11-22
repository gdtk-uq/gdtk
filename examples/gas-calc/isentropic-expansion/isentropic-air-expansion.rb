# Step through a steady isentropic expansion,
# from stagnation condition to sonic condition.
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ ruby isentropic-air-expansion.rb
#
# Ruby port, PJ, 2019-11-21
#
$LOAD_PATH << '~/dgdinst/lib'
require 'eilmer/gas'

gmodel = GasModel.new('ideal-air-gas-model.lua')
q = GasState.new(gmodel)
q.p = 500e3 # Pa
q.T = 300.0 # K
q.update_thermo_from_pT()
# Compute enthalpy and entropy at stagnation conditions
h0 = q.enthalpy
s0 = q.entropy
# Set up for stepping process
dp = 1.0 # Pa, use 1 Pa as pressure step size
q.p = q.p - dp
mach_tgt = 1.0
# Begin stepping until Mach = mach_tgt
while true do
  q.update_thermo_from_ps(s0)
  h1 = q.enthalpy
  v1 = Math.sqrt(2*(h0 - h1))
  q.update_sound_speed()
  m1 = v1/q.a
  if m1 >= mach_tgt then
    puts "Stopping at Mach=%g" % [m1]
    break
  end
  q.p = q.p - dp
end

puts "Gas properties at sonic point are:"
puts "p=%g T=%g" % [q.p, q.T]
