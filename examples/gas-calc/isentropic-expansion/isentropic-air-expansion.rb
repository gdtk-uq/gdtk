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
gs = GasState.new(gmodel)
gs.p = 500e3 # Pa
gs.T = 300.0 # K
gs.update_thermo_from_pT()
# Compute enthalpy and entropy at stagnation conditions
h0 = gs.enthalpy
s0 = gs.entropy
# Set up for stepping process
dp = 1.0 # Pa, use 1 Pa as pressure step size
gs.p = gs.p - dp
mach_tgt = 1.0
# Begin stepping until Mach = mach_tgt
while true do
  gs.update_thermo_from_ps(s0)
  h1 = gs.enthalpy
  v1 = Math.sqrt(2*(h0 - h1))
  gs.update_sound_speed()
  m1 = v1/gs.a
  if m1 >= mach_tgt then
    puts "Stopping at Mach=%g" % [m1]
    break
  end
  gs.p = gs.p - dp
end

puts "Gas properties at sonic point are:"
puts "p=%g T=%g" % [gs.p, gs.T]
