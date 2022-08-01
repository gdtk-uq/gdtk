# A script to output Cp and h for O2 over temperature range 200--20000 K.
#
# Author: Peter J. and Rowan J. Gollan
# Date: 2019-11-21
#
# To run this script:
# $ prep-gas O2.inp O2-gas-model.lua
# $ ruby thermo-curves-for-O2.rb
#
$LOAD_PATH << '~/dgdinst/lib'
require 'gdtk/gas'

gasModelFile = 'O2-gas-model.lua'
gmodel = GasModel.new(gasModelFile)

gs = GasState.new(gmodel)
gs.p = 1.0e5 # Pa
gs.massf = {"O2"=>1.0}

outputFile = 'O2-thermo.dat'
puts "Opening file for writing: %s" % outputFile
f = open(outputFile, "w")
f.write("#  1:T[K]      2:Cp[J/kg/K]     3:h[J/kg]\n")

lowT = 200.0
dT = 100.0

(0..198).each do |i|
  gs.T = dT*i + lowT
  gs.update_thermo_from_pT()
  f.write(" %12.6e %12.6e %12.6e\n" % [gs.T, gs.Cp, gs.enthalpy])
end

f.close()
puts "File closed. Done."
