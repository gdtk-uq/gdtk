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
require 'gasmodule'

gasModelFile = 'O2-gas-model.lua'
gmodel = GasModel.new(gasModelFile)

Q = GasState.new(gmodel)
Q.p = 1.0e5 # Pa
Q.massf = {"O2"=>1.0}

outputFile = 'O2-thermo.dat'
puts "Opening file for writing: %s" % outputFile
f = open(outputFile, "w")
f.write("#  1:T[K]      2:Cp[J/kg/K]     3:h[J/kg]\n")

Tlow = 200.0
dT = 100.0

(0..198).each do |i|
  Q.T = dT*i + Tlow
  gmodel.update_thermo_from_pT(Q)
  f.write(" %12.6e %12.6e %12.6e\n" % [Q.T, Q.Cp, Q.enthalpy])
end

f.close()
puts "File closed. Done."
