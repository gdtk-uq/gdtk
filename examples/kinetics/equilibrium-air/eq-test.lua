-- Author: Nick Gibbons
-- Date: 2020, Jan
-- Check equilibrium air gas model is behaving appropriately
-- Usage instructions:
--    $ prep-gas air-5sp-eq.inp air-5sp-eq-gas-model.lua
--    $ gas-calc eq-test.lua

spFile = "air-5sp-eq-gas-model.lua"
p = 0.1*101.35e3 -- Pa
T = 2500 -- K
X0 = {N2=0.76, O2=0.23, N=0.0, O=0.0, NO=0.0}
Yst = {N2=0.7321963, O2=0.23281198, N=0.0, O=0.01160037,  NO=0.02339135}

function main()
   print("Begin therm_perf_gas_equilibrium test...")
   print("    Test type (pt), T=", T, "p=", p)
   local gm = GasModel:new{spFile}

    speciesNames = {}
    for i=1,gm:nSpecies() do
       speciesNames[i] = gm:speciesName(i-1)
    end

   local Q = gm:createGasState()
   Q.p = p
   Q.T = T
   Q.massf = gm:molef2massf(X0)
   print("    Calling updateThermoFromPT...")
   gm:updateThermoFromPT(Q)
   print("    Done\n")

   print("    Species   Target      Value")
   print("    -------------------------------")
   for i,name in ipairs(speciesNames) do
       print(string.format("    %2s:     %10.8f  %10.8f", name, Yst[name], Q.massf[name]))
   end

   print("... Done\n")
end

main()

