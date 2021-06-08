-- Author: Nick Gibbons
-- Date: 2021, June
-- Check equilibrium air gas model is behaving appropriately
-- Usage instructions:
--    $ prep-gas air-5sp-1T.inp air-5sp-1T.lua
--    $ gas-calc eq-test.lua

spFile = "air-5sp-1T.lua"
p = 0.1*101.35e3 -- Pa
T = 2500 -- K
X0 = {N2=0.76, O2=0.23, N=0.0, O=0.0, NO=0.0}
Yst = {N2=0.7321963, O2=0.23281198, N=0.0, O=0.01160037,  NO=0.02339135}

function main()
   print("Begin therm_perf_gas_equilibrium test...")
   print("    Test type (pt), T=", T, "p=", p)
   local gasmodel = GasModel:new{spFile}
   local eqcalc = EquilibriumCalculator:new{filename=spFile}

   local Q = gasmodel:createGasState()
   Q.p = p
   Q.T = T
   Q.massf = gasmodel:molef2massf(X0)

   print("    Calling set_massf_from_pT...")
   eqcalc:set_massf_from_pT(Q, gasmodel)
   print("    Done\n")

   print("    Species   Target      Value     Diff (%)")
   print("    ----------------------------------------")
   for name,Ytarget in pairs(Yst) do
       Y = Q.massf[name]
       diff = (Y-Ytarget)/(1e-9+Ytarget)*100
       print(string.format("    %2s:     %10.8f  %10.8f  %6.2f", name, Ytarget, Y, diff))
   end
   print("... Done\n")

   print("Checking ps solve ...")
   s_target = 10000.0
   eqcalc:set_massf_and_T_from_ps(Q, s_target, gasmodel)

   s = eqcalc:get_entropy(Q, gasmodel)
   print(string.format("Entropy: %6.2f should be %6.2f", s, s_target))
end

main()

