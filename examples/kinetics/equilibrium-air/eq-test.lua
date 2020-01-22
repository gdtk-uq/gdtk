-- Authors: Nick Gibbons
-- Date: 2020, Jan
-- Check equilibrium air gas model is behaving appropriately

spFile = "air-5sp-eq-gas-model.lua"
p = 50e3 -- Pa
T = 1500 -- K
X0 = {N2=0.76, O2=0.23, N=0.0, O=0.0, NO=0.0}

function main()
   print("Making eq gas function")
   local gm = GasModel:new{spFile}
   print("Done")
end

main()

