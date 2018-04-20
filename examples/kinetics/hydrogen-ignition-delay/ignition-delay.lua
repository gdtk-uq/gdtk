spFile = "Stanford-2011-gas-model.lua"
reacFile = "Stanford-2011-reac-file.lua"
outFile = "Stanford-ignition-delay.dat"

tFinal = 1500.0e-6 -- s
pInit = P_atm
Tlow = 900.0 -- K
Thigh = 1300.0 -- K
dT = 10.0

igCriteria = 5.0e-3 -- mol/m^3 : OH

function ignition_delay(T, gm, chemUpdate)
   local Q = gm:createGasState()
   Q.p = pInit
   Q.T = T
   local total = 2 + 1 + 3.76
   local molef = {H2=2/total, O2=1/total, N2=3.76/total}
   Q.massf = gm:molef2massf(molef)
   gm:updateThermoFromPT(Q)

   local t = 0.0
   local dt = 1.0e-6
   local dtSuggest = 1.0e-11
   while t <= tFinal do
      dtSuggest = chemUpdate:updateState(Q, dt, dtSuggest, gm)
      t = t + dt
      dt = dtSuggest
      gm:updateThermoFromRHOE(Q)
      local conc = gm:massf2conc(Q)
      if conc.OH > igCriteria then
	 return t
      end
   end
   return false
end

function main()
   local gm = GasModel:new{spFile}
   local chemUpdate = ChemistryUpdate:new{filename=reacFile, gasmodel=gm}

   local f = assert(io.open(outFile, 'w'))
   f:write('# 1:T(K)  2:t(s)\n')

   for T=Tlow,Thigh,dT do
      local tIg = ignition_delay(T, gm, chemUpdate)
      if tIg then
	 f:write(string.format("%20.12e %20.12e\n", T, tIg))
      else
	 print("No ignition at T= ", T)
      end
   end

end

main()
