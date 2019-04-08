spFile = "C2H2-C2H4-model-gas.lua"
reacFile = "C2H2-C2H4-model-reactions-eilmer.lua"
outFile = "ignition-delay.dat"
outFile2 = "temp-history.dat"
tFinal = 3500e-6 -- s
--pInit = 60.795e03 --0.6atm SHock mixture 25
--pInit = 101325*3 --SHOCK MIXTURE 7
pInit = 3*101325 --SHOCK MIXTURE 10
Tlow = 1100.0 -- K
Thigh = 1700.0 -- K
dTemp = 50

function find_ignition_time(gm, chemUpdate, T)
   local Q = gm:createGasState()
   Q.p = pInit
   Q.T = T
   local init_temp = T
   tIg = 0.0
   --local molef = {C2H2=0.02, O2=0.025, Ar=0.955} --SHOCK MIXTURE 25
   --local molef = {C2H4=0.01, O2=0.015, Ar=0.975} --SHOCK MIXTURE 7
   local molef = {C2H4=0.0025, O2=0.015, Ar=0.9825} --SHOCK MIXTURE 25
   Q.massf = gm:molef2massf(molef)
   gm:updateThermoFromPT(Q)
   local t = 0.0
   local dt = 1.0e-7
   local dtSuggest = 1.0e-11
   local i = 1
   time = {}
   pressure = {}
   localTemp = {}
   if T > 1400 then
      tFinal = 750e-06
   end
   while t <= tFinal do
      dtSuggest = chemUpdate:updateState(Q, dt, dtSuggest, gm)
      dt = dtSuggest
      gm:updateThermoFromRHOU(Q)
      local conc = gm:massf2conc(Q)
      localTemp[i] = Q.T
      time[i] = t
      pressure[i] = Q.p
      t=t+dt
      i=i+1
   end    
   return time, pressure, localTemp
   
end

function main()
   local gm = GasModel:new{spFile}
   local chemUpdate = ChemistryUpdate:new{filename=reacFile, gasmodel=gm}
   local f = assert(io.open(outFile, 'w'))
   local ft = assert(io.open(outFile2, 'w'))
   print("Begin Ignition-delay-time Simulation")
   for T=Tlow,Thigh,dTemp do 
      pc = ((T-Tlow)/(Thigh-Tlow))*100
      print("Percentage Completed = ", pc, "%")
      local final = find_ignition_time(gm, chemUpdate, T)
      for z=1, 500000 do
         if not time[z] then
            break
         else
         --500000 is just a large aribitary number. I assume that the array is smaller than this.
         f:write(string.format("%20.12e , %20.12e ,\n", time[z], pressure[z]))
         ft:write(string.format("%20.12e , %20.12e ,\n", time[z], localTemp[z]))
         --print(string.format("%20.12e , %20.12e , %20.12e,\n", time[z], pressure[z], localTemp[z]))
         end
      end
      f:write(string.format("Temperature Completed at: %20.12e \n", T))
      ft:write(string.format("Temperature Completed at: %20.12e \n", T))
      --print(string.format("Temperature Completed at: %20.12e \n", T))

   end
  
end
main()
