db.CO2 = {}
db.CO2.atomicConstituents = {}
db.CO2.charge = 0
db.CO2.M = {
   value = 0.04401,
   units = 'kg/mol',
}
db.CO2.gamma = {
   value = 1.3,
   note = "valid at low temperatures"
}
db.CO2.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.CO2.sutherlandVisc = {
   mu_ref = 14.8e-6, 
   T_ref = 293.15,
   S = 240.0,
   reference = "Crane Company (1988) - Flow of fluids through valves, fittings and pipes"
}
db.CO2.sutherlandThermCond = {
   T_ref = 273.0, --these have not been updated
   k_ref = 0.0241, --these have not been updated
   S = 194.0,--these have not been updated
   reference = "Table 1-3, White (2006)"
}
db.CO2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper =  1000.0,
       coeffs = { 1.009950160e+04,--these have not been updated either
                 -1.968275610e+02,
                  5.009155110e+00,
	         -5.761013730e-03,
                  1.066859930e-05,
                 -7.940297970e-09,
                  2.185231910e-12,
                 -1.767967310e+02,
                 -3.921504225e+00 },
   },
   segment1 = { 
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {  2.415214430e+05,
                 -1.257874600e+03,
                  5.144558670e+00,
                 -2.138541790e-04,
                  7.065227840e-08,
                 -1.071483490e-11,
                  6.577800150e-16,
                  6.462263190e+03,
                 -8.147411905e+00 }
  }
}


