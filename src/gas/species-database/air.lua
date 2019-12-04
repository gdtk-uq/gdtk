db.air = {}
db.air.atomicConstituents = {}
db.air.charge = 0
db.air.M = {
   value = 28.96e-3,
   units = 'kg/mol',
}
db.air.gamma = {
   value = 1.4,
   note = "valid at low temperatures"
}
db.air.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.air.sutherlandVisc = {
   mu_ref = 1.716e-5, 
   T_ref = 273.0,
   S = 111.0,
   reference = "Table 1-2, White (2006)"
}
db.air.sutherlandThermCond = {
   T_ref = 273.0, 
   k_ref = 0.0241, 
   S = 194.0,
   reference = "Table 1-3, White (2006)"
}
db.air.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      1.009950160e+04,
     -1.968275610e+02,
      5.009155110e+00,
     -5.761013730e-03,
      1.066859930e-05,
     -7.940297970e-09,
      2.185231910e-12,
     -1.767967310e+02,
     -3.921504225e+00
   },
   segment1 = { 
      2.415214430e+05,
     -1.257874600e+03,
      5.144558670e+00,
     -2.138541790e-04,
      7.065227840e-08,
     -1.071483490e-11,
      6.577800150e-16,
      6.462263190e+03,
     -8.147411905e+00
   }
}

