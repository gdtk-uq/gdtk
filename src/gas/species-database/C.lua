db.C = {}
db.C.atomicConstituents = {C=1,}
db.C.charge = 0
db.C.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C.gamma = { 
   value = 1.664,
   units = 'non-dimensional',
    description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.C.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.55423955E+00,
         -3.21537724E-04,
          7.33792245E-07,
         -7.32234889E-10,
          2.66521446E-13,
          8.54438832E+04,
          4.53130848E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.49266888E+00,
          4.79889284E-05,
         -7.24335020E-08,
          3.74291029E-11,
         -4.87277893E-15,
          8.54512953E+04,
          4.80150373E+00,
      }
   }
}
