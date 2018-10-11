db.CH2 = {}
db.CH2.atomicConstituents = {C=1,H=2,}
db.CH2.charge = 0
db.CH2.M = {
   value = 14.0265800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH2.gamma = {
   value = 1.311,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH2.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2.Lewis = {
   value = 1.023
}
db.CH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.76267867E+00,
          9.68872143E-04,
          2.79489841E-06,
         -3.85091153E-09,
          1.68741719E-12,
          4.60040401E+04,
          1.56253185E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.87410113E+00,
          3.65639292E-03,
         -1.40894597E-06,
          2.60179549E-10,
         -1.87727567E-14,
          4.62636040E+04,
          6.17119324E+00,
      }
   }
}
