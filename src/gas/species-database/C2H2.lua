db.C2H2 = {}
db.C2H2.atomicConstituents = {C=2,H=2,}
db.C2H2.charge = 0
db.C2H2.M = {
   value = 26.037280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H2.gamma = {
   value = 1.2321e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H2.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H2.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          8.08681094E-01,
          2.33615629E-02,
         -3.55171815E-05,
          2.80152437E-08,
         -8.50072974E-12,
          2.64289807E+04,
          1.39397051E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          4.14756964E+00,
          5.96166664E-03,
         -2.37294852E-06,
          4.67412171E-10,
         -3.61235213E-14,
          2.59359992E+04,
         -1.23028121E+00,
      }
   }
}
