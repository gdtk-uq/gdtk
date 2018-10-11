db.C2H5 = {}
db.C2H5.atomicConstituents = {C=2,H=5,}
db.C2H5.charge = 0
db.C2H5.M = {
   value = 29.061100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H5.gamma = {
   value = 1.1963e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H5.sigma = {
   value = 4.302,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H5.epsilon = {
   value = 252.300,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H5.Lewis = {
   value = 1.551
}
db.C2H5.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.30646568E+00,
         -4.18658892E-03,
          4.97142807E-05,
         -5.99126606E-08,
          2.30509004E-11,
          1.28416265E+04,
          4.70720924E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          1.95465642E+00,
          1.73972722E-02,
         -7.98206668E-06,
          1.75217689E-09,
         -1.49641576E-13,
          1.28575200E+04,
          1.34624343E+01,
      }
   }
}
