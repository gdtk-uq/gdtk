db.C2H6 = {}
db.C2H6.atomicConstituents = {C=2,H=6,}
db.C2H6.charge = 0
db.C2H6.M = {
   value = 30.069040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H6.gamma = {
   value = 1.1872e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H6.sigma = {
   value = 4.302,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H6.epsilon = {
   value = 252.300,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H6.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.29142492E+00,
         -5.50154270E-03,
          5.99438288E-05,
         -7.08466285E-08,
          2.68685771E-11,
         -1.15222055E+04,
          2.66682316E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          1.07188150E+00,
          2.16852677E-02,
         -1.00256067E-05,
          2.21412001E-09,
         -1.90002890E-13,
         -1.14263932E+04,
          1.51156107E+01,
      }
   }
}
