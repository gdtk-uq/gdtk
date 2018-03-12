db.CH2CO = {}
db.CH2CO.atomicConstituents = {C=2,H=2,O=1,}
db.CH2CO.charge = 0
db.CH2CO.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CO.gamma = {
   value = 1.1908e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.13583630E+00,
          1.81188721E-02,
         -1.73947474E-05,
          9.34397568E-09,
         -2.01457615E-12,
         -7.04291804E+03,
          1.22156480E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          4.51129732E+00,
          9.00359745E-03,
         -4.16939635E-06,
          9.23345882E-10,
         -7.94838201E-14,
         -7.55105311E+03,
          6.32247205E-01,
      }
   }
}
