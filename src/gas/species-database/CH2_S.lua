db.CH2_S = {}
db.CH2_S.atomicConstituents = {C=1,H=2,}
db.CH2_S.charge = 0
db.CH2_S.M = {
   value = 14.026580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2_S.gamma = {
   value = 1.3263e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2_S.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2_S.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2_S.Lewis = {
   value = 1.022
}
db.CH2_S.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.19860411E+00,
     -2.36661419E-03,
      8.23296220E-06,
     -6.68815981E-09,
      1.94314737E-12,
      5.04968163E+04,
     -7.69118967E-01,
   },
   segment1 = {
      0,
      0,
      2.29203842E+00,
      4.65588637E-03,
     -2.01191947E-06,
      4.17906000E-10,
     -3.39716365E-14,
      5.09259997E+04,
      8.62650169E+00,
   }
}
