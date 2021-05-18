db.CH2OH = {}
db.CH2OH.atomicConstituents = {C=1,H=3,O=1,}
db.CH2OH.charge = 0
db.CH2OH.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2OH.gamma = {
   value = 1.2070e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2OH.sigma = {
   value = 3.690,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2OH.epsilon = {
   value = 417.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.86388918E+00,
      5.59672304E-03,
      5.93271791E-06,
     -1.04532012E-08,
      4.36967278E-12,
     -3.19391367E+03,
      5.47302243E+00,
   },
   segment1 = {
      0,
      0,
      3.69266569E+00,
      8.64576797E-03,
     -3.75101120E-06,
      7.87234636E-10,
     -6.48554201E-14,
     -3.24250627E+03,
      5.81043215E+00,
   }
}
db.CH2OH.Hf = {
   value = -17800.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
