db.C3H8 = {}
db.C3H8.atomicConstituents = {C=3,H=8,}
db.C3H8.charge = 0
db.C3H8.M = {
   value = 44.095620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H8.gamma = {
   value = 1.1267e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H8.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H8.epsilon = {
   value = 266.80,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H8.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.93355381E+00,
      0.26424579E-01,
      0.61059727E-05,
     -0.21977499E-07,
      0.95149253E-11,
     -0.13958520E+05,
      0.19201691E+02,
   },
   segment1 = {
      0,
      0,
      0.75341368E+01,
      0.18872239E-01,
     -0.62718491E-05,
      0.91475649E-09,
     -0.47838069E-13,
     -0.16467516E+05,
     -0.17892349E+02,
   }
}
db.C3H8.Hf = {
   value = -104680.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
