db.CH3OH = {}
db.CH3OH.atomicConstituents = {C=1,H=4,O=1,}
db.CH3OH.charge = 0
db.CH3OH.M = {
   value = 32.041860e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3OH.gamma = {
   value = 1.2320e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3OH.sigma = {
   value = 3.626,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3OH.epsilon = {
   value = 481.8,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      5.71539582E+00,
     -1.52309129E-02,
      6.52441155E-05,
     -7.10806889E-08,
      2.61352698E-11,
     -2.56427656E+04,
     -1.50409823E+00,
   },
   segment1 = {
      0,
      0,
      1.78970791E+00,
      1.40938292E-02,
     -6.36500835E-06,
      1.38171085E-09,
     -1.17060220E-13,
     -2.53748747E+04,
      1.45023623E+01,
   }
}
db.CH3OH.Hf = {
   value = -200940.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
