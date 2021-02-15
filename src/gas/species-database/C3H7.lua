db.C3H7 = {}
db.C3H7.atomicConstituents = {C=3,H=7,}
db.C3H7.charge = 0
db.C3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H7.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H7.epsilon = {
   value = 266.8,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H7.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.10515518E+01,
      0.25991980E-01,
      0.23800540E-05,
     -0.19609569E-07,
      0.93732470E-11,
      0.10631863E+05,
      0.21122559E+02,
   },
   segment1 = {
      0,
      0,
      0.77026987E+01,
      0.16044203E-01,
     -0.52833220E-05,
      0.76298590E-09,
     -0.39392284E-13,
      0.82984336E+04,
     -0.15480180E+02,
   }
}
