db.HNCO = {}
db.HNCO.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HNCO.charge = 0
db.HNCO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HNCO.gamma = {
   value = 1.2173e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HNCO.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNCO.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1478.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.63096317E+00,
      7.30282357E-03,
     -2.28050003E-06,
     -6.61271298E-10,
      3.62235752E-13,
     -1.55873636E+04,
      6.19457727E+00,
   },
   segment1 = {
      0,
      0,
      6.22395134E+00,
      3.17864004E-03,
     -1.09378755E-06,
      1.70735163E-10,
     -9.95021955E-15,
     -1.66599344E+04,
     -8.38224741E+00,
   }
}
db.HNCO.Hf = {
   value = -118056.529,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
