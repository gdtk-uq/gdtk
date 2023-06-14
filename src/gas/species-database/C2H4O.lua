db.C2H4O = {}
db.C2H4O.atomicConstituents = {C=2,H=4,O=1}
db.C2H4O.charge = 0
db.C2H4O.M = {
   value = 44.0525600e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2 thermo.inp'
}
db.C2H4O.gamma = {
   value = 1.1762e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'taken from isomer: CH3CHO'
}
db.C2H4O.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'taken from isomer: CH3CHO'
}
db.C2H4O.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'taken from isomer: CH3CHO'
}
db.C2H4O.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.728233345e+05,
      3.816678800e+03,
     -2.629851977e+01,
      1.014103162e-01,
     -1.240578373e-04,
      8.034040350e-08,
     -2.120942544e-11,
     -2.437519333e+04,
      1.654885056e+02
   },
   segment1 = {
      3.151809957e+06,
     -1.423646316e+04,
      2.708080476e+01,
     -2.606238456e-03,
      4.853891930e-07,
     -4.852144760e-11,
      2.011778721e-15,
      7.662561440e+04,
     -1.563952401e+02
   }
}

