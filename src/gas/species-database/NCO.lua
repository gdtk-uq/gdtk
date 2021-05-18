db.NCO = {}
db.NCO.atomicConstituents = {C=1,N=1,O=1,}
db.NCO.charge = 0
db.NCO.M = {
   value = 42.016800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NCO.gamma = {
   value = 1.2609e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NCO.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NCO.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.28269308E+01,
      0.88051688E-02,
     -0.83866134E-05,
      0.48016964E-08,
     -0.13313595E-11,
      0.14682477E+05,
      0.95504646E+01,
   },
   segment1 = {
      0,
      0,
      0.51521845E+01,
      0.23051761E-02,
     -0.88033153E-06,
      0.14789098E-09,
     -0.90977996E-14,
      0.14004123E+05,
     -0.25442660E+01,
   }
}
db.NCO.Hf = {
   value = 131847.241,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
