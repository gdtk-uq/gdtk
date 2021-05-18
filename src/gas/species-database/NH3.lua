db.NH3 = {}
db.NH3.atomicConstituents = {N=1,H=3,}
db.NH3.charge = 0
db.NH3.M = {
   value = 17.030520e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH3.gamma = {
   value = 1.3036e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH3.sigma = {
   value = 2.92,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH3.epsilon = {
   value = 481.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.42860274E+01,
     -0.46605230E-02,
      0.21718513E-04,
     -0.22808887E-07,
      0.82638046E-11,
     -0.67417285E+04,
     -0.62537277E+00,
   },
   segment1 = {
      0,
      0,
      0.26344521E+01,
      0.56662560E-02,
     -0.17278676E-05,
      0.23867161E-09,
     -0.12578786E-13,
     -0.65446958E+04,
      0.65662928E+01,
   }
}
db.NH3.Hf = {
   value = -45940.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
