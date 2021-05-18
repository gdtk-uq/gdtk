db.NH2 = {}
db.NH2.atomicConstituents = {N=1,H=2,}
db.NH2.charge = 0
db.NH2.M = {
   value = 16.022580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH2.gamma = {
   value = 1.3254e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH2.sigma = {
   value = 2.65,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH2.epsilon = {
   value = 80.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.42040029E+01,
     -0.21061385E-02,
      0.71068348E-05,
     -0.56115197E-08,
      0.16440717E-11,
      0.21885910E+05,
     -0.14184248E+00,
   },
   segment1 = {
      0,
      0,
      0.28347421E+01,
      0.32073082E-02,
     -0.93390804E-06,
      0.13702953E-09,
     -0.79206144E-14,
      0.22171957E+05,
      0.65204163E+01,
   }
}
db.NH2.Hf = {
   value = 189134.713,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
