db.H2CN = {}
db.H2CN.atomicConstituents = {C=1,H=2,N=1,}
db.H2CN.charge = 0
db.H2CN.M = {
   value = 28.033280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.H2CN.gamma = {
   value = 1.2769e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.H2CN.sigma = {
   value = 3.63,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2CN.epsilon = {
   value = 569.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2CN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   T_break_points = {300.0, 1000.0, 4000.0},
   T_blend_ranges = {400.0},
   nsegments = 2, 
   segment0 ={
      0,
      0,
      0.28516610E+01,
      0.56952331E-02,
      0.10711400E-05,
     -0.16226120E-08,
     -0.23511081E-12,
      0.28637820E+05,
      0.89927511E+01,
   },
   segment1 = {
      0,
      0,
      0.52097030E+01,
      0.29692911E-02,
     -0.28555891E-06,
     -0.16355500E-09,
      0.30432589E-13,
      0.27677109E+05,
     -0.44444780E+01,
   }
}
