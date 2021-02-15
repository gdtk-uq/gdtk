db.HCNO = {}
db.HCNO.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HCNO.charge = 0
db.HCNO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCNO.gamma = {
   value = 1.2154e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCNO.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCNO.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCNO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1382.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.64727989E+00,
      1.27505342E-02,
     -1.04794236E-05,
      4.41432836E-09,
     -7.57521466E-13,
      1.92990252E+04,
      1.07332972E+01,
   },
   segment1 = {
      0,
      0,
      6.59860456E+00,
      3.02778626E-03,
     -1.07704346E-06,
      1.71666528E-10,
     -1.01439391E-14,
      1.79661339E+04,
     -1.03306599E+01,
   }
}
