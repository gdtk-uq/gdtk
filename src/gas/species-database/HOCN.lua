db.HOCN = {}
db.HOCN.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HOCN.charge = 0
db.HOCN.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HOCN.gamma = {
   value = 1.2185e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HOCN.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HOCN.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HOCN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1368.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.78604952E+00,
      6.88667922E-03,
     -3.21487864E-06,
      5.17195767E-10,
      1.19360788E-14,
     -2.82698400E+03,
      5.63292162E+00,
   },
   segment1 = {
      0,
      0,
      5.89784885E+00,
      3.16789393E-03,
     -1.11801064E-06,
      1.77243144E-10,
     -1.04339177E-14,
     -3.70653331E+03,
     -6.18167825E+00,
   }
}
