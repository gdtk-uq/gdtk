db.HCCOH = {}
db.HCCOH.atomicConstituents = {C=2,H=2,O=1,}
db.HCCOH.charge = 0
db.HCCOH.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCCOH.gamma = {
   value = 1.1656e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCCOH.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCOH.epsilon = {
   value = 436.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCOH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.12423733E+01,
      0.31072201E-01,
     -0.50866864E-04,
      0.43137131E-07,
     -0.14014594E-10,
      0.80316143E+04,
      0.13874319E+02,
   },
   segment1 = {
      0,
      0,
      0.59238291E+01,
      0.67923600E-02,
     -0.25658564E-05,
      0.44987841E-09,
     -0.29940101E-13,
      0.72646260E+04,
     -0.76017742E+01,
   }
}
