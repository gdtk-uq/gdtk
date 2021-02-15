db.NNH = {}
db.NNH.atomicConstituents = {N=2,H=1,}
db.NNH.charge = 0
db.NNH.M = {
   value = 29.021340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NNH.gamma = {
   value = 1.3152e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NNH.sigma = {
   value = 3.798,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NNH.epsilon = {
   value = 71.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NNH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.43446927E+01,
     -0.48497072E-02,
      0.20059459E-04,
     -0.21726464E-07,
      0.79469539E-11,
      0.28791973E+05,
      0.29779410E+01,
   },
   segment1 = {
      0,
      0,
      0.37667544E+01,
      0.28915082E-02,
     -0.10416620E-05,
      0.16842594E-09,
     -0.10091896E-13,
      0.28650697E+05,
      0.44705067E+01,
   }
}
