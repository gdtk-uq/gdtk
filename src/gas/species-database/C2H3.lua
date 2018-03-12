db.C2H3 = {}
db.C2H3.atomicConstituents = {C=2,H=3,}
db.C2H3.charge = 0
db.C2H3.M = {
   value = 27.045220e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H3.gamma = {
   value = 1.2408e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H3.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H3.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.21246645E+00,
          1.51479162E-03,
          2.59209412E-05,
         -3.57657847E-08,
          1.47150873E-11,
          3.48598468E+04,
          8.51054025E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.01672400E+00,
          1.03302292E-02,
         -4.68082349E-06,
          1.01763288E-09,
         -8.62607041E-14,
          3.46128739E+04,
          7.78732378E+00,
      }
   }
}
