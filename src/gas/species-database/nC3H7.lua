db.nC3H7 = {}
db.nC3H7.atomicConstituents = {C=3,H=7,}
db.nC3H7.charge = 0
db.nC3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.nC3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.nC3H7.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.nC3H7.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.nC3H7.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.10491173E+01,
          0.26008973E-01,
          0.23542516E-05,
         -0.19595132E-07,
          0.93720207E-11,
          0.10312346E+05,
          0.21136034E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3000.0,
      coeffs = {
         0,
         0,
          0.77097479E+01,
          0.16031485E-01,
         -0.52720238E-05,
          0.75888352E-09,
         -0.38862719E-13,
          0.79762236E+04,
         -0.15515297E+02,
      }
   }
}
