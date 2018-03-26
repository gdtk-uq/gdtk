db.C3H6 = {}
db.C3H6.atomicConstituents = {C=3,H=6,}
db.C3H6.charge = 0
db.C3H6.M = {
   value = 42.079740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H6.gamma = {
   value = 1.1474e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H6.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.C3H6.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.C3H6.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.01493307E+02,
          0.02092518E+00,
          0.04486794E-04,
         -0.01668912E-06,
          0.07158146E-10,
          0.01074826E+05,
          0.01614534E+03,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.06732257E+02,
          0.01490834E+00,
         -0.04949899E-04,
          0.07212022E-08,
         -0.03766204E-12,
         -0.09235703E+04,
         -0.01331335E+03,
      }
   }
}
