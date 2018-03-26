db.aC3H5 = {}
db.aC3H5.atomicConstituents = {C=3,H=5,}
db.aC3H5.charge = 0
db.aC3H5.M = {
   value = 41.071800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.aC3H5.gamma = {
   value = 1.1502e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.aC3H5.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.aC3H5.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.aC3H5.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.13631835E+01,
          0.19813821E-01,
          0.12497060E-04,
         -0.33355555E-07,
          0.15846571E-10,
          0.19245629E+05,
          0.17173214E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3000.0,
      coeffs = {
         0,
         0,
          0.65007877E+01,
          0.14324731E-01,
         -0.56781632E-05,
          0.11080801E-08,
         -0.90363887E-13,
          0.17482449E+05,
         -0.11243050E+02,
      }
   }
}
