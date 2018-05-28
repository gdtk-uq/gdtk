db.H2CC = {}
db.H2CC.atomicConstituents = {H=2,C=2,}
db.H2CC.charge = 0
db.H2CC.M = {
   value = 26.0372800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2CC.gamma = {
   value = 1.242,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2CC.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file'
}
db.H2CC.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file'
}
db.H2CC.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.32815483E+01,
          0.69764791E-02,
         -0.23855244E-05,
         -0.12104432E-08,
          0.98189545E-12,
          0.48621794E+05,
          0.59203910E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.42780340E+01,
          0.47562804E-02,
         -0.16301009E-05,
          0.25462806E-09,
         -0.14886379E-13,
          0.48316688E+05,
          0.64023701E+00,
      }
   }
}


