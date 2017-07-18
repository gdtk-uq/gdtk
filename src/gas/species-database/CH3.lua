db.CH3 = {}
db.CH3.atomicConstituents = {C=1,H=3,}
db.CH3.charge = 0
db.CH3.M = {
   value = 15.0345200e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH3.gamma = {
   value = 1.276,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH3.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.67359040E+00,
          2.01095175E-03,
          5.73021856E-06,
         -6.87117425E-09,
          2.54385734E-12,
          1.64449988E+04,
          1.60456433E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.28571772E+00,
          7.23990037E-03,
         -2.98714348E-06,
          5.95684644E-10,
         -4.67154394E-14,
          1.67755843E+04,
          8.48007179E+00,
      }
   }
}
