db.NH = {}
db.NH.atomicConstituents = {N=1,H=1,}
db.NH.charge = 0
db.NH.M = {
   value = 15.0146400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NH.gamma = {
   value = 1.398,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NH.sigma = {
   value = 2.65,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH.epsilon = {
   value = 80.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.34929085E+01,
      0.31179198E-03,
     -0.14890484E-05,
      0.24816442E-08,
     -0.10356967E-11,
      0.41880629E+05,
      0.18483278E+01,
   },
   segment1 = {
      0,
      0,
      0.27836928E+01,
      0.13298430E-02,
     -0.42478047E-06,
      0.78348501E-10,
     -0.55044470E-14,
      0.42120848E+05,
      0.57407799E+01,
   }
}
db.NH.Hf = {
   value = 357032.001,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
