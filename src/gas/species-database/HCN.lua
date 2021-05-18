db.HCN = {}
db.HCN.atomicConstituents = {C=1,H=1,N=1,}
db.HCN.charge = 0
db.HCN.M = {
   value = 27.0253400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HCN.gamma = {
   value = 1.301,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HCN.sigma = {
   value = 3.63,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCN.epsilon = {
   value = 569.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.22589886E+01,
      0.10051170E-01,
     -0.13351763E-04,
      0.10092349E-07,
     -0.30089028E-11,
      0.14712633E+05,
      0.89164419E+01,
   },
   segment1 = {
      0,
      0,
      0.38022392E+01,
      0.31464228E-02,
     -0.10632185E-05,
      0.16619757E-09,
     -0.97997570E-14,
      0.14407292E+05,
      0.15754601E+01,
   }
}
db.HCN.Hf = {
   value = 133082.46,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
