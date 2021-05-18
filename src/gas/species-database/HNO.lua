db.HNO = {}
db.HNO.atomicConstituents = {H=1,N=1,O=1}
db.HNO.M = {
   value = 31.01404e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO.charge = 0
db.HNO.gamma = {
   value = 1.325,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO.sigma = {
   value = 3.492,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNO.epsilon = {
   value = 116.7,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNO.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -6.854764860e+04,
      9.551627200e+02,
     -6.000720210e-01,
      7.995176750e-03,
     -6.547079160e-07,
     -3.670513400e-09,
      1.783392519e-12,
      6.435351260e+03,
      3.048166179e+01,
   },
   segment1 = {
     -5.795614980e+06,
      1.945457427e+04,
     -2.152568374e+01,
      1.797428992e-02,
     -4.976040670e-06,
      6.397924170e-10,
     -3.142619368e-14,
     -1.104192372e+05,
      1.818650338e+02,
   }
}
db.HNO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.45334916E+01,
     -0.56696171E-02,
      0.18473207E-04,
     -0.17137094E-07,
      0.55454573E-11,
      0.11548297E+05,
      0.17498417E+01,
   },
   segment1 = {
      0,
      0,
      0.29792509E+01,
      0.34944059E-02,
     -0.78549778E-06,
      0.57479594E-10,
     -0.19335916E-15,
      0.11750582E+05,
      0.86063728E+01,
   }
}
db.HNO.Hf = {
   value = 102032.725,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
