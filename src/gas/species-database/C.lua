db.C = {}
db.C.atomicConstituents = {C=1,}
db.C.type = "atom"
db.C.charge = 0
db.C.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C.gamma = { 
   value = 1.664,
   units = 'non-dimensional',
    description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.C.sigma = {
   value = 3.298,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C.epsilon = {
   value = 71.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      6.495031470e+02,
     -9.649010860e-01,
      2.504675479e+00,
     -1.281448025e-05,
      1.980133654e-08,
     -1.606144025e-11,
      5.314483411e-15,
      8.545763110e+04,
      4.747924288e+00
   },
   segment1 = {
     -1.289136472e+05,
      1.719528572e+02,
      2.646044387e+00,
     -3.353068950e-04,
      1.742092740e-07,
     -2.902817829e-11,
      1.642182385e-15,
      8.410597850e+04,
      4.130047418e+00
   },
   segment2 = {
      4.432528010e+08,
     -2.886018412e+05,
      7.737108320e+01,
     -9.715281890e-03,
      6.649595330e-07,
     -2.230078776e-11,
      2.899388702e-16,
      2.355273444e+06,
     -6.405123160e+02
   }
}
db.C.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0,
      0,
      2.55423955E+00,
     -3.21537724E-04,
      7.33792245E-07,
     -7.32234889E-10,
      2.66521446E-13,
      8.54438832E+04,
      4.53130848E+00,
   },
   segment1 = {
      0,
      0,
      2.49266888E+00,
      4.79889284E-05,
     -7.24335020E-08,
      3.74291029E-11,
     -4.87277893E-15,
      8.54512953E+04,
      4.80150373E+00,
   }
}
db.C.Hf = {
   value = 716680.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
