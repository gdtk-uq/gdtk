db.H2O = {}
db.H2O.atomicConstituents = {O=1,H=2,}
db.H2O.charge = 0
db.H2O.M = {
   value = 0.01801528,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2O.gamma = {
   value = 1.32900000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2O.sigma = {
   value = 2.605,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O.epsilon = {
   value = 572.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O.Lewis = {
   value = 0.854
}
db.H2O.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -3.947960830e+04,
      5.755731020e+02,
      9.317826530e-01,
      7.222712860e-03,
     -7.342557370e-06,
      4.955043490e-09,
     -1.336933246e-12,
     -3.303974310e+04,
      1.724205775e+01,
   },
   segment1 = {
      1.034972096e+06,
     -2.412698562e+03,
      4.646110780e+00,
      2.291998307e-03,
     -6.836830480e-07,
      9.426468930e-11,
     -4.822380530e-15,
     -1.384286509e+04,
     -7.978148510e+00,
   },
}
db.H2O.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  5.0019557e-01,
      B = -6.9712796e+02,
      C =  8.8163892e+04,
      D =  3.0836508e+00
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  5.8988538e-01,
      B = -5.3769814e+02,
      C =  5.4263513e+04,
      D =  2.3386375e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4330087e-01,
      B = -9.5668913e+01,
      C = -3.7742283e+05,
      D =  1.8125190e+00
   },
}
db.H2O.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  1.0966389e+00,
      B = -5.5513429e+02,
      C =  1.0623408e+05,
      D = -2.4664550e-01
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  3.9367933e-01,
      B = -2.2524226e+03,
      C =  6.1217458e+05,
      D =  5.8011317e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.1858737e-01,
      B = -1.4096649e+04,
      C =  1.9179190e+07,
      D =  1.4345613e+01
   },
}
db.H2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.19864056E+00,
     -2.03643410E-03,
      6.52040211E-06,
     -5.48797062E-09,
      1.77197817E-12,
     -3.02937267E+04,
     -8.49032208E-01,
   },
   segment1 = {
      0,
      0,
      3.03399249E+00,
      2.17691804E-03,
     -1.64072518E-07,
     -9.70419870E-11,
      1.68200992E-14,
     -3.00042971E+04,
      4.96677010E+00,
   }
}

db.H2O.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.596254204075e+01,
      B = 1.678935618571e-01,
      C = 1.929065548334e-01,
      D = -1.356763840380e-02,
   }
}
db.H2O.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = 5.416062110310e+00,
      B = -5.968633913268e+00,
      C = 1.091216246422e+00,
      D = -5.536357184347e-02,
   }
}

db.H2O.Hf = {
   value = -241826.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
