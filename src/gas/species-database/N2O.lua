db.N2O = {}
db.N2O.atomicConstituents = {N=2,O=1,}
db.N2O.charge = 0
db.N2O.M = {
   value = 44.012800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.N2O.gamma = {
   value = 1.2735e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.N2O.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2O.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2O.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      4.288225970e+04,
     -6.440118440e+02,
      6.034351430e+00,
      2.265394436e-04,
      3.472782850e-06,
     -3.627748640e-09,
      1.137969552e-12,
      1.179405506e+04,
     -1.003128570e+01,
   },
   segment1 = {
      3.438448040e+05,
     -2.404557558e+03,
      9.125636220e+00,
     -5.401667930e-04,
      1.315124031e-07,
     -1.414215100e-11,
      6.381066870e-16,
      2.198632638e+04,
     -3.147805016e+01,
   }
}
db.N2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.22571502E+01,
      0.11304728E-01,
     -0.13671319E-04,
      0.96819806E-08,
     -0.29307182E-11,
      0.87417744E+04,
      0.10757992E+02,
   },
   segment1 = {
      0,
      0,
      0.48230729E+01,
      0.26270251E-02,
     -0.95850874E-06,
      0.16000712E-09,
     -0.97752303E-14,
      0.80734048E+04,
     -0.22017207E+01,
   }
}
db.N2O.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58959112e+00,
      B = -0.15565178e+03,
      C =  0.37630431e+04,
      D =  0.21223853e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.64571469e+00,
      B = -0.88806585e+01,
      C = -0.41560559e+05,
      D = 0.16332498e+01
   }
}
db.N2O.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.65165376e+00,
      B = -0.34373058e+03,
      C =  0.15090399e+05,
      D =  0.24242359e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.64720604e+00,
      B = -0.78680195e+02,
      C = -0.11965729e+06,
      D = 0.23246569e+01
   }
}
db.N2O.Hf = {
   value = 81600.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}