db.NO2 = {}
db.NO2.atomicConstituents = {O=2,N=1}
db.NO2.charge = 0
db.NO2.M = {
   value = 0.0460055000,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NO2.gamma = {
   value = 1.287,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NO2.sigma = {
   value = 3.5,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NO2.epsilon = {
   value = 200.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NO2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -5.642038780e+04,
      9.633085720e+02,
     -2.434510974e+00,
      1.927760886e-02,
     -1.874559328e-05,
      9.145497730e-09,
     -1.777647635e-12,
     -1.547925037e+03,
      4.067851210e+01,
   },
   segment1 = {
      7.213001570e+05,
     -3.832615200e+03,
      1.113963285e+01,
     -2.238062246e-03,
      6.547723430e-07,
     -7.611335900e-11,
      3.328361050e-15,
      2.502497403e+04,
     -4.305130040e+01,
   }
}
db.NO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.39440312E+01,
     -0.15854290E-02,
      0.16657812E-04,
     -0.20475426E-07,
      0.78350564E-11,
      0.28966179E+04,
      0.63119917E+01,
   },
   segment1 = {
      0,
      0,
      0.48847542E+01,
      0.21723956E-02,
     -0.82806906E-06,
      0.15747510E-09,
     -0.10510895E-13,
      0.23164983E+04,
     -0.11741695E+00,
   }
}
db.NO2.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A =  0.57379100e+00,
      B = -0.12636034e+03,
      C =  0.21566823e+04,
      D =  0.22287492e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.64239645e+00,
      B =  0.60012144e+00,
      C = -0.27020876e+05,
      D = 0.16570566e+01
   }
}
db.NO2.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A =  0.48574998e+00,
      B = -0.50702110e+03,
      C =  0.46605820e+05,
      D =  0.36444556e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.97660465e+00,
      B =  0.72760751e+03,
      C = -0.32527989e+06,
      D = -0.60899123e+00
   }
}
db.NO2.Hf = {
   value = 34193.019,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
