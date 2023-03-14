db.N = {}
db.N.type = "atom"
db.N.atomicConstituents = {N=1}
db.N.charge = 0
db.N.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   reference = 'CEA2::thermo.inp'
}
db.N.gamma =  {
   value = 5./3.,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.N.sigma = {
   value = 3.298,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Svehla (1962)'
}
db.N.epsilon = {
   value = 71.4, 
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Svehla (1962)'
}
db.N.electronic_levels = {
  Te = {
    value = {0.0, 19227.95, 28839.18, 83335.6, 86192.79, 88132.45, 93581.55},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {4, 10, 6, 12, 6, 12, 2},
    units = 'NA',
    description = 'degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db.N.ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      5.61061063e+04,
      4.19425139e+00
   },
   segment1 = {
     -2.27073277e+05,
      8.14052944e+02,
      1.32705137e+00,
      8.62721731e-04,
     -3.35747089e-07,
      6.29010687e-11,
     -3.90674587e-15,
      5.10943141e+04,
      1.22823713e+01
   },
   segment2 = {
     -2.04738994e+09,
      1.45842847e+06,
     -4.18833824e+02,
      6.25994407e-02,
     -4.96542822e-06,
      1.98247052e-10,
     -3.05470194e-15,
     -1.12727730e+07,
      3.58487417e+03
   },
   segment3 = {
      5.74291902e+11,
     -1.29039294e+08,
      1.15381467e+04,
     -5.25078568e-01,
      1.29219090e-05,
     -1.63974231e-10,
      8.41878585e-16,
      1.15261836e+09,
     -1.11649232e+05
   },
   reference="Johnston et al. AIAA Paper 2015-3110"
}
db.N.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.25000000E+01,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
      0.56104637E+05,
      0.41939087E+01,
   },
   segment1 = {
      0,
      0,
      0.24159429E+01,
      0.17489065E-03,
     -0.11902369E-06,
      0.30226245E-10,
     -0.20360982E-14,
      0.56133773E+05,
      0.46496096E+01,
   }
}
db.N.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83724737e+00,
      B = 0.43997150e+03,
      C = -0.17450753e+06,
      D = 0.10365689e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.89986588e+00,
      B = 0.14112801e+04,
      C = -0.18200478e+07,
      D = -0.55811716e+00
   },
   reference = 'from CEA2::trans.inp which cites Levin et al (1990)'
}
db.N.ceaThermCond = {
   nsegments = 2, 
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83771661e+00,
      B = 0.44243270e+03,
      C = -0.17578446e+06,
      D = 0.89942915e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.90001710e+00,
      B = 0.14141175e+04,
      C = -0.18262403e+07,
      D = 0.24048513e+00
   },
   reference = 'from CEA2::trans.inp which cites Levin et al (1990)'
}

db.N.Hf = {
   value = 472680.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
