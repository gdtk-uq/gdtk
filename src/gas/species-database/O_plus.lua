db['O+'] = {}
db['O+'].type = "atom"
db['O+'].atomicConstituents = {O=1,}
db['O+'].charge = 1
db['O+'].M = {
   value = 15.9988514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['O+'].gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['O+'].electronic_levels = {
  Te = {
    value = {0.0, 91.25, 61871.81, 61903.47, 61944.19, 107807.11},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {2, 4, 2, 4, 6, 6,},
    units = 'NA',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db['O+'].ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      0.00000000e+00,
      0.00000000e+00,
      2.50000016e+00,
      0.00000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      1.87935284e+05,
      4.39337767e+00
   },
   segment1 = {
     -3.14817544e+05,
      9.52751354e+02,
      1.37960003e+00,
      6.50225275e-04,
     -1.93781966e-07,
      2.73210815e-11,
     -1.29823418e-15,
      1.81897151e+05,
      1.23661521e+01
   },
   segment2 = {
     -1.89172344e+08,
      1.32258091e+05,
     -3.33573784e+01,
      4.61791665e-03,
     -2.81231154e-07,
      8.25139893e-12,
     -9.46296238e-17,
     -8.44265776e+05,
      3.12572886e+02
   },
   segment3 = {
     -1.25053902e+09,
      1.95862039e+05,
     -7.04252592e+00,
      2.94936789e-04,
     -7.50877176e-09,
      1.46522532e-13,
     -1.07948513e-18,
     -1.65453653e+06,
      1.03404000e+02
   },
   reference="Johnston et al. AIAA Paper 2015-3110"
}
-- No CEA transport data, just use O
db['O+'].ceaViscosity = db.O.ceaViscosity
db['O+'].ceaThermCond = db.O.ceaThermCond 

db['O+'].Hf = {
   value = 1568787.228,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
