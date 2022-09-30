db['N+'] = {}
db['N+'].atomicConstituents = {N=1}
db['N+'].type = "atom"
db['N+'].charge = 1
db['N+'].M = {
   value = 14.0061514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['N+'].gamma =  {
   value = 1.641, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
db['N+'].electronic_levels = {
  Te = {
    value = {0.0, 70.07, 188.19, 22036.47, 47031.63, 67312.23},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {1, 3, 5, 5, 1, 5},
    units = 'NA',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db['N+'].ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      5.23707921e+03,
      2.29995832e+00,
      2.48237856e+00,
      2.73749076e-05,
     -3.13444758e-08,
      1.85011133e-11,
     -4.44735098e-15,
      2.25626950e+05,
      5.04771460e+00
   },
   segment1 = {
      3.51764922e+05,
     -1.04468568e+03,
      3.68498994e+00,
     -6.41320819e-04,
      1.66764737e-07,
     -1.82509101e-11,
      7.37893962e-16,
      2.32355850e+05,
     -3.48540450e+00
   },
   segment2 = {
      2.08931831e+07,
     -1.36213868e+04,
      5.54565345e+00,
     -2.68374643e-04,
      1.46878484e-08,
     -4.22829994e-13,
      5.70569781e-18,
      3.33626095e+05,
     -2.20294391e+01
   },
   segment3 = {
      1.48351337e+09,
     -2.69100020e+05,
      1.91340133e+01,
     -2.13380850e-04,
     -1.51589694e-08,
      5.96512853e-13,
     -5.46224064e-18,
      2.70008612e+06,
     -1.64949269e+02
   },
   reference="Johnston et al. AIAA Paper 2015-3110"
}
-- No CEA transport data for N+, just use N
db['N+'].ceaViscosity = db.N.ceaViscosity 
db['N+'].ceaThermCond = db.N.ceaThermCond 


db['N+'].Hf = {
   value = 1882127.624,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
