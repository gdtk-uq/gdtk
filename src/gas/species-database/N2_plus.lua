db['N2+'] = {}
db['N2+'].type = "molecule"
db['N2+'].molecule_type = "linear"
db['N2+'].theta_v = {
   value = 3393.44,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Copied from N2'
}
db['N2+'].electronic_levels = {
  Te = {
    value = {0.0, 9167.46, 25461.11, 51663.2, 64609.03},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {2, 4, 2, 4, 2},
    units = 'cm^(-1)',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST',
  }
}
db['N2+'].atomicConstituents = {N=2}
db['N2+'].charge = 1
db['N2+'].M = {
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db['N2+'].gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

db['N2+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -3.474047470e+04,
      2.696222703e+02,
      3.164916370e+00,
     -2.132239781e-03,
      6.730476400e-06,
     -5.637304970e-09,
      1.621756000e-12,
      1.790004424e+05,
      6.832974166e+00
   },
   segment1 = {
     -2.845599002e+06,
      7.058893030e+03,
     -2.884886385e+00,
      3.068677059e-03,
     -4.361652310e-07,
      2.102514545e-11, 
      5.411996470e-16,
      1.340388483e+05,
      5.090897022e+01
   },
   segment2 = {
     -3.712829770e+08,
      3.139287234e+05,
     -9.603518050e+01, 
      1.571193286e-02,
     -1.175065525e-06,
      4.144441230e-11,
     -5.621893090e-16,
     -2.217361867e+06,
      8.436270947e+02
   },
   ref="from CEA2::thermo.inp"
}
-- No CEA transport data for N2+, just use N2
db['N2+'].ceaViscosity = db.N2.ceaViscosity
db['N2+'].ceaThermCond = db.N2.ceaThermCond



db['N2+'].Hf = {
   value = 1509508.424,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
