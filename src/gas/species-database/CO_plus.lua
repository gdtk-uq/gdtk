db['CO+'] = {}
db['CO+'].type = 'molecule'
db['CO+'].molecule_type = 'linear'
db['CO+'].theta_v = {
   value = 3186.0,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Capitelli et al (2005), Table 15. omega_e in ground state converted to K'
}
db['CO+'].atomicConstituents = {C=1,O=1,}
db['CO+'].charge = 1
db['CO+'].M = {
   value = 28.095514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db['CO+'].gamma = {
   value = 1.3992e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db['CO+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -2.178786658e+04,
      1.288857032e+02, 
      3.769057550e+00,
     -3.431730130e-03, 
      8.193945750e-06,
     -6.463814690e-09,
      1.803727574e-12,
      1.482345898e+05,
      3.990547070e+00
   },
   segment1 = { 
      2.316847506e+05,
     -1.057646148e+03,
      4.554257780e+00,
      4.495520320e-04,
     -2.489507047e-07,
      5.267566420e-11,
     -3.289510270e-15,
      1.555050724e+05,
     -3.873462640e+00
   },
   segment2 = {
     -3.035604054e+08,
      2.393118392e+05,
     -7.034999240e+01,
      1.139551440e-02,
     -8.315173100e-07,
      2.863705515e-11,
     -3.803269410e-16,
     -1.688617704e+06,
      6.291980420e+02
   },
}
-- No CEA transport data for CO+, just use CO
db['CO+'].ceaViscosity = db.CO.ceaViscosity 
db['CO+'].ceaThermCond = db.CO.ceaThermCond 


db['CO+'].Hf = {
   value = 1247789.204,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
