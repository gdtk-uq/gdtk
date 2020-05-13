db['N+'] = {}
db['N+'].atomicConstituents = {N=1}
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
db['N+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {1000.0, 1000.0},
   segment0 = {
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
   segment1 = {
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
   segment2 = {
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


