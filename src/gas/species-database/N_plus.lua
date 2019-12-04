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
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      5.237079210e+03,
      2.299958315e+00,
      2.487488821e+00,
      2.737490756e-05,
     -3.134447576e-08,
      1.850111332e-11,
     -4.447350984e-15,
      2.256284738e+05,
      5.076830786e+00
   },
   segment1 = {
      2.904970374e+05,
     -8.557908610e+02,
      3.477389290e+00,
     -5.288267190e-04,
      1.352350307e-07,
     -1.389834122e-11,
      5.046166279e-16,
      2.310809984e+05,
     -1.994146545e+00
   },
   segment2 = {
      1.646092148e+07,
     -1.113165218e+04,
      4.976986640e+00,
     -2.005393583e-04,
      1.022481356e-08,
     -2.691430863e-13,
      3.539931593e-18,
      3.136284696e+05,
     -1.706646380e+01
   },
   reference="from CEA2::thermo.inp"
}
-- No CEA transport data for N+, just use N
db['N+'].ceaViscosity = db.N.ceaViscosity 
db['N+'].ceaThermCond = db.N.ceaThermCond 


