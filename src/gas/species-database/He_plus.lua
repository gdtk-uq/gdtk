-- Helium_plus ported from Daniel's He_plus.lua file in the cfcfd3 collection
-- Yu Liu, 2018-06-27
db['He+'] = {}
db['He+'].atomicConstituents = {He=1,}
db['He+'].charge = 1
db['He+'].M = {
   value = 4.0020534e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db['He+'].gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['He+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      0.000000000e+00,  
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,  
      0.000000000e+00,
      0.000000000e+00, 
      2.853233739e+05, 
      1.621665557e+00
   },
   segment1 = {
      0.000000000e+00,  
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,  
      0.000000000e+00,
      0.000000000e+00, 
      2.853233739e+05, 
      1.621665557e+00
   },
   segment2 = {
      0.000000000e+00,  
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,  
      0.000000000e+00,
      0.000000000e+00, 
      2.853233739e+05, 
      1.621665557e+00
   },
   reference = 'cea2::thermo.inp'
}
-- No CEA transport data for He+, just use He
db['He+'].ceaViscosity = db.He.ceaViscosity 
db['He+'].ceaThermCond = db.He.ceaThermCond 
db['He+'].Hf = {
   value = 2378521.473,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
