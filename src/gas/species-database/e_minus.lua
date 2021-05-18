db['e-'] = {}
db['e-'].type = "electron"
db['e-'].atomicConstituents = {}
db['e-'].charge = -1
db['e-'].M = {
   value = 0.000548579903e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['e-'].gamma = {
   value = 1.667, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}

db['e-'].ceaThermoCoeffs = {
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
     -7.453750000e+02,
     -1.172081224e+01
   },
   segment1 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
     -7.453750000e+02,
     -1.172081224e+01
   },
   segment2 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
     -7.453750000e+02,
     -1.172081224e+01
   },
}


db['e-'].Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
