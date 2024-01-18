db['C+'] = {}
db['C+'].type = 'atom'
db['C+'].atomicConstituents = {C=1,}
db['C+'].charge = 1
db['C+'].M = { 
   value = 12.0101514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db['C+'].gamma = { 
   value = 1.664,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'copied from C.lua'
}
db['C+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      2.258535929e+03,
     -1.574575687e+00,
      2.503637730e+00,
     -5.202878370e-06,
      4.516908390e-09,
     -2.181431053e-12,
      4.495047033e-16,
      2.168951913e+05,
      4.345699505e+00
   },
   segment1 = {
      1.255112551e+04,
     -3.411874670e+01,
      2.543383218e+00,
     -2.805120849e-05,
      9.751641970e-09,
     -1.736855394e-12,
      1.246191931e-16,
      2.171001786e+05, 
      4.063913515e+00
   },
   segment2 = {
      5.618135320e+05,
     -6.047058900e+03,
      5.884541470e+00,
     -7.211894530e-04,
      6.823484110e-08,
     -2.599878590e-12,
      3.633868358e-17,
      2.581370458e+05,
     -2.280019759e+01
   }
}
db['C+'].Hf = {
   value = 1809444.482,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
