db.Kr = {}
db.Kr.atomicConstituents = {Kr=1,}
db.Kr.charge = 0
db.Kr.M = {
   value = 0.0838,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.Kr.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Kr.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
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
      5.490956510e+00
   },
   segment1 = {
      2.643639057e+02,
     -7.910050820e-01,
      2.500920585e+00,
     -5.328164110e-07,
      1.620730161e-10,
     -2.467898017e-14,
      1.478585040e-18,
     -7.403488940e+02,
      5.484398150e+00
   },
   segment2 = {
     -1.375531087e+09,
      9.064030530e+05,
     -2.403481435e+02,
      3.378312030e-02,
     -2.563103877e-06,
      9.969787790e-11,
     -1.521249677e-15,
     -7.111667370e+06,
      2.086866326e+03
   },
}
db.Kr.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58597795e+00,
      B = -0.12924832e+03,
      C =  0.47495759e+04,
      D =  0.25793650e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.68985845e+00,
      B =  0.56296306e+02,
      C = -0.36082600e+05,
      D =  0.17170715e+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  0.76582939e+00,
      B =  0.68610377e+03,
      C = -0.82537190e+06,
      D =  0.97565128e+00
   },
}
db.Kr.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58008139e+00,
      B = -0.13792556e+03,
      C =  0.60771460e+04,
      D =  0.16420039e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.68859431e+00,
      B =  0.51765647e+02,
      C = -0.34512131e+05,
      D =  0.74332130e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  0.76365443e+00,
      B =  0.65175847e+03,
      C = -0.73589800e+06,
      D =  0.12112126e-01
   },
}
db.Kr.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
