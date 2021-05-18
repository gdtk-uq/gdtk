db.HI = {}
db.HI.atomicConstituents = {I=1,H=1,}
db.HI.charge = 0
db.HI.M = {
   value = 0.12791241,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.HI.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.HI.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      1.872881730e+04,
     -3.431788840e+02,
      5.956712430e+00,
     -8.543439600e-03,
      1.454780274e-05,
     -1.049104164e-08,
      2.839734003e-12,
      3.682950720e+03,
     -8.149756090e+00,
   },
   segment1 = {
      4.724921450e+05,
     -1.923465741e+03,
      5.758048970e+00,
     -4.066266380e-04,
      9.474332050e-08,
     -1.033534431e-11,
      4.611614790e-16,
      1.394857037e+04,
     -1.182487652e+01,
   },
}
db.HI.Hf = {
   value = 26359.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
