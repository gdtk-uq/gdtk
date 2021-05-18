db.I2 = {}
db.I2.atomicConstituents = {I=2,}
db.I2.charge = 0
db.I2.M = {
   value = 0.25380894,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.I2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.I2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -5.087968770e+03,
     -1.249585210e+01,
      4.504219090e+00,
      1.370962533e-04,
     -1.390523014e-07,
      1.174813853e-10,
     -2.337541043e-14,
      6.213469810e+03,
      5.583836940e+00,
   },
   segment1 = {
     -5.632594160e+06,
      1.793961560e+04,
     -1.723055169e+01,
      1.244214080e-02,
     -3.332768580e-06,
      4.125477940e-10,
     -1.960461713e-14,
     -1.068505292e+05,
      1.600531883e+02,
   },
}
db.I2.Hf = {
   value = 62420.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
