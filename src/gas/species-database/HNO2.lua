db.HNO2 = {}
db.HNO2.atomicConstituents = {H=1,N=1,O=2}
db.HNO2.charge = 0
db.HNO2.M = {
   value = 0.04701344,
   unit = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO2.gamma = {
   value = 1.218,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      8.591985060e+03,
      1.203644046e+02,
      9.412979120e-01,
      1.942891839e-02,
     -2.253174194e-05,
      1.384587594e-08,
     -3.473550460e-12,
     -1.106337202e+04,
      2.073967331e+01,
   },
   segment1 = {
      8.787904130e+05,
     -3.990455030e+03,
      1.187349269e+01,
     -4.881900610e-04,
      7.133636790e-08,
     -5.376303340e-12,
      1.581778986e-16,
      1.246343241e+04,
     -4.608874688e+01,
   }
}
db.HNO2.Hf = {
   value = -78451.922,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
