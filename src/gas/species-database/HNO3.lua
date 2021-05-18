db.HNO3 = {}
db.HNO3.atomicConstituents = {H=1,N=1,O=3}
db.HNO3.M = {
   value = 63.01284e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO3.charge = 0
db.HNO3.gamma = {
   value = 1.181,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO3.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      9.202869010e+03,
      1.093774496e+02,
     -4.521042450e-01,
      2.984914503e-02,
     -3.190635500e-05,
      1.720931528e-08,
     -3.782649830e-12,
     -1.764048507e+04,
      2.746644879e+01,
   },
   segment1 = {
     -9.497809640e+04,
     -2.733024468e+03,
      1.449426995e+01,
     -7.821868050e-04,
      1.702693665e-07,
     -1.930543961e-11,
      8.870455120e-16,
     -4.882517780e+03,
     -5.928392985e+01,
   }
}
db.HNO3.Hf = {
   value = -133912.869,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
