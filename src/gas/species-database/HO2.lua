db.HO2 = {}
db.HO2.atomicConstituents = {O=2,H=1,}
db.HO2.charge = 0
db.HO2.M = {
   value = 0.03300674,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HO2.gamma = {
   value = 1.31200000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HO2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -7.598882540e+04,
          1.329383918e+03,
         -4.677388240e+00,
          2.508308202e-02,
         -3.006551588e-05,
          1.895600056e-08,
         -4.828567390e-12,
         -5.873350960e+03,
          5.193602140e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         -1.810669724e+06,
          4.963192030e+03,
         -1.039498992e+00,
          4.560148530e-03,
         -1.061859447e-06,
          1.144567878e-10,
         -4.763064160e-15,
         -3.200817190e+04,
          4.066850920e+01,
      }
   },
}
