db.CH3O = {}
db.CH3O.atomicConstituents = {C=1,H=3,O=1,}
db.CH3O.charge = 0
db.CH3O.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3O.gamma = {
   value = 1.2802e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.02106204E+02,
          0.07216595E-01,
          0.05338472E-04,
         -0.07377636E-07,
          0.02075610E-10,
          0.09786011E+04,
          0.13152177E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3000.0,
      coeffs = {
         0,
         0,
          0.03770799E+02,
          0.07871497E-01,
         -0.02656384E-04,
          0.03944431E-08,
         -0.02112616E-12,
          0.12783252E+03,
          0.02929575E+02,
      }
   }
}
