db.C2H = {}
db.C2H.atomicConstituents = {C=2,H=1,}
db.C2H.charge = 0
db.C2H.M = {
   value = 25.029340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H.gamma = {
   value = 1.2465e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.88965733E+00,
          1.34099611E-02,
         -2.84769501E-05,
          2.94791045E-08,
         -1.09331511E-11,
          6.68393932E+04,
          6.22296438E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.16780652E+00,
          4.75221902E-03,
         -1.83787077E-06,
          3.04190252E-10,
         -1.77232770E-14,
          6.71210650E+04,
          6.63589475E+00,
      }
   }
}
