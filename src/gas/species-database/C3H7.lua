db.C3H7 = {}
db.C3H7.atomicConstituents = {C=3,H=7,}
db.C3H7.charge = 0
db.C3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H7.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.10515518E+01,
          0.25991980E-01,
          0.23800540E-05,
         -0.19609569E-07,
          0.93732470E-11,
          0.10631863E+05,
          0.21122559E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.77026987E+01,
          0.16044203E-01,
         -0.52833220E-05,
          0.76298590E-09,
         -0.39392284E-13,
          0.82984336E+04,
         -0.15480180E+02,
      }
   }
}
