db.CH2CHO = {}
db.CH2CHO.atomicConstituents = {C=2,H=3,O=1,}
db.CH2CHO.charge = 0
db.CH2CHO.M = {
   value = 43.044620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CHO.gamma = {
   value = 1.1776e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.03409062E+02,
          0.10738574E-01,
          0.01891492E-04,
         -0.07158583E-07,
          0.02867385E-10,
          0.15214766E+04,
          0.09558290E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.05975670E+02,
          0.08130591E-01,
         -0.02743624E-04,
          0.04070304E-08,
         -0.02176017E-12,
          0.04903218E+04,
         -0.05045251E+02,
      }
   }
}
