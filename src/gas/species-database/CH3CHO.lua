db.CH3CHO = {}
db.CH3CHO.atomicConstituents = {C=2,H=4,O=1,}
db.CH3CHO.charge = 0
db.CH3CHO.M = {
   value = 44.052560e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3CHO.gamma = {
   value = 1.1762e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.47294595E+01,
         -0.31932858E-02,
          0.47534921E-04,
         -0.57458611E-07,
          0.21931112E-10,
         -0.21572878E+05,
          0.41030159E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.54041108E+01,
          0.11723059E-01,
         -0.42263137E-05,
          0.68372451E-09,
         -0.40984863E-13,
         -0.22593122E+05,
         -0.34807917E+01,
      }
   }
}
