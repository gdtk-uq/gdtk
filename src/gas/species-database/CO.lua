db.CO = {}
db.CO.atomicConstituents = {C=1,O=1,}
db.CO.charge = 0
db.CO.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CO.gamma = {
   value = 1.3992e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.57953347E+00,
         -6.10353680E-04,
          1.01681433E-06,
          9.07005884E-10,
         -9.04424499E-13,
         -1.43440860E+04,
          3.50840928E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.71518561E+00,
          2.06252743E-03,
         -9.98825771E-07,
          2.30053008E-10,
         -2.03647716E-14,
         -1.41518724E+04,
          7.81868772E+00,
      }
   }
}
