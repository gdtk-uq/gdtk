db.N2O = {}
db.N2O.atomicConstituents = {N=2,O=1,}
db.N2O.charge = 0
db.N2O.M = {
   value = 44.012800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.N2O.gamma = {
   value = 1.2735e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.N2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.22571502E+01,
          0.11304728E-01,
         -0.13671319E-04,
          0.96819806E-08,
         -0.29307182E-11,
          0.87417744E+04,
          0.10757992E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.48230729E+01,
          0.26270251E-02,
         -0.95850874E-06,
          0.16000712E-09,
         -0.97752303E-14,
          0.80734048E+04,
         -0.22017207E+01,
      }
   }
}
