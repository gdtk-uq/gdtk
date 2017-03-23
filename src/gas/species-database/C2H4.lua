db.C2H4 = {}
db.C2H4.atomicConstituents = {C=2,H=4,}
db.C2H4.charge = 0
db.C2H4.M = {
   value = 28.053160e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H4.gamma = {
   value = 1.2393e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H4.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.95920148E+00,
         -7.57052247E-03,
          5.70990292E-05,
         -6.91588753E-08,
          2.69884373E-11,
          5.08977593E+03,
          4.09733096E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.03611116E+00,
          1.46454151E-02,
         -6.71077915E-06,
          1.47222923E-09,
         -1.25706061E-13,
          4.93988614E+03,
          1.03053693E+01,
      }
   }
}
