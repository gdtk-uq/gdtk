db.NH2 = {}
db.NH2.atomicConstituents = {N=1,H=2,}
db.NH2.charge = 0
db.NH2.M = {
   value = 16.022580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH2.gamma = {
   value = 1.3254e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.42040029E+01,
         -0.21061385E-02,
          0.71068348E-05,
         -0.56115197E-08,
          0.16440717E-11,
          0.21885910E+05,
         -0.14184248E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.28347421E+01,
          0.32073082E-02,
         -0.93390804E-06,
          0.13702953E-09,
         -0.79206144E-14,
          0.22171957E+05,
          0.65204163E+01,
      }
   }
}
