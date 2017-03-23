db.HCCO = {}
db.HCCO.atomicConstituents = {C=2,H=1,O=1,}
db.HCCO.charge = 0
db.HCCO.M = {
   value = 41.028740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCCO.gamma = {
   value = 1.2067e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.22517214E+01,
          0.17655021E-01,
         -0.23729101E-04,
          0.17275759E-07,
         -0.50664811E-11,
          0.20059449E+05,
          0.12490417E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 4000.0,
      coeffs = {
         0,
         0,
          0.56282058E+01,
          0.40853401E-02,
         -0.15934547E-05,
          0.28626052E-09,
         -0.19407832E-13,
          0.19327215E+05,
         -0.39302595E+01,
      }
   }
}
