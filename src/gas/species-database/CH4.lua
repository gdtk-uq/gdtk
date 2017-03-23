db.CH4 = {}
db.CH4.atomic_constituents = {C=1,H=4,}
db.CH4.charge = 0
db.CH4.M = {
   value = 16.0424600e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH4.gamma = {
   value = 1.303,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH4.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          5.14987613E+00,
         -1.36709788E-02,
          4.91800599E-05,
         -4.84743026E-08,
          1.66693956E-11,
         -1.02466476E+04,
         -4.64130376E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          7.48514950E-02,
          1.33909467E-02,
         -5.73285809E-06,
          1.22292535E-09,
         -1.01815230E-13,
         -9.46834459E+03,
          1.84373180E+01,
      }
   }
}
