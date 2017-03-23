db.NH = {}
db.NH.atomicConstituents = {N=1,H=1,}
db.NH.charge = 0
db.NH.M = {
   value = 15.0146400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NH.gamma = {
   value = 1.398,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.34929085E+01,
          0.31179198E-03,
         -0.14890484E-05,
          0.24816442E-08,
         -0.10356967E-11,
          0.41880629E+05,
          0.18483278E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.27836928E+01,
          0.13298430E-02,
         -0.42478047E-06,
          0.78348501E-10,
         -0.55044470E-14,
          0.42120848E+05,
          0.57407799E+01,
      }
   }
}
