db.CN = {}
db.CN.atomicConstituents = {C=1,N=1,}
db.CN.charge = 0
db.CN.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CN.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.36129351E+01,
         -0.95551327E-03,
          0.21442977E-05,
         -0.31516323E-09,
         -0.46430356E-12,
          0.51708340E+05,
          0.39804995E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.37459805E+01,
          0.43450775E-04,
          0.29705984E-06,
         -0.68651806E-10,
          0.44134173E-14,
          0.51536188E+05,
          0.27867601E+01,
      }
   }
}
