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
db.CN.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         3.949148570e+03,
        -1.391590572e+02,
         4.930835320e+00,
        -6.304670510e-03,
         1.256836472e-05,
        -9.878300500e-09,
         2.843137221e-12,
         5.228455380e+04,
        -2.763115585e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
        -2.228006270e+06,
         5.040733390e+03,
        -2.121897722e-01,
         1.354901134e-03,
         1.325929798e-07,
        -6.937006370e-11,
         5.494952270e-15,
         1.784496132e+04,
         3.282563919e+01
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
        -1.794798118e+08,
         1.054346069e+05,
        -1.729624170e+01,
         2.194895530e-03,
        -8.508938030e-08,
         9.318692990e-13,
         6.358139930e-18,
        -7.962594120e+05,
         1.913139639e+02
      }
   }
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
