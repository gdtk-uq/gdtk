db.Ar = {}
db.Ar.atomicConstituents = {Ar=1,}
db.Ar.charge = 0
db.Ar.M = {
   value = 0.03994800,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.Ar.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Ar.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
          0.000000000e+00,
          0.000000000e+00,
          2.500000000e+00,
          0.000000000e+00,
          0.000000000e+00,
          0.000000000e+00,
          0.000000000e+00,
         -7.453750000e+02,
          4.379674910e+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          2.010538475e+01,
         -5.992661070e-02,
          2.500069401e+00,
         -3.992141160e-08,
          1.205272140e-11,
         -1.819015576e-15,
          1.078576636e-19,
         -7.449939610e+02,
          4.379180110e+00,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
         -9.951265080e+08,
          6.458887260e+05,
         -1.675894697e+02,
          2.319933363e-02,
         -1.721080911e-06,
          6.531938460e-11,
         -9.740147729e-16,
         -5.078300340e+06,
          1.465298484e+03,
      }
   },
}
db.Ar.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.1205763e-01,
      B = -6.7714354e+01,
      C =  1.9040660e+02,
      D =  2.1588272e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.9357334e-01,
      B =  7.0953943e+01,
      C = -2.8386007e+04,
      D =  1.4856447e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.6608935e-01,
      B =  6.7867215e+02,
      C = -8.4991417e+05,
      D =  7.7935167e-01
   },
}
db.Ar.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0968928e-01,
      B = -7.0892249e+01,
      C =  5.8420624e+02,
      D =  1.9337152e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.9075463e-01,
      B =  6.2676058e+01,
      C = -2.5667413e+04,
      D =  1.2664189e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.6269502e-01,
      B =  6.2341752e+02,
      C = -7.1899552e+05,
      D =  5.6927918e-01
   },
}
