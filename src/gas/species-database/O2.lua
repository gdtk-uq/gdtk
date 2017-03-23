db.O2 = {}
db.O2.atomicConstituents = {O=2,}
db.O2.charge = 0
db.O2.M = {
   value = 0.03199880,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.O2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.O2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -3.425563420e+04,
          4.847000970e+02,
          1.119010961e+00,
          4.293889240e-03,
         -6.836300520e-07,
         -2.023372700e-09,
          1.039040018e-12,
         -3.391454870e+03,
          1.849699470e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         -1.037939022e+06,
          2.344830282e+03,
          1.819732036e+00,
          1.267847582e-03,
         -2.188067988e-07,
          2.053719572e-11,
         -8.193467050e-16,
         -1.689010929e+04,
          1.738716506e+01,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
          4.975294300e+08,
         -2.866106874e+05,
          6.690352250e+01,
         -6.169959020e-03,
          3.016396027e-07,
         -7.421416600e-12,
          7.278175770e-17,
          2.293554027e+06,
         -5.530621610e+02,
      }
   },
}
db.O2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.78245636E+00,
         -2.99673416E-03,
          9.84730201E-06,
         -9.68129509E-09,
          3.24372837E-12,
         -1.06394356E+03,
          3.65767573E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.28253784E+00,
          1.48308754E-03,
         -7.57966669E-07,
          2.09470555E-10,
         -2.16717794E-14,
         -1.08845772E+03,
          5.45323129E+00,
      }
   }
}
db.O2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0916180e-01,
      B = -5.2244847e+01,
      C = -5.9974009e+02,
      D =  2.0410801e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.2216486e-01,
      B =  1.7550839e+02,
      C = -5.7974816e+04,
      D =  1.0901044e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.3981127e-01,
      B =  3.9194906e+02,
      C = -3.7833168e+05,
      D =  9.0931780e-01
   },
}
db.O2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.7229167e-01,
      B =  6.8463210e+00,
      C = -5.8933377e+03,
      D =  1.2210365e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.0917351e-01,
      B =  2.9124182e+02,
      C = -7.9650171e+04,
      D =  6.4851631e-02
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.1218262e+00,
      B = -1.9286378e+04,
      C =  2.3295011e+07,
      D =  2.0342043e+01
   },
}
