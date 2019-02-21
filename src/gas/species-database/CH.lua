db.CH = {}
db.CH.atomicConstituents = {C=1,H=1,}
db.CH.charge = 0
db.CH.M = {
   value = 13.0186400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH.sigma = {
   value = 2.750,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH.epsilon = {
   value = 80.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         2.220590133e+04,
        -3.405411530e+02,
         5.531452290e+00,
        -5.794964260e-03,
         7.969554880e-06,
        -4.465911590e-09,
         9.596338320e-13,
         7.240783270e+04,
        -9.107673050e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         2.060763440e+06,
        -5.396206660e+03,
         7.856293850e+00,
        -7.965907450e-04,
         1.764308305e-07,
        -1.976386267e-11,
         5.030429510e-16,
         1.062236592e+05,
        -3.154757439e+01
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
        -8.068368690e+08,
         4.575450540e+05,
        -9.843975080e+01,
         1.235244098e-02,
        -8.485608570e-07,
         3.040410624e-11,
        -4.400315170e-16,
        -3.595851590e+06,
         8.953477440e+02
      }
   }
}
db.CH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.48981665E+00,
          3.23835541E-04,
         -1.68899065E-06,
          3.16217327E-09,
         -1.40609067E-12,
          7.07972934E+04,
          2.08401108E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.87846473E+00,
          9.70913681E-04,
          1.44445655E-07,
         -1.30687849E-10,
          1.76079383E-14,
          7.10124364E+04,
          5.48497999E+00,
      }
   }
}
