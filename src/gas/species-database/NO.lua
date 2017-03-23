db.NO = {}
db.NO.atomicConstituents = {O=1,N=1,}
db.NO.charge = 0
db.NO.M = {
   value = 0.03000610,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NO.gamma = {
   value = 1.38600000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NO.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -1.143916503e+04,
          1.536467592e+02,
          3.431468730e+00,
         -2.668592368e-03,
          8.481399120e-06,
         -7.685111050e-09,
          2.386797655e-12,
          9.098214410e+03,
          6.728725490e+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          2.239018716e+05,
         -1.289651623e+03,
          5.433936030e+00,
         -3.656034900e-04,
          9.880966450e-08,
         -1.416076856e-11,
          9.380184620e-16,
          1.750317656e+04,
         -8.501669090e+00,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
         -9.575303540e+08,
          5.912434480e+05,
         -1.384566826e+02,
          1.694339403e-02,
         -1.007351096e-06,
          2.912584076e-11,
         -3.295109350e-16,
         -4.677501240e+06,
          1.242081216e+03,
      }
   },
}
db.NO.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0262029e-01,
      B = -6.2017783e+01,
      C = -1.3954524e+02,
      D =  2.0268332e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.8009050e-01,
      B =  3.0486891e+02,
      C = -9.4847722e+04,
      D =  5.2873381e-01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.0580582e-01,
      B =  6.2427878e+02,
      C = -5.7879210e+05,
      D =  2.6516450e-01
   },
}
db.NO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.42184763E+01,
         -0.46389760E-02,
          0.11041022E-04,
         -0.93361354E-08,
          0.28035770E-11,
          0.98446230E+04,
          0.22808464E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.32606056E+01,
          0.11911043E-02,
         -0.42917048E-06,
          0.69457669E-10,
         -0.40336099E-14,
          0.99209746E+04,
          0.63693027E+01,
      }
   }
}
db.NO.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  9.5028758e-01,
      B =  7.6667058e+01,
      C = -9.9894764e+03,
      D = -6.2776717e-03
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.6215238e-01,
      B =  4.4568223e+02,
      C = -2.3856466e+05,
      D =  4.6209876e-01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.0377865e+00,
      B = -3.4486864e+04,
      C =  6.7451187e+07,
      D =  2.0928749e+01
   },
}
