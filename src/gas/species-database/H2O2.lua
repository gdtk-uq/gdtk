db.H2O2 = {}
db.H2O2.atomicConstituents = {O=2,H=2,}
db.H2O2.charge = 0
db.H2O2.M = {
   value = 0.03401468,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2O2.gamma = {
   value = 1.24400000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2O2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -9.279533580e+04,
          1.564748385e+03,
         -5.976460140e+00,
          3.270744520e-02,
         -3.932193260e-05,
          2.509255235e-08,
         -6.465045290e-12,
         -2.494004728e+04,
          5.877174180e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          1.489428027e+06,
         -5.170821780e+03,
          1.128204970e+01,
         -8.042397790e-05,
         -1.818383769e-08,
          6.947265590e-12,
         -4.827831900e-16,
          1.418251038e+04,
         -4.650855660e+01,
      }
   },
}
db.H2O2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.27611269E+00,
         -5.42822417E-04,
          1.67335701E-05,
         -2.15770813E-08,
          8.62454363E-12,
         -1.77025821E+04,
          3.43505074E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          4.16500285E+00,
          4.90831694E-03,
         -1.90139225E-06,
          3.71185986E-10,
         -2.87908305E-14,
         -1.78617877E+04,
          2.91615662E+00,
      }
   }
}
