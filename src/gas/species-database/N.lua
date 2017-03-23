db.N = {}
db.N.atomicConstituents = {N=1}
db.N.charge = 0
db.N.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   reference = 'CEA2::thermo.inp'
}
db.N.gamma =  {
   value = 5./3.,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.N.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = { 0.000000000e+00,
                 0.000000000e+00,
                 2.500000000e+00,
   		 0.000000000e+00,
		 0.000000000e+00,
		 0.000000000e+00,
		 0.000000000e+00,
		 5.610463780e+04,
		 4.193905036e+00 }		 
  },
  segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = { 8.876501380e+04,
       	        -1.071231500e+02,
                 2.362188287e+00,
   		 2.916720081e-04,
	        -1.729515100e-07,
                 4.012657880e-11,
		-2.677227571e-15,
		 5.697351330e+04,
		 4.865231506e+00 }
  },
  segment2 = {
    T_lower  = 6000.0,
    T_upper= 20000.0,
    coeffs = { 5.475181050e+08,
              -3.107574980e+05,
	       6.916782740e+01,
	      -6.847988130e-03,
	       3.827572400e-07,
	      -1.098367709e-11,
	       1.277986024e-16,
	       2.550585618e+06,
	      -5.848769753e+02 }
   },
   reference="from CEA2::thermo.inp"
}
db.N.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.25000000E+01,
          0.00000000E+00,
          0.00000000E+00,
          0.00000000E+00,
          0.00000000E+00,
          0.56104637E+05,
          0.41939087E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.24159429E+01,
          0.17489065E-03,
         -0.11902369E-06,
          0.30226245E-10,
         -0.20360982E-14,
          0.56133773E+05,
          0.46496096E+01,
      }
   }
}
db.N.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83724737e+00,
      B = 0.43997150e+03,
      C = -0.17450753e+06,
      D = 0.10365689e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.89986588e+00,
      B = 0.14112801e+04,
      C = -0.18200478e+07,
      D = -0.55811716e+00
   },
   reference = 'from CEA2::trans.inp which cites Levin et al (1990)'
}
db.N.ceaThermCond = {
   nsegments = 2, 
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83771661e+00,
      B = 0.44243270e+03,
      C = -0.17578446e+06,
      D = 0.89942915e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.90001710e+00,
      B = 0.14141175e+04,
      C = -0.18262403e+07,
      D = 0.24048513e+00
   },
   reference = 'from CEA2::trans.inp which cites Levin et al (1990)'
}

