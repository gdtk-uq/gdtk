-- Curator: Rowan J. Gollan
-- Date: 08-Mar-2015
--
-- History:
--   08-Mar-2015 -- first cut.
--               -- Experiment to see if this form has longevity.
--

db = {}

db.air = {}
db.air.atomicConstituents = {}
db.air.charge = 0
db.air.M = {
   value = 28.96e-3,
   units = 'kg/mol',
}
db.air.gamma = {
   value = 1.4,
   note = "valid at low temperatures"
}
db.air.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.air.sutherlandVisc = {
   mu_ref = 1.716e-5, 
   T_ref = 273.0,
   S = 111.0,
   reference = "Table 1-2, White (2006)"
}
db.air.sutherlandThermCond = {
   T_ref = 273.0, 
   k_ref = 0.0241, 
   S = 194.0,
   reference = "Table 1-3, White (2006)"
}
db.air.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper =  1000.0,
       coeffs = { 1.009950160e+04,
                 -1.968275610e+02,
                  5.009155110e+00,
	         -5.761013730e-03,
                  1.066859930e-05,
                 -7.940297970e-09,
                  2.185231910e-12,
                 -1.767967310e+02,
                 -3.921504225e+00 },
   },
   segment1 = { 
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {  2.415214430e+05,
                 -1.257874600e+03,
                  5.144558670e+00,
                 -2.138541790e-04,
                  7.065227840e-08,
                 -1.071483490e-11,
                  6.577800150e-16,
                  6.462263190e+03,
                 -8.147411905e+00 }
  }
}


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

db.N2 = {}
db.N2.atomicConstituents = {N=2}
db.N2.charge = 0
db.N2.M = { 
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.N2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.N2.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.N2.sutherlandVisc = {
   mu_ref = 1.663e-05,
   T_ref = 273.0,
   S = 107.0,
   reference = "Table 1-2, White (2006)"
}
db.N2.sutherlandThermCond = { 
   k_ref = 0.0242,
   T_ref = 273.0,
   S = 150.0,
   reference = "Table 1-3, White (2006)"
}
db.N2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = { 2.210371497e+04,
		-3.818461820e+02,
		 6.082738360e+00,
		-8.530914410e-03,
		 1.384646189e-05,
	        -9.625793620e-09,
		 2.519705809e-12,
		 7.108460860e+02,
		-1.076003744e+01
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = { 5.877124060e+05,
                -2.239249073e+03,
                 6.066949220e+00,
		-6.139685500e-04,
                 1.491806679e-07,
                -1.923105485e-11,
                 1.061954386e-15,
                 1.283210415e+04,
                -1.586640027e+01
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = { 8.310139160e+08,
                -6.420733540e+05,
                 2.020264635e+02,
		-3.065092046e-02,
		 2.486903333e-06,
		-9.705954110e-11,
		 1.437538881e-15,
		 4.938707040e+06,
		-1.672099740e+03
      }
   },
   ref="from CEA2::thermo.inp"
}
db.N2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200,
      T_upper = 1000,
      A = 0.62526577,
      B = -31.779652,
      C = -1640.7983,
      D = 1.7454992,
    },
   segment1 = {
      T_lower = 1000,
      T_upper = 5000,
      A = 0.87395209,
      B = 561.52222,
      C = -173948.09,
      D = -0.39335958,
    },
   segment2 = {
      T_lower = 5000,
      T_upper = 15000,
      A = 0.88503551,
      B = 909.02171,
      C = -731290.61,
      D = -0.53503838,
    },
}
db.N2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200,
      T_upper = 1000,
      A = 0.85439436,
      B = 105.73224,
      C = -12347.848,
      D = 0.47793128,
    },
    segment1 =  {
      T_lower = 1000,
      T_upper = 5000,
      A = 0.88407146,
      B = 133.57293,
      C = -11429.64,
      D = 0.24417019,
    },
    segment2 = {
      T_lower = 5000,
      T_upper = 15000,
      A = 2.4176185,
      B = 8047.7749,
      C = 3105580.2,
      D = -14.517761,
    },
}

db.CO2 = {}
db.CO2.atomicConstituents = {}
db.CO2.charge = 0
db.CO2.M = {
   value = 0.04401,
   units = 'kg/mol',
}
db.CO2.gamma = {
   value = 1.3,
   note = "valid at low temperatures"
}
db.CO2.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.CO2.sutherlandVisc = {
   mu_ref = 1.716e-5, 
   T_ref = 273.0,
   S = 111.0,
   reference = "Table 1-2, White (2006)"
}
db.CO2.sutherlandThermCond = {
   T_ref = 273.0, 
   k_ref = 0.0241, 
   S = 194.0,
   reference = "Table 1-3, White (2006)"
}
db.CO2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper =  1000.0,
       coeffs = { 1.009950160e+04,
                 -1.968275610e+02,
                  5.009155110e+00,
	         -5.761013730e-03,
                  1.066859930e-05,
                 -7.940297970e-09,
                  2.185231910e-12,
                 -1.767967310e+02,
                 -3.921504225e+00 },
   },
   segment1 = { 
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {  2.415214430e+05,
                 -1.257874600e+03,
                  5.144558670e+00,
                 -2.138541790e-04,
                  7.065227840e-08,
                 -1.071483490e-11,
                  6.577800150e-16,
                  6.462263190e+03,
                 -8.147411905e+00 }
  }
}
