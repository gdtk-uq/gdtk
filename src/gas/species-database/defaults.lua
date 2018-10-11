-- Curator: Rowan J. Gollan
-- Date: 08-Mar-2015
--
-- History:
--   08-Mar-2015 -- first cut.
--               -- Experiment to see if this form has longevity.
--   24-Nov-2015 -- Split into separate files.
--               -- Introduced a "defaults" table which is just
--               -- the properties for air.


db = {}

db.default = {}
db.default.atomicConstituents = {}
db.default.charge = 0
db.default.M = {
   value = 28.96e-3,
   units = 'kg/mol',
}
db.default.gamma = {
   value = 1.4,
   note = "valid at low temperatures for diatomic molecules"
}
db.default.sigma = {
   value = 3.621,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'taken from N2.'
}
db.default.epsilon = {
   value = 97.530,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'taken from N2.'
}
db.default.Lewis = {
   value = 1.152,
   reference = 'taken from N2'
}
db.default.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.default.sutherlandVisc = {
   mu_ref = 1.716e-5, 
   T_ref = 273.0,
   S = 111.0,
   reference = "Table 1-2, White (2006)"
}
db.default.sutherlandThermCond = {
   T_ref = 273.0, 
   k_ref = 0.0241, 
   S = 194.0,
   reference = "Table 1-3, White (2006)"
}
db.default.ceaThermoCoeffs = {
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
-- CEA viscosity uses N2 values
db.default.ceaViscosity = {
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
-- CEA thermal conductivity uses N2 values
db.default.ceaThermCond = {
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


