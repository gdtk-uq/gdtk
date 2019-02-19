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
db.Ar.sigma = {
   value = 3.330,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.Ar.epsilon = {
   value = 136.500,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.Ar.Lewis = {
   value = 1.173
}
db.Ar.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
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
db.Ar.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.02500000E+02,
          0.00000000E+00,
          0.00000000E+00,
          0.00000000E+00,
          0.00000000E+00,
         -0.07453750E+04,
          0.04366000E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.02500000E+02,
          0.00000000E+00,
          0.00000000E+00,
          0.00000000E+00,
          0.00000000E+00,
         -0.07453750E+04,
          0.04366000E+02,
      }
   }
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
db['Ar+'] = {}
db['Ar+'].atomicConstituents = {Ar=1,}
db['Ar+'].charge = 1
db['Ar+'].M = {
   value = 0.0399474514,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db['Ar+'].gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
--The following values are taken from the Ar.lua file.
db['Ar+'].sigma = {
   value = 3.330,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db['Ar+'].epsilon = {
   value = 136.500,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db['Ar+'].entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
--THe following is from CAE ThermaoCoefficient tables.
db['Ar+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.15,
      T_upper = 1000.0007,
      coeffs = {
         -5.731209170e+04,
          7.930791470e+02,
         -1.717121217e+00,
          1.044184018e-02,
         -1.180207501e-05,
          6.528134780e-09,
         -1.447558130e-12,
          1.790572230e+05,
          2.949150950e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         -3.835965400e+05,
          8.162019700e+02,
          2.301342628e+00,
         -4.952983770e-06,
          1.205108477e-08,
         -2.185050286e-12,
          1.265493898e-16,
          1.771811455e+05,
          7.947507480e+00,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
          1.006884827e+07,
         -6.624361280e+03,
          4.446908200e+00,
         -3.017567664e-04,
          2.612882069e-08,
         -1.201637769e-12,
          2.299206903e-17,
          2.349504137e+05,
         -1.032262257e+01,
      }
   },
   ref="from CEA2::thermo.inp"
}

--No CEA transport data for Ar_plus, just use Ar
db['Ar+'].ceaViscosity = db.Ar.ceaViscosity
db['Ar+'].ceaThermCond = db.Ar.ceaThermCond

db.C = {}
db.C.atomicConstituents = {C=1,}
db.C.charge = 0
db.C.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C.gamma = { 
   value = 1.664,
   units = 'non-dimensional',
    description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.C.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.55423955E+00,
         -3.21537724E-04,
          7.33792245E-07,
         -7.32234889E-10,
          2.66521446E-13,
          8.54438832E+04,
          4.53130848E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.49266888E+00,
          4.79889284E-05,
         -7.24335020E-08,
          3.74291029E-11,
         -4.87277893E-15,
          8.54512953E+04,
          4.80150373E+00,
      }
   }
}
db.C2 = {}
db.C2.atomicConstituents = {C=2,}
db.C2.charge = 0
db.C2.M = { 
   value = 24.0214000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'diatom -- assumed 1.4 at low temperatures'
}
db.C2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = {  5.559634510e+05,
                 -9.980126440e+03,
                  6.681620370e+01,
                 -1.743432724e-01,
                  2.448523051e-04,
                 -1.703467580e-07,
                  4.684527730e-11,
                  1.445869634e+05,
                 -3.448229700e+02
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = { -9.689267930e+05,
                  3.561092990e+03,
                 -5.064138930e-01,
                  2.945154879e-03,
                 -7.139441190e-07,
                  8.670657250e-11,
                 -4.076906810e-15,
                  7.681796830e+04,
                  3.339985240e+01
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = {  6.315145920e+06,
                  1.365420661e+04,
                 -3.996903670e+00,
                  1.937561376e-03,
                 -1.584446580e-07,
                  5.520861660e-12,
                 -7.253735340e-17,
                  9.387024990e+03,
                  6.614329920e+01
      }
   },
   ref="from CEA2::thermo.inp"
}
db.C2.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 300,
      T_upper = 1000,
      A =  0.62126764e+00,
      B = -0.19814414e+02,
      C = -0.16506365e+04,
      D =  0.15582169e+01
    },
   segment1 = {
      T_lower = 1000,
      T_upper = 5000,
      A =  0.64809340e+00,
      B =  0.36749201e+01,
      C = -0.24685282e+04,
      D =  0.13505925e+01
    },
   ref="McBride et al (1983)"
}
db.C2.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 300,
      T_upper = 1000,
      A =  0.11782197e+01,
      B =  0.51596967e+03,
      C = -0.42793543e+05,
      D = -0.20201745e+01
    },
    segment1 =  {
      T_lower = 1000,
      T_upper = 5000,
      A =  0.84536557e+00,
      B =  0.16283010e+03,
      C = -0.21960714e+05,
      D =  0.60979956e+00
    },
    ref="McBride et al (1983)"
}


db.C2H = {}
db.C2H.atomicConstituents = {C=2,H=1,}
db.C2H.charge = 0
db.C2H.M = {
   value = 25.029340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H.gamma = {
   value = 1.2465e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.88965733E+00,
          1.34099611E-02,
         -2.84769501E-05,
          2.94791045E-08,
         -1.09331511E-11,
          6.68393932E+04,
          6.22296438E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.16780652E+00,
          4.75221902E-03,
         -1.83787077E-06,
          3.04190252E-10,
         -1.77232770E-14,
          6.71210650E+04,
          6.63589475E+00,
      }
   }
}
db.C2H2 = {}
db.C2H2.atomicConstituents = {C=2,H=2,}
db.C2H2.charge = 0
db.C2H2.M = {
   value = 26.037280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H2.gamma = {
   value = 1.2321e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H2.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H2.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          8.08681094E-01,
          2.33615629E-02,
         -3.55171815E-05,
          2.80152437E-08,
         -8.50072974E-12,
          2.64289807E+04,
          1.39397051E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          4.14756964E+00,
          5.96166664E-03,
         -2.37294852E-06,
          4.67412171E-10,
         -3.61235213E-14,
          2.59359992E+04,
         -1.23028121E+00,
      }
   }
}
db.C2H3 = {}
db.C2H3.atomicConstituents = {C=2,H=3,}
db.C2H3.charge = 0
db.C2H3.M = {
   value = 27.045220e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H3.gamma = {
   value = 1.2408e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H3.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H3.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.21246645E+00,
          1.51479162E-03,
          2.59209412E-05,
         -3.57657847E-08,
          1.47150873E-11,
          3.48598468E+04,
          8.51054025E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.01672400E+00,
          1.03302292E-02,
         -4.68082349E-06,
          1.01763288E-09,
         -8.62607041E-14,
          3.46128739E+04,
          7.78732378E+00,
      }
   }
}
db.C2H4 = {}
db.C2H4.atomicConstituents = {C=2,H=4,}
db.C2H4.charge = 0
db.C2H4.M = {
   value = 28.053160e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H4.gamma = {
   value = 1.2393e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H4.sigma = {
   value = 3.971,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H4.epsilon = {
   value = 280.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H4.Lewis = {
   value = 1.402
}
db.C2H4.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.95920148E+00,
         -7.57052247E-03,
          5.70990292E-05,
         -6.91588753E-08,
          2.69884373E-11,
          5.08977593E+03,
          4.09733096E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.03611116E+00,
          1.46454151E-02,
         -6.71077915E-06,
          1.47222923E-09,
         -1.25706061E-13,
          4.93988614E+03,
          1.03053693E+01,
      }
   }
}
db.C2H5 = {}
db.C2H5.atomicConstituents = {C=2,H=5,}
db.C2H5.charge = 0
db.C2H5.M = {
   value = 29.061100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H5.gamma = {
   value = 1.1963e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H5.sigma = {
   value = 4.302,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H5.epsilon = {
   value = 252.300,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H5.Lewis = {
   value = 1.551
}
db.C2H5.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.30646568E+00,
         -4.18658892E-03,
          4.97142807E-05,
         -5.99126606E-08,
          2.30509004E-11,
          1.28416265E+04,
          4.70720924E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          1.95465642E+00,
          1.73972722E-02,
         -7.98206668E-06,
          1.75217689E-09,
         -1.49641576E-13,
          1.28575200E+04,
          1.34624343E+01,
      }
   }
}
db.C2H6 = {}
db.C2H6.atomicConstituents = {C=2,H=6,}
db.C2H6.charge = 0
db.C2H6.M = {
   value = 30.069040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H6.gamma = {
   value = 1.1872e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H6.sigma = {
   value = 4.302,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H6.epsilon = {
   value = 252.300,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H6.Lewis = {
   value = 1.546
}
db.C2H6.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.29142492E+00,
         -5.50154270E-03,
          5.99438288E-05,
         -7.08466285E-08,
          2.68685771E-11,
         -1.15222055E+04,
          2.66682316E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          1.07188150E+00,
          2.16852677E-02,
         -1.00256067E-05,
          2.21412001E-09,
         -1.90002890E-13,
         -1.14263932E+04,
          1.51156107E+01,
      }
   }
}
db.C3H6 = {}
db.C3H6.atomicConstituents = {C=3,H=6,}
db.C3H6.charge = 0
db.C3H6.M = {
   value = 42.079740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H6.gamma = {
   value = 1.1474e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H6.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.C3H6.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.C3H6.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.01493307E+02,
          0.02092518E+00,
          0.04486794E-04,
         -0.01668912E-06,
          0.07158146E-10,
          0.01074826E+05,
          0.01614534E+03,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.06732257E+02,
          0.01490834E+00,
         -0.04949899E-04,
          0.07212022E-08,
         -0.03766204E-12,
         -0.09235703E+04,
         -0.01331335E+03,
      }
   }
}
db.C3H7 = {}
db.C3H7.atomicConstituents = {C=3,H=7,}
db.C3H7.charge = 0
db.C3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H7.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.10515518E+01,
          0.25991980E-01,
          0.23800540E-05,
         -0.19609569E-07,
          0.93732470E-11,
          0.10631863E+05,
          0.21122559E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.77026987E+01,
          0.16044203E-01,
         -0.52833220E-05,
          0.76298590E-09,
         -0.39392284E-13,
          0.82984336E+04,
         -0.15480180E+02,
      }
   }
}
db.C3H8 = {}
db.C3H8.atomicConstituents = {C=3,H=8,}
db.C3H8.charge = 0
db.C3H8.M = {
   value = 44.095620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H8.gamma = {
   value = 1.1267e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H8.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.93355381E+00,
          0.26424579E-01,
          0.61059727E-05,
         -0.21977499E-07,
          0.95149253E-11,
         -0.13958520E+05,
          0.19201691E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.75341368E+01,
          0.18872239E-01,
         -0.62718491E-05,
          0.91475649E-09,
         -0.47838069E-13,
         -0.16467516E+05,
         -0.17892349E+02,
      }
   }
}
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
db.CH2 = {}
db.CH2.atomicConstituents = {C=1,H=2,}
db.CH2.charge = 0
db.CH2.M = {
   value = 14.0265800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH2.gamma = {
   value = 1.311,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH2.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2.Lewis = {
   value = 1.023
}
db.CH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.76267867E+00,
          9.68872143E-04,
          2.79489841E-06,
         -3.85091153E-09,
          1.68741719E-12,
          4.60040401E+04,
          1.56253185E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.87410113E+00,
          3.65639292E-03,
         -1.40894597E-06,
          2.60179549E-10,
         -1.87727567E-14,
          4.62636040E+04,
          6.17119324E+00,
      }
   }
}
db.CH2_S = {}
db.CH2_S.atomicConstituents = {C=1,H=2,}
db.CH2_S.charge = 0
db.CH2_S.M = {
   value = 14.026580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2_S.gamma = {
   value = 1.3263e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2_S.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2_S.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2_S.Lewis = {
   value = 1.022
}
db.CH2_S.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.19860411E+00,
         -2.36661419E-03,
          8.23296220E-06,
         -6.68815981E-09,
          1.94314737E-12,
          5.04968163E+04,
         -7.69118967E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.29203842E+00,
          4.65588637E-03,
         -2.01191947E-06,
          4.17906000E-10,
         -3.39716365E-14,
          5.09259997E+04,
          8.62650169E+00,
      }
   }
}
db['CH2*'] = {}
db['CH2*'].atomicConstituents = {C=1,H=2,}
db['CH2*'].charge = 0
db['CH2*'].M = {
   value = 14.026580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db['CH2*'].gamma = {
   value = 1.3263e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db['CH2*'].sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db['CH2*'].epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}

db['CH2*'].grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.19860411E+00,
         -2.36661419E-03,
          8.23296220E-06,
         -6.68815981E-09,
          1.94314737E-12,
          5.04968163E+04,
         -7.69118967E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.29203842E+00,
          4.65588637E-03,
         -2.01191947E-06,
          4.17906000E-10,
         -3.39716365E-14,
          5.09259997E+04,
          8.62650169E+00,
      }
   }
}
db.CH2CHO = {}
db.CH2CHO.atomicConstituents = {C=2,H=3,O=1,}
db.CH2CHO.charge = 0
db.CH2CHO.M = {
   value = 43.044620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CHO.gamma = {
   value = 1.1776e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CHO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CHO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.03409062E+02,
          0.10738574E-01,
          0.01891492E-04,
         -0.07158583E-07,
          0.02867385E-10,
          0.15214766E+04,
          0.09558290E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.05975670E+02,
          0.08130591E-01,
         -0.02743624E-04,
          0.04070304E-08,
         -0.02176017E-12,
          0.04903218E+04,
         -0.05045251E+02,
      }
   }
}
db.CH2CO = {}
db.CH2CO.atomicConstituents = {C=2,H=2,O=1,}
db.CH2CO.charge = 0
db.CH2CO.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CO.gamma = {
   value = 1.1908e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.13583630E+00,
          1.81188721E-02,
         -1.73947474E-05,
          9.34397568E-09,
         -2.01457615E-12,
         -7.04291804E+03,
          1.22156480E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          4.51129732E+00,
          9.00359745E-03,
         -4.16939635E-06,
          9.23345882E-10,
         -7.94838201E-14,
         -7.55105311E+03,
          6.32247205E-01,
      }
   }
}
db.CH2O = {}
db.CH2O.atomicConstituents = {C=1,H=2,O=1}
db.CH2O.charge = 0
db.CH2O.M = {
   value = 30.025980e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2O.gamma = {
   value = 1.3065e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2O.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2O.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2O.Lewis = {
   value = 1.329
}
db.CH2O.ceaThermoCoeffs = {
   notes = 'converted data from Eilmer 3, origin. Chemkin Thermo. Database Kee et al. (1993)',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
     coeffs = {
   	 0.000000000e+00, 
	 0.000000000e+00, 
	 4.793723150e+00, 
	-9.908333690e-03, 
	 3.732200080e-05, 
	-3.792852610e-08, 
	 1.317726520e-11, 
	-1.430895670e+04, 
	 6.028129000e-01, 
      }
  },
  segment1 = {
     T_lower = 1000.0,
     T_upper = 3500.0,
     coeffs = {
   	 0.000000000e+00, 
	 0.000000000e+00, 
	 1.760690080e+00, 
	 9.200000820e-03, 
	-4.422588130e-06, 
	 1.006412120e-09, 
	-8.838556400e-14, 
	-1.399583230e+04, 
	 1.365632300e+01, 
      }
   }
}
db.CH2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.79372315E+00,
         -9.90833369E-03,
          3.73220008E-05,
         -3.79285261E-08,
          1.31772652E-11,
         -1.43089567E+04,
          6.02812900E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          1.76069008E+00,
          9.20000082E-03,
         -4.42258813E-06,
          1.00641212E-09,
         -8.83855640E-14,
         -1.39958323E+04,
          1.36563230E+01,
      }
   }
}
db.CH2OH = {}
db.CH2OH.atomicConstituents = {C=1,H=3,O=1,}
db.CH2OH.charge = 0
db.CH2OH.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2OH.gamma = {
   value = 1.2070e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2OH.sigma = {
   value = 3.690,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2OH.epsilon = {
   value = 417.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.86388918E+00,
          5.59672304E-03,
          5.93271791E-06,
         -1.04532012E-08,
          4.36967278E-12,
         -3.19391367E+03,
          5.47302243E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.69266569E+00,
          8.64576797E-03,
         -3.75101120E-06,
          7.87234636E-10,
         -6.48554201E-14,
         -3.24250627E+03,
          5.81043215E+00,
      }
   }
}
db.CH3 = {}
db.CH3.atomicConstituents = {C=1,H=3,}
db.CH3.charge = 0
db.CH3.M = {
   value = 15.0345200e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH3.gamma = {
   value = 1.276,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH3.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.Lewis = {
   value = 1.049
}
db.CH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.67359040E+00,
          2.01095175E-03,
          5.73021856E-06,
         -6.87117425E-09,
          2.54385734E-12,
          1.64449988E+04,
          1.60456433E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.28571772E+00,
          7.23990037E-03,
         -2.98714348E-06,
          5.95684644E-10,
         -4.67154394E-14,
          1.67755843E+04,
          8.48007179E+00,
      }
   }
}
db.CH3CHO = {}
db.CH3CHO.atomicConstituents = {C=2,H=4,O=1,}
db.CH3CHO.charge = 0
db.CH3CHO.M = {
   value = 44.052560e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3CHO.gamma = {
   value = 1.1762e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3CHO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3CHO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.47294595E+01,
         -0.31932858E-02,
          0.47534921E-04,
         -0.57458611E-07,
          0.21931112E-10,
         -0.21572878E+05,
          0.41030159E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.54041108E+01,
          0.11723059E-01,
         -0.42263137E-05,
          0.68372451E-09,
         -0.40984863E-13,
         -0.22593122E+05,
         -0.34807917E+01,
      }
   }
}
db.CH3O = {}
db.CH3O.atomicConstituents = {C=1,H=3,O=1,}
db.CH3O.charge = 0
db.CH3O.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3O.gamma = {
   value = 1.2802e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3O.sigma = {
   value = 3.690,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3O.epsilon = {
   value = 417.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3O.Lewis = {
   value = 1.360
}
db.CH3O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.02106204E+02,
          0.07216595E-01,
          0.05338472E-04,
         -0.07377636E-07,
          0.02075610E-10,
          0.09786011E+04,
          0.13152177E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3000.0,
      coeffs = {
         0,
         0,
          0.03770799E+02,
          0.07871497E-01,
         -0.02656384E-04,
          0.03944431E-08,
         -0.02112616E-12,
          0.12783252E+03,
          0.02929575E+02,
      }
   }
}
db.CH3OH = {}
db.CH3OH.atomicConstituents = {C=1,H=4,O=1,}
db.CH3OH.charge = 0
db.CH3OH.M = {
   value = 32.041860e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3OH.gamma = {
   value = 1.2320e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          5.71539582E+00,
         -1.52309129E-02,
          6.52441155E-05,
         -7.10806889E-08,
          2.61352698E-11,
         -2.56427656E+04,
         -1.50409823E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          1.78970791E+00,
          1.40938292E-02,
         -6.36500835E-06,
          1.38171085E-09,
         -1.17060220E-13,
         -2.53748747E+04,
          1.45023623E+01,
      }
   }
}
db.CH4 = {}
db.CH4.atomicConstituents = {C=1,H=4,}
db.CH4.charge = 0
db.CH4.M = {
   value = 16.0424600e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH4.gamma = {
   value = 1.303,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH4.sigma = {
   value = 3.746,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH4.epsilon = {
   value = 141.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH4.Lewis = {
   value = 1.043
}
db.CH4.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          5.14987613E+00,
         -1.36709788E-02,
          4.91800599E-05,
         -4.84743026E-08,
          1.66693956E-11,
         -1.02466476E+04,
         -4.64130376E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          7.48514950E-02,
          1.33909467E-02,
         -5.73285809E-06,
          1.22292535E-09,
         -1.01815230E-13,
         -9.46834459E+03,
          1.84373180E+01,
      }
   }
}
db.CH4.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = { 
	-1.766850998e+05,  
	 2.786181020e+03, 
	-1.202577850e+01,
	 3.917619290e-02, 
	-3.619054430e-05,  
	 2.026853043e-08,
	-4.976705490e-12, 
	-2.331314360e+04,  
	 8.904322750e+01
      }
   },
   segment1 = { 
      T_lower  = 1000.0,
      T_upper = 6000.0,
      coeffs = {  
	 3.730042760e+06, 
	-1.383501485e+04,  
	 2.049107091e+01,
	-1.961974759e-03,
	 4.727313040e-07, 
	-3.728814690e-11,
	 1.623737207e-15,
	 7.532066910e+04, 
	 -1.219124889e+02,
      }  
   }
}

db.CH4.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower =200.0,
      T_upper =1000.0,
      A= 0.57643622e+00, 
      B=-0.93704079e+02, 
      C=0.86992395e+03, 
      D=0.17333347e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A= 0.66400044e+00, 
      B=0.10860843e+02, 
      C=-0.76307841e+04, 
      D=0.10323984e+01
      -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
}

db.CH4.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper =1000.0, 
      A = 0.10238177e+01, 
      B=-0.31092375e+03, 
      C=0.32944309e+05, 
      D=0.67787437e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A = 0.77485028e+00, 
      B = -0.40089627e+03, 
      C = -0.46551082e+05, 
      D = 0.25671481e+01
   },
      --ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
}


db.CHO = {}
db.CHO.atomicConstituents = {C=1,H=1,O=1,}
db.CHO.charge = 0
db.CHO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp',
}
db.CHO.gamma = {
   value = 1.316,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = "Gokel (2004), Dean's Handbook of Organic Chemistry",
}
db.CHO.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CHO.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CHO.Lewis = {
   value = 1.314
}
db.CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.22118584E+00,
         -3.24392532E-03,
          1.37799446E-05,
         -1.33144093E-08,
          4.33768865E-12,
          3.83956496E+03,
          3.39437243E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.77217438E+00,
          4.95695526E-03,
         -2.48445613E-06,
          5.89161778E-10,
         -5.33508711E-14,
          4.01191815E+03,
          9.79834492E+00,
      }
   }
}
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
db.CO = {}
db.CO.atomicConstituents = {C=1,O=1,}
db.CO.charge = 0
db.CO.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CO.gamma = {
   value = 1.3992e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CO.sigma = {
   value = 3.650,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO.epsilon = {
   value = 98.100,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO.Lewis = {
   value = 1.171
}
db.CO.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = { 
	 1.489045326e+04,
	-2.922285939e+02,
	 5.724527170e+00,
	-8.176235030e-03,
	 1.456903469e-05,
	-1.087746302e-08,
	 3.027941827e-12,
	-1.303131878e+04,
	-7.859241350e+00
      }
   },
   segment1 = { 
      T_lower  = 1000.0,
      T_upper = 6000.0,
      coeffs = {  
	 4.619197250e+05,
	-1.944704863e+03,
	 5.916714180e+00,
	-5.664282830e-04,
	 1.398814540e-07,
	-1.787680361e-11,
	 9.620935570e-16,
	-2.466261084e+03,
	-1.387413108e+01
      }  
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = { 
	 8.868662960e+08,
	-7.500377840e+05,
	 2.495474979e+02,
	-3.956351100e-02,
	 3.297772080e-06,
	-1.318409933e-10,
	 1.998937948e-15,
	 5.701421130e+06,
	-2.060704786e+03
      }
   },
}



db.CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.57953347E+00,
         -6.10353680E-04,
          1.01681433E-06,
          9.07005884E-10,
         -9.04424499E-13,
         -1.43440860E+04,
          3.50840928E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.71518561E+00,
          2.06252743E-03,
         -9.98825771E-07,
          2.30053008E-10,
         -2.03647716E-14,
         -1.41518724E+04,
          7.81868772E+00,
      }
   }
}


db.CO.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower =200.0,
      T_upper =1000.0,
      A= 0.62526577e+00, 
      B=-0.31779652e+02, 
      C=-0.16407983e+04, 
      D=0.17454992e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A= 0.87395209e+00,
      B= 0.56152222e+03,
      C= -0.17394809e+06,
      D= -0.39335958e+00
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A= 0.88503551e+00,
      B= 0.90902171e+03,
      C= -0.73129061e+06,
      D= -0.53503838e+00
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
}

db.CO.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper =1000.0, 
      A = 0.85439436e+00,
      B = 0.10573224e+03,
      C = -0.12347848e+05,
      D = 0.47793128e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A = 0.88407146e+00,
      B = 0.13357293e+03,
      C = -0.11429640e+05,
      D = 0.24417019E+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A = 0.24175411e+01,
      B = 0.80462671e+04,
      C = 0.31090740e+07,
      D = -0.14516932e+02
   },
      --ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
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
db.CO2.sigma = {
   value = 3.763,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO2.epsilon = {
   value = 244.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO2.Lewis = {
   value = 1.404
}
db.CO2.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.CO2.sutherlandVisc = {
   mu_ref = 14.8e-6, 
   T_ref = 293.15,
   S = 240.0,
   reference = "Crane Company (1988) - Flow of fluids through valves, fittings and pipes"
}
db.CO2.sutherlandThermCond = {
   T_ref = 273.0, --these have not been updated
   k_ref = 0.0241, --these have not been updated
   S = 194.0,--these have not been updated
   reference = "Table 1-3, White (2006)"
}

db.CO2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper =  1000.0,
      coeffs = {
	 4.943650540e+04,
	-6.264116010e+02,
	 5.301725240e+00,
	 2.503813816e-03,
	-2.127308728e-07,
	-7.689988780e-10,
	 2.849677801e-13,
	-4.528198460e+04,
	-7.048279440e+00
      }
   },
   segment1 = { 
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 1.176962419e+05,
	-1.788791477e+03,
	 8.291523190e+00,
	-9.223156780e-05,
	 4.863676880e-09,
	-1.891053312e-12,
	 6.330036590e-16,
	-3.908350590e+04,
	-2.652669281e+01
      }
   },
   segment2 = { 
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	-1.544423287e+09,
	 1.016847056e+06,
	-2.561405230e+02,
	 3.369401080e-02,
	-2.181184337e-06,
	 6.991420840e-11,
	-8.842351500e-16,
	-8.043214510e+06,
	 2.254177493e+03
      }
   } -- from thermo.inp Gurvich, 1991 pt1 p211 pt2 p200
}

db.CO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.35677352E+00,
          8.98459677E-03,
         -7.12356269E-06,
          2.45919022E-09,
         -1.43699548E-13,
         -4.83719697E+04,
          9.90105222E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.85746029E+00,
          4.41437026E-03,
         -2.21481404E-06,
          5.23490188E-10,
         -4.72084164E-14,
         -4.87591660E+04,
          2.27163806E+00,
      }
   }
}

db.CO2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower =200.0,
      T_upper =1000.0,
      A= 0.51137258e+00,
      B= -0.22951321e+03,
      C= 0.13710678e+05,
      D= 0.27075538e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A= 0.63978285e+00,
      B= -0.42637076e+02,
      C= -0.15522605e+05,
      D= 0.16628843e+01
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A= 0.72150912e+00,
      B= 0.75012895e+03,
      C= -0.11825507e+07,
      D= 0.85493645e+00
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   }
}

db.CO2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper =1000.0, 
      A = 0.48056568e+00,
      B = -0.50786720e+03,
      C = 0.35088811e+05,
      D = 0.36747794e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A = 0.69857277e+00,
      B = -0.11830477e+03,
      C = -0.50688859e+05,
      D =  0.18650551e+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A = 0.10518358e+01,
      B = -0.42555944e+04,
      C = 0.14288688e+08,
      D = -0.88950473e+00
   }
      --ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
}

db['e-'] = {}
db['e-'].atomicConstituents = {}
db['e-'].charge = -1
db['e-'].M = {
   value = 0.000548579903e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['e-'].gamma = {
   value = 1.667, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
-- The coefficients are the same over all temperature ranges
-- as given by CEA (which makes sense because there is no
-- internal energy storage in an electron), so we can just
-- use one range.
db['e-'].ceaThermoCoeffs = {
   nsegments = 1,
   segment0 = {
      T_lower = 298.15,
      T_upper = 20000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	-7.453750000e+02,
	-1.172081224e+01
      }
   }
}


db.H = {}
db.H.atomicConstituents = {H=1,}
db.H.charge = 0
db.H.M = {
   value = 0.00100794,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.H.sigma = {
   value = 2.050,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H.epsilon = {
   value = 145.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H.Lewis = {
   value = 0.189
}
db.H.ceaThermoCoeffs = {
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
          2.547370801e+04,
         -4.466828530e-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          6.078774250e+01,
         -1.819354417e-01,
          2.500211817e+00,
         -1.226512864e-07,
          3.732876330e-11,
         -5.687744560e-15,
          3.410210197e-19,
          2.547486398e+04,
         -4.481917770e-01,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
          2.173757694e+08,
         -1.312035403e+05,
          3.399174200e+01,
         -3.813999680e-03,
          2.432854837e-07,
         -7.694275540e-12,
          9.644105630e-17,
          1.067638086e+06,
         -2.742301051e+02,
      }
   },
}
db.H.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.4226149e-01,
      B = -4.0132865e+02,
      C =  1.8554165e+05,
      D =  4.6741844e-02
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7486623e-01,
      B = -2.5022902e+03,
      C =  7.0955048e+06,
      D = -9.3888455e-01
   },
}
db.H.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.4166119e-01,
      B = -4.0487203e+02,
      C =  1.8775642e+05,
      D =  3.4843121e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7447639e-01,
      B = -2.5089452e+03,
      C =  7.1081294e+06,
      D =  2.4970991e+00
   },
}
db.H.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.50000000E+00,
          7.05332819E-13,
         -1.99591964E-15,
          2.30081632E-18,
         -9.27732332E-22,
          2.54736599E+04,
         -4.46682853E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.50000001E+00,
         -2.30842973E-11,
          1.61561948E-14,
         -4.73515235E-18,
          4.98197357E-22,
          2.54736599E+04,
         -4.46682914E-01,
      }
   }
}
db['H+'] = {}
db['H+'].atomicConstituents = {H=1}
db['H+'].charge = 1
db['H+'].M = {
   value = 1.0073914e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
--Use Gamma value for H.
db['H+'].gamma = { 
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

db['H+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 298.15,
      T_upper = 1000.0,
      coeffs = { 
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.840214877e+05,
	-1.140646644e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.840214877e+05,
	-1.140646644e+00
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = { 
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.840214877e+05,
	-1.140646644e+00
      }
   },
   ref="from CEA2::thermo.inp"
}

-- No CEA transport data for H_plus, just use H
db['H+'].ceaViscosity = db.H.ceaViscosity
db['H+'].ceaThermCond = db.H.ceaThermCond


db.H2 = {}
db.H2.atomicConstituents = {H=2,}
db.H2.charge = 0
db.H2.M = {
   value = 0.00201588,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.H2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.H2.sigma = {
   value = 2.920,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2.epsilon = {
   value = 38.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2.Lewis = {
   value = 0.317
}
db.H2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
          4.078323210e+04,
         -8.009186040e+02,
          8.214702010e+00,
         -1.269714457e-02,
          1.753605076e-05,
         -1.202860270e-08,
          3.368093490e-12,
          2.682484665e+03,
         -3.043788844e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          5.608128010e+05,
         -8.371504740e+02,
          2.975364532e+00,
          1.252249124e-03,
         -3.740716190e-07,
          5.936625200e-11,
         -3.606994100e-15,
          5.339824410e+03,
         -2.202774769e+00,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
          4.966884120e+08,
         -3.147547149e+05,
          7.984121880e+01,
         -8.414789210e-03,
          4.753248350e-07,
         -1.371873492e-11,
          1.605461756e-16,
          2.488433516e+06,
         -6.695728110e+02,
      }
   },
}
db.H2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.4553182e-01,
      B =  4.3555109e+01,
      C = -3.2579340e+03,
      D =  1.3556243e-01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.6730605e-01,
      B =  6.7931897e+02,
      C = -2.1025179e+05,
      D = -1.8251697e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  1.0126129e+00,
      B =  1.4973739e+03,
      C = -1.4428484e+06,
      D = -2.3254928e+00
   },
}
db.H2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  1.0059461e+00,
      B =  2.7951262e+02,
      C = -2.9792018e+04,
      D =  1.1996252e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  1.0582450e+00,
      B =  2.4875372e+02,
      C =  1.1736907e+04,
      D =  8.2758695e-01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -2.2364420e-01,
      B = -6.9650442e+03,
      C = -7.7771313e+04,
      D =  1.3189369e+01
   },
}
db.H2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          2.34433112E+00,
          7.98052075E-03,
         -1.94781510E-05,
          2.01572094E-08,
         -7.37611761E-12,
         -9.17935173E+02,
          6.83010238E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.33727920E+00,
         -4.94024731E-05,
          4.99456778E-07,
         -1.79566394E-10,
          2.00255376E-14,
         -9.50158922E+02,
         -3.20502331E+00,
      }
   }
}
db.H2CC = {}
db.H2CC.atomicConstituents = {H=2,C=2,}
db.H2CC.charge = 0
db.H2CC.M = {
   value = 26.0372800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2CC.gamma = {
   value = 1.242,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2CC.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file'
}
db.H2CC.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file'
}
db.H2CC.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.32815483E+01,
          0.69764791E-02,
         -0.23855244E-05,
         -0.12104432E-08,
          0.98189545E-12,
          0.48621794E+05,
          0.59203910E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.42780340E+01,
          0.47562804E-02,
         -0.16301009E-05,
          0.25462806E-09,
         -0.14886379E-13,
          0.48316688E+05,
          0.64023701E+00,
      }
   }
}


db.H2CN = {}
db.H2CN.atomicConstituents = {C=1,H=2,N=1,}
db.H2CN.charge = 0
db.H2CN.M = {
   value = 28.033280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.H2CN.gamma = {
   value = 1.2769e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.H2CN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.28516610E+01,
          0.56952331E-02,
          0.10711400E-05,
         -0.16226120E-08,
         -0.23511081E-12,
          0.28637820E+05,
          0.89927511E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 4000.0,
      coeffs = {
         0,
         0,
          0.52097030E+01,
          0.29692911E-02,
         -0.28555891E-06,
         -0.16355500E-09,
          0.30432589E-13,
          0.27677109E+05,
         -0.44444780E+01,
      }
   }
}
db.H2O = {}
db.H2O.atomicConstituents = {O=1,H=2,}
db.H2O.charge = 0
db.H2O.M = {
   value = 0.01801528,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2O.gamma = {
   value = 1.32900000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2O.sigma = {
   value = 2.605,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O.epsilon = {
   value = 572.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O.Lewis = {
   value = 0.854
}
db.H2O.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -3.947960830e+04,
          5.755731020e+02,
          9.317826530e-01,
          7.222712860e-03,
         -7.342557370e-06,
          4.955043490e-09,
         -1.336933246e-12,
         -3.303974310e+04,
          1.724205775e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          1.034972096e+06,
         -2.412698562e+03,
          4.646110780e+00,
          2.291998307e-03,
         -6.836830480e-07,
          9.426468930e-11,
         -4.822380530e-15,
         -1.384286509e+04,
         -7.978148510e+00,
      }
   },
}
db.H2O.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  5.0019557e-01,
      B = -6.9712796e+02,
      C =  8.8163892e+04,
      D =  3.0836508e+00
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  5.8988538e-01,
      B = -5.3769814e+02,
      C =  5.4263513e+04,
      D =  2.3386375e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4330087e-01,
      B = -9.5668913e+01,
      C = -3.7742283e+05,
      D =  1.8125190e+00
   },
}
db.H2O.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  1.0966389e+00,
      B = -5.5513429e+02,
      C =  1.0623408e+05,
      D = -2.4664550e-01
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  3.9367933e-01,
      B = -2.2524226e+03,
      C =  6.1217458e+05,
      D =  5.8011317e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.1858737e-01,
      B = -1.4096649e+04,
      C =  1.9179190e+07,
      D =  1.4345613e+01
   },
}
db.H2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.19864056E+00,
         -2.03643410E-03,
          6.52040211E-06,
         -5.48797062E-09,
          1.77197817E-12,
         -3.02937267E+04,
         -8.49032208E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.03399249E+00,
          2.17691804E-03,
         -1.64072518E-07,
         -9.70419870E-11,
          1.68200992E-14,
         -3.00042971E+04,
          4.96677010E+00,
      }
   }
}
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
db.H2O2.sigma = {
   value = 3.458,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O2.epsilon = {
   value = 107.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
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
db.HCN = {}
db.HCN.atomicConstituents = {C=1,H=1,N=1,}
db.HCN.charge = 0
db.HCN.M = {
   value = 27.0253400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HCN.gamma = {
   value = 1.301,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HCN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.22589886E+01,
          0.10051170E-01,
         -0.13351763E-04,
          0.10092349E-07,
         -0.30089028E-11,
          0.14712633E+05,
          0.89164419E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.38022392E+01,
          0.31464228E-02,
         -0.10632185E-05,
          0.16619757E-09,
         -0.97997570E-14,
          0.14407292E+05,
          0.15754601E+01,
      }
   }
}
db.HCNN = {}
db.HCNN.atomicConstituents = {C=1,H=1,N=2,}
db.HCNN.charge = 0
db.HCNN.M = {
   value = 41.032040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCNN.gamma = {
   value = 1.2032e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCNN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.25243194E+01,
          0.15960619E-01,
         -0.18816354E-04,
          0.12125540E-07,
         -0.32357378E-11,
          0.54261984E+05,
          0.11675870E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.58946362E+01,
          0.39895959E-02,
         -0.15982380E-05,
          0.29249395E-09,
         -0.20094686E-13,
          0.53452941E+05,
         -0.51030502E+01,
      }
   }
}
db.HCNO = {}
db.HCNO.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HCNO.charge = 0
db.HCNO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCNO.gamma = {
   value = 1.2154e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCNO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1382.0,
      coeffs = {
         0,
         0,
          2.64727989E+00,
          1.27505342E-02,
         -1.04794236E-05,
          4.41432836E-09,
         -7.57521466E-13,
          1.92990252E+04,
          1.07332972E+01,
      }
   },
   segment1 = {
      T_lower = 1382.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          6.59860456E+00,
          3.02778626E-03,
         -1.07704346E-06,
          1.71666528E-10,
         -1.01439391E-14,
          1.79661339E+04,
         -1.03306599E+01,
      }
   }
}
db.HCO = {}
db.HCO.atomicConstituents = {C=1,H=1,O=1,}
db.HCO.charge = 0
db.HCO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp',
}
db.HCO.gamma = {
   value = 1.316,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = "Gokel (2004), Dean's Handbook of Organic Chemistry",
}
db.HCO.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCO.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.22118584E+00,
         -3.24392532E-03,
          1.37799446E-05,
         -1.33144093E-08,
          4.33768865E-12,
          3.83956496E+03,
          3.39437243E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.77217438E+00,
          4.95695526E-03,
         -2.48445613E-06,
          5.89161778E-10,
         -5.33508711E-14,
          4.01191815E+03,
          9.79834492E+00,
      }
   }
}
db.HCCO = {}
db.HCCO.atomicConstituents = {C=2,H=1,O=1,}
db.HCCO.charge = 0
db.HCCO.M = {
   value = 41.028740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCCO.gamma = {
   value = 1.2067e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCCO.sigma = {
   value = 2.500,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCO.epsilon = {
   value = 150.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.22517214E+01,
          0.17655021E-01,
         -0.23729101E-04,
          0.17275759E-07,
         -0.50664811E-11,
          0.20059449E+05,
          0.12490417E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 4000.0,
      coeffs = {
         0,
         0,
          0.56282058E+01,
          0.40853401E-02,
         -0.15934547E-05,
          0.28626052E-09,
         -0.19407832E-13,
          0.19327215E+05,
         -0.39302595E+01,
      }
   }
}
db.HCCOH = {}
db.HCCOH.atomicConstituents = {C=2,H=2,O=1,}
db.HCCOH.charge = 0
db.HCCOH.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCCOH.gamma = {
   value = 1.1656e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCCOH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.12423733E+01,
          0.31072201E-01,
         -0.50866864E-04,
          0.43137131E-07,
         -0.14014594E-10,
          0.80316143E+04,
          0.13874319E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.59238291E+01,
          0.67923600E-02,
         -0.25658564E-05,
          0.44987841E-09,
         -0.29940101E-13,
          0.72646260E+04,
         -0.76017742E+01,
      }
   }
}
db.HI = {}
db.HI.atomicConstituents = {I=1,H=1,}
db.HI.charge = 0
db.HI.M = {
   value = 0.12791241,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.HI.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.HI.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
          1.872881730e+04,
         -3.431788840e+02,
          5.956712430e+00,
         -8.543439600e-03,
          1.454780274e-05,
         -1.049104164e-08,
          2.839734003e-12,
          3.682950720e+03,
         -8.149756090e+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          4.724921450e+05,
         -1.923465741e+03,
          5.758048970e+00,
         -4.066266380e-04,
          9.474332050e-08,
         -1.033534431e-11,
          4.611614790e-16,
          1.394857037e+04,
         -1.182487652e+01,
      }
   },
}
db.HNCO = {}
db.HNCO.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HNCO.charge = 0
db.HNCO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HNCO.gamma = {
   value = 1.2173e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HNCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1478.0,
      coeffs = {
         0,
         0,
          3.63096317E+00,
          7.30282357E-03,
         -2.28050003E-06,
         -6.61271298E-10,
          3.62235752E-13,
         -1.55873636E+04,
          6.19457727E+00,
      }
   },
   segment1 = {
      T_lower = 1478.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          6.22395134E+00,
          3.17864004E-03,
         -1.09378755E-06,
          1.70735163E-10,
         -9.95021955E-15,
         -1.66599344E+04,
         -8.38224741E+00,
      }
   }
}
db.HNO = {}
db.HNO.atomicConstituents = {H=1,N=1,O=1}
db.HNO.M = {
   value = 31.01404e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO.charge = 0
db.HNO.gamma = {
   value = 1.325,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
	 -6.854764860e+04,
	  9.551627200e+02,
	 -6.000720210e-01,
	  7.995176750e-03,
	 -6.547079160e-07,
	 -3.670513400e-09,
	  1.783392519e-12,
	  6.435351260e+03,
	  3.048166179e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 -5.795614980e+06,
	  1.945457427e+04,
	 -2.152568374e+01,
	  1.797428992e-02,
	 -4.976040670e-06,
	  6.397924170e-10,
	 -3.142619368e-14,
	 -1.104192372e+05,
	  1.818650338e+02,
      }
   }
}
db.HNO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.45334916E+01,
         -0.56696171E-02,
          0.18473207E-04,
         -0.17137094E-07,
          0.55454573E-11,
          0.11548297E+05,
          0.17498417E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.29792509E+01,
          0.34944059E-02,
         -0.78549778E-06,
          0.57479594E-10,
         -0.19335916E-15,
          0.11750582E+05,
          0.86063728E+01,
      }
   }
}
db.HNO2 = {}
db.HNO2.atomicConstituents = {H=1,N=1,O=2}
db.HNO2.charge = 0
db.HNO2.M = {
   value = 0.04701344,
   unit = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO2.gamma = {
   value = 1.218,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
	  8.591985060e+03,
	  1.203644046e+02,
	  9.412979120e-01,
	  1.942891839e-02,
	 -2.253174194e-05,
	  1.384587594e-08,
	 -3.473550460e-12,
	 -1.106337202e+04,
	  2.073967331e+01,
      }
   },
   segment1 = {
      T_lower  = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	  8.787904130e+05,
	 -3.990455030e+03,
	  1.187349269e+01,
	 -4.881900610e-04,
	  7.133636790e-08,
	 -5.376303340e-12,
	  1.581778986e-16,
	  1.246343241e+04,
	 -4.608874688e+01,
      }
   }
}
db.HNO3 = {}
db.HNO3.atomicConstituents = {H=1,N=1,O=3}
db.HNO3.M = {
   value = 63.01284e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO3.charge = 0
db.HNO3.gamma = {
   value = 1.181,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO3.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
	  9.202869010e+03,
	  1.093774496e+02,
	 -4.521042450e-01,
	  2.984914503e-02,
	 -3.190635500e-05,
	  1.720931528e-08,
	 -3.782649830e-12,
	 -1.764048507e+04,
	  2.746644879e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 -9.497809640e+04,
	 -2.733024468e+03,
	  1.449426995e+01,
	 -7.821868050e-04,
	  1.702693665e-07,
	 -1.930543961e-11,
	  8.870455120e-16,
	 -4.882517780e+03,
	 -5.928392985e+01,
      }
   }
}
db.HO2 = {}
db.HO2.atomicConstituents = {O=2,H=1,}
db.HO2.charge = 0
db.HO2.M = {
   value = 0.03300674,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HO2.gamma = {
   value = 1.31200000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HO2.sigma = {
   value = 3.458,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HO2.epsilon = {
   value = 107.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HO2.Lewis = {
   value = 1.079
}
db.HO2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -7.598882540e+04,
          1.329383918e+03,
         -4.677388240e+00,
          2.508308202e-02,
         -3.006551588e-05,
          1.895600056e-08,
         -4.828567390e-12,
         -5.873350960e+03,
          5.193602140e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         -1.810669724e+06,
          4.963192030e+03,
         -1.039498992e+00,
          4.560148530e-03,
         -1.061859447e-06,
          1.144567878e-10,
         -4.763064160e-15,
         -3.200817190e+04,
          4.066850920e+01,
      }
   },
}
db.HO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.30179801E+00,
         -4.74912051E-03,
          2.11582891E-05,
         -2.42763894E-08,
          9.29225124E-12,
          2.94808040E+02,
          3.71666245E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          4.01721090E+00,
          2.23982013E-03,
         -6.33658150E-07,
          1.14246370E-10,
         -1.07908535E-14,
          1.11856713E+02,
          3.78510215E+00,
      }
   }
}
db.HOCN = {}
db.HOCN.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HOCN.charge = 0
db.HOCN.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HOCN.gamma = {
   value = 1.2185e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HOCN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1368.0,
      coeffs = {
         0,
         0,
          3.78604952E+00,
          6.88667922E-03,
         -3.21487864E-06,
          5.17195767E-10,
          1.19360788E-14,
         -2.82698400E+03,
          5.63292162E+00,
      }
   },
   segment1 = {
      T_lower = 1368.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          5.89784885E+00,
          3.16789393E-03,
         -1.11801064E-06,
          1.77243144E-10,
         -1.04339177E-14,
         -3.70653331E+03,
         -6.18167825E+00,
      }
   }
}
-- Helium ported from Rowan's He.lua file in the cfcfd3 collection
-- PJ, 2017-05-24
db.He = {}
db.He.atomicConstituents = {He=1,}
db.He.charge = 0
db.He.M = {
   value = 4.0026020e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.He.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.He.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 6000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00, 
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00, 
	 0.000000000e+00,
	-7.453750000e+02, 
	 9.287239740e-01
      }
   },
   segment1 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
         3.396845420e+06,
	-2.194037652e+03,
	 3.080231878e+00,
	-8.068957550e-05,
	 6.252784910e-09,
	-2.574990067e-13,
	 4.429960218e-18,
	 1.650518960e+04,
	-4.048814390e+00
      }
   },
   reference = 'cea2::thermo.inp'
   -- NOTE: Data for first 2 ranges has been combined and is same as Chemkin data
}
db.He.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A = 0.75015944e+00,
      B = 0.35763243e+02,
      C =-0.22121291e+04,
      D = 0.92126352e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83394166e+00,
      B = 0.22082656e+03,
      C =-0.52852591e+05,
      D = 0.20809361e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.86316349e+00,
      B = 0.96205176e+03,
      C =-0.12498705e+07,
      D =-0.14115714e+00
   },
   reference = "from CEA2::trans.inp which cites Bich et al. (1990)"
}
db.He.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A = 0.75007833e+00,
      B = 0.36577987e+02,
      C =-0.23636600e+04,
      D =0.29766475e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83319259e+00,
      B = 0.22157417e+03,
      C =-0.53304530e+05,
      D = 0.22684592e+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.85920953e+00,
      B = 0.89873206e+03,
      C =-0.11069262e+07,
      D = 0.19535742e+01
   },
   reference = "from CEA2::trans.inp which cites Bich et al. (1990)"
}
db.I2 = {}
db.I2.atomicConstituents = {I=2,}
db.I2.charge = 0
db.I2.M = {
   value = 0.25380894,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.I2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.I2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -5.087968770e+03,
         -1.249585210e+01,
          4.504219090e+00,
          1.370962533e-04,
         -1.390523014e-07,
          1.174813853e-10,
         -2.337541043e-14,
          6.213469810e+03,
          5.583836940e+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         -5.632594160e+06,
          1.793961560e+04,
         -1.723055169e+01,
          1.244214080e-02,
         -3.332768580e-06,
          4.125477940e-10,
         -1.960461713e-14,
         -1.068505292e+05,
          1.600531883e+02,
      }
   },
}
db.Kr = {}
db.Kr.atomicConstituents = {Kr=1,}
db.Kr.charge = 0
db.Kr.M = {
   value = 0.0838,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.Kr.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Kr.ceaThermoCoeffs = {
   nsegments = 2,
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
	 5.490956510e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 2.643639057e+02,
	-7.910050820e-01,
	 2.500920585e+00,
	-5.328164110e-07,
	 1.620730161e-10,
	-2.467898017e-14,
	 1.478585040e-18,
	-7.403488940e+02,
	 5.484398150e+00
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	-1.375531087e+09,
	 9.064030530e+05,
	-2.403481435e+02,
	 3.378312030e-02,
	-2.563103877e-06,
	 9.969787790e-11,
	-1.521249677e-15,
	-7.111667370e+06,
	 2.086866326e+03
      }
   },
}
db.Kr.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58597795e+00,
      B = -0.12924832e+03,
      C =  0.47495759e+04,
      D =  0.25793650e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.68985845e+00,
      B =  0.56296306e+02,
      C = -0.36082600e+05,
      D =  0.17170715e+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  0.76582939e+00,
      B =  0.68610377e+03,
      C = -0.82537190e+06,
      D =  0.97565128e+00
   },
}
db.Kr.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58008139e+00,
      B = -0.13792556e+03,
      C =  0.60771460e+04,
      D =  0.16420039e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.68859431e+00,
      B =  0.51765647e+02,
      C = -0.34512131e+05,
      D =  0.74332130e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  0.76365443e+00,
      B =  0.65175847e+03,
      C = -0.73589800e+06,
      D =  0.12112126e-01
   },
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
db.N.sigma = {
   value = 3.298,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Svehla (1962)'
}
db.N.epsilon = {
   value = 71.4, 
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Svehla (1962)'
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

db['N+'] = {}
db['N+'].atomicConstituents = {N=1}
db['N+'].charge = 1
db['N+'].M = {
   value = 14.0061514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['N+'].gamma =  {
   value = 1.641, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
db['N+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 298.15,
      T_upper = 1000.0,
      coeffs = {
	 5.237079210e+03,
	 2.299958315e+00,
	 2.487488821e+00,
	 2.737490756e-05,
	-3.134447576e-08,
	 1.850111332e-11,
	-4.447350984e-15,
	 2.256284738e+05,
	 5.076830786e+00
      }
  },
  segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 2.904970374e+05,
	-8.557908610e+02,
	 3.477389290e+00,
	-5.288267190e-04,
	 1.352350307e-07,
	-1.389834122e-11,
	 5.046166279e-16,
	 2.310809984e+05,
	-1.994146545e+00
      }
  },
  segment2 = {
      T_lower  = 6000.0,
      T_upper= 20000.0,
      coeffs = { 
	 1.646092148e+07,
	-1.113165218e+04,
	 4.976986640e+00,
	-2.005393583e-04,
	 1.022481356e-08,
	-2.691430863e-13,
	 3.539931593e-18,
	 3.136284696e+05,
	-1.706646380e+01
      }
   },
   reference="from CEA2::thermo.inp"
}
-- No CEA transport data for N+, just use N
db['N+'].ceaViscosity = db.N.ceaViscosity 
db['N+'].ceaThermCond = db.N.ceaThermCond 


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
db.N2.sigma = {
   value = 3.621,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2.epsilon = {
   value = 97.530,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2.Lewis = {
   value = 1.152
}
db.N2.entropyRefValues = {
   s1 = 6839.91,
   T1 = 298.15,
   p1 = 101.325e3,
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
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
db.N2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.03298677E+02,
          0.14082404E-02,
         -0.03963222E-04,
          0.05641515E-07,
         -0.02444854E-10,
         -0.10208999E+04,
          0.03950372E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      coeffs = {
         0,
         0,
          0.02926640E+02,
          0.14879768E-02,
         -0.05684760E-05,
          0.10097038E-09,
         -0.06753351E-13,
         -0.09227977E+04,
          0.05980528E+02,
      }
   }
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


db['N2+'] = {}
db['N2+'].atomicConstituents = {N=2}
db['N2+'].charge = 1
db['N2+'].M = {
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db['N2+'].gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

db['N2+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 298.15,
      T_upper = 1000.0,
      coeffs = { 
	-3.474047470e+04,
	 2.696222703e+02,
	 3.164916370e+00,
	-2.132239781e-03,
	 6.730476400e-06,
	-5.637304970e-09,
	 1.621756000e-12,
	 1.790004424e+05,
	 6.832974166e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	-2.845599002e+06,
	 7.058893030e+03,
	-2.884886385e+00,
	 3.068677059e-03,
	-4.361652310e-07,
	 2.102514545e-11, 
	 5.411996470e-16,
	 1.340388483e+05,
	 5.090897022e+01
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = { 
	-3.712829770e+08,
	 3.139287234e+05,
	-9.603518050e+01, 
	 1.571193286e-02,
	-1.175065525e-06,
	 4.144441230e-11,
	-5.621893090e-16,
	-2.217361867e+06,
	 8.436270947e+02
      }
   },
   ref="from CEA2::thermo.inp"
}
-- No CEA transport data for N2+, just use N2
db['N2+'].ceaViscosity = db.N2.ceaViscosity
db['N2+'].ceaThermCond = db.N2.ceaThermCond



db.N2O = {}
db.N2O.atomicConstituents = {N=2,O=1,}
db.N2O.charge = 0
db.N2O.M = {
   value = 44.012800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.N2O.gamma = {
   value = 1.2735e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.N2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.22571502E+01,
          0.11304728E-01,
         -0.13671319E-04,
          0.96819806E-08,
         -0.29307182E-11,
          0.87417744E+04,
          0.10757992E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.48230729E+01,
          0.26270251E-02,
         -0.95850874E-06,
          0.16000712E-09,
         -0.97752303E-14,
          0.80734048E+04,
         -0.22017207E+01,
      }
   }
}
db.NCO = {}
db.NCO.atomicConstituents = {C=1,N=1,O=1,}
db.NCO.charge = 0
db.NCO.M = {
   value = 42.016800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NCO.gamma = {
   value = 1.2609e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.28269308E+01,
          0.88051688E-02,
         -0.83866134E-05,
          0.48016964E-08,
         -0.13313595E-11,
          0.14682477E+05,
          0.95504646E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.51521845E+01,
          0.23051761E-02,
         -0.88033153E-06,
          0.14789098E-09,
         -0.90977996E-14,
          0.14004123E+05,
         -0.25442660E+01,
      }
   }
}
db.NH = {}
db.NH.atomicConstituents = {N=1,H=1,}
db.NH.charge = 0
db.NH.M = {
   value = 15.0146400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NH.gamma = {
   value = 1.398,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.34929085E+01,
          0.31179198E-03,
         -0.14890484E-05,
          0.24816442E-08,
         -0.10356967E-11,
          0.41880629E+05,
          0.18483278E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.27836928E+01,
          0.13298430E-02,
         -0.42478047E-06,
          0.78348501E-10,
         -0.55044470E-14,
          0.42120848E+05,
          0.57407799E+01,
      }
   }
}
db.NH2 = {}
db.NH2.atomicConstituents = {N=1,H=2,}
db.NH2.charge = 0
db.NH2.M = {
   value = 16.022580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH2.gamma = {
   value = 1.3254e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.42040029E+01,
         -0.21061385E-02,
          0.71068348E-05,
         -0.56115197E-08,
          0.16440717E-11,
          0.21885910E+05,
         -0.14184248E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.28347421E+01,
          0.32073082E-02,
         -0.93390804E-06,
          0.13702953E-09,
         -0.79206144E-14,
          0.22171957E+05,
          0.65204163E+01,
      }
   }
}
db.NH3 = {}
db.NH3.atomicConstituents = {N=1,H=3,}
db.NH3.charge = 0
db.NH3.M = {
   value = 17.030520e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH3.gamma = {
   value = 1.3036e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.42860274E+01,
         -0.46605230E-02,
          0.21718513E-04,
         -0.22808887E-07,
          0.82638046E-11,
         -0.67417285E+04,
         -0.62537277E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.26344521E+01,
          0.56662560E-02,
         -0.17278676E-05,
          0.23867161E-09,
         -0.12578786E-13,
         -0.65446958E+04,
          0.65662928E+01,
      }
   }
}
db.NNH = {}
db.NNH.atomicConstituents = {N=2,H=1,}
db.NNH.charge = 0
db.NNH.M = {
   value = 29.021340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NNH.gamma = {
   value = 1.3152e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NNH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.43446927E+01,
         -0.48497072E-02,
          0.20059459E-04,
         -0.21726464E-07,
          0.79469539E-11,
          0.28791973E+05,
          0.29779410E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.37667544E+01,
          0.28915082E-02,
         -0.10416620E-05,
          0.16842594E-09,
         -0.10091896E-13,
          0.28650697E+05,
          0.44705067E+01,
      }
   }
}
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
db['NO+'] = {}
db['NO+'].atomicConstituents = {O=1,N=1,}
db['NO+'].charge = 1
db['NO+'].M = {
   value = 30.0055514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['NO+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.15,
      T_upper = 1000.0,
      coeffs = {
	 1.398106635e+03,
	-1.590446941e+02,
	 5.122895400e+00,
	-6.394388620e-03,
	 1.123918342e-05,
	-7.988581260e-09,
	 2.107383677e-12,
	 1.187495132e+05,
	-4.398433810e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 6.069876900e+05,
	-2.278395427e+03,
	 6.080324670e+00,
	-6.066847580e-04,
	 1.432002611e-07,
	-1.747990522e-11,
	 8.935014060e-16,
	 1.322709615e+05,
	-1.519880037e+01
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	 2.676400347e+09,
	-1.832948690e+06,
	 5.099249390e+02, 
	-7.113819280e-02,
	 5.317659880e-06,
	-1.963208212e-10,
	 2.805268230e-15,
	 1.443308939e+07,
	-4.324044462e+03 
      }
   },
}
-- No CEA transport data for NO+, just use NO
db['NO+'].ceaViscosity = db.NO.ceaViscosity  
db['NO+'].ceaThermCond = db.NO.ceaThermCond  
db.NO2 = {}
db.NO2.atomicConstituents = {O=2,N=1}
db.NO2.charge = 0
db.NO2.M = {
   value = 0.0460055000,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NO2.gamma = {
   value = 1.287,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NO2.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
	 -5.642038780e+04,
	  9.633085720e+02,
	 -2.434510974e+00,
	  1.927760886e-02,
	 -1.874559328e-05,
	  9.145497730e-09,
	 -1.777647635e-12,
	 -1.547925037e+03,
	  4.067851210e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	  7.213001570e+05,
	 -3.832615200e+03,
	  1.113963285e+01,
         -2.238062246e-03,
	  6.547723430e-07,
	 -7.611335900e-11,
	  3.328361050e-15,
	  2.502497403e+04,
	 -4.305130040e+01,
      }
   }
}
db.NO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.39440312E+01,
         -0.15854290E-02,
          0.16657812E-04,
         -0.20475426E-07,
          0.78350564E-11,
          0.28966179E+04,
          0.63119917E+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
         0,
         0,
          0.48847542E+01,
          0.21723956E-02,
         -0.82806906E-06,
          0.15747510E-09,
         -0.10510895E-13,
          0.23164983E+04,
         -0.11741695E+00,
      }
   }
}
db.NO2.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A =  0.57379100e+00,
      B = -0.12636034e+03,
      C =  0.21566823e+04,
      D =  0.22287492e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.64239645e+00,
      B =  0.60012144e+00,
      C = -0.27020876e+05,
      D = 0.16570566e+01
   }
}
db.NO2.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A =  0.48574998e+00,
      B = -0.50702110e+03,
      C =  0.46605820e+05,
      D =  0.36444556e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.97660465e+00,
      B =  0.72760751e+03,
      C = -0.32527989e+06,
      D = -0.60899123e+00
   }
}
db.O = {}
db.O.atomicConstituents = {O=1,}
db.O.charge = 0
db.O.M = {
   value = 0.01599940,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.O.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.O.sigma = {
   value = 2.750,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O.epsilon = {
   value = 80.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O.Lewis = {
   value = 0.712
}
db.O.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -7.953611300e+03,
          1.607177787e+02,
          1.966226438e+00,
          1.013670310e-03,
         -1.110415423e-06,
          6.517507500e-10,
         -1.584779251e-13,
          2.840362437e+04,
          8.404241820e+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          2.619020262e+05,
         -7.298722030e+02,
          3.317177270e+00,
         -4.281334360e-04,
          1.036104594e-07,
         -9.438304330e-12,
          2.725038297e-16,
          3.392428060e+04,
         -6.679585350e-01,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
          1.779004264e+08,
         -1.082328257e+05,
          2.810778365e+01,
         -2.975232262e-03,
          1.854997534e-07,
         -5.796231540e-12,
          7.191720164e-17,
          8.890942630e+05,
         -2.181728151e+02,
      }
   },
}
db.O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.16826710E+00,
         -3.27931884E-03,
          6.64306396E-06,
         -6.12806624E-09,
          2.11265971E-12,
          2.91222592E+04,
          2.05193346E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.56942078E+00,
         -8.59741137E-05,
          4.19484589E-08,
         -1.00177799E-11,
          1.22833691E-15,
          2.92175791E+04,
          4.78433864E+00,
      }
   }
}
db.O.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7269241e-01,
      B =  8.3842977e+01,
      C = -5.8502098e+04,
      D =  8.5100827e-01
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7669586e-01,
      B =  1.0158420e+03,
      C = -1.0884566e+06,
      D = -1.8001077e-01
   },
}
db.O.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7271664e-01,
      B =  8.3989100e+01,
      C = -5.8580966e+04,
      D =  1.5179900e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7676666e-01,
      B =  1.0170744e+03,
      C = -1.0906690e+06,
      D =  4.8644232e-01
   },
}
db['O+'] = {}
db['O+'].atomicConstituents = {O=1,}
db['O+'].charge = 1
db['O+'].M = {
   value = 15.9988514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['O+'].gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['O+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.15,
      T_upper = 1000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.879352842e+05,
	 4.393376760e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	-2.166513208e+05,
	 6.665456150e+02,
	 1.702064364e+00,  
	 4.714992810e-04,
	-1.427131823e-07,
	 2.016595903e-11,
	-9.107157762e-16,
	 1.837191966e+05,
	 1.005690382e+01
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	-2.143835383e+08,
	 1.469518523e+05,
	-3.680864540e+01,
	 5.036164540e-03,
	-3.087873854e-07,
	 9.186834870e-12,
	-1.074163268e-16,
	-9.614208960e+05,
	 3.426193080e+02
      }
   },
}
-- No CEA transport data, just use O
db['O+'].ceaViscosity = db.O.ceaViscosity
db['O+'].ceaThermCond = db.O.ceaThermCond 

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
db.O2.sigma = {
   value = 3.458,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O2.epsilon = {
   value = 107.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O2.Lewis = {
   value = 1.086
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
db['O2+'] = {}
db['O2+'].atomicConstituents = {O=2,}
db['O2+'].charge = 1
db['O2+'].M = {
   value = 31.9982514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['O2+'].gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db['O2+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.15,
      T_upper = 1000.0,
      coeffs = {
	-8.607205450e+04,
	 1.051875934e+03,
	-5.432380470e-01,
	 6.571166540e-03,
	-3.274263750e-06,
	 5.940645340e-11,
	 3.238784790e-13,
	 1.345544668e+05,
	 2.902709750e+01
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 7.384654880e+04,
	-8.459559540e+02,
	 4.985164160e+00,
	-1.611010890e-04,
	 6.427083990e-08,
	-1.504939874e-11,
	 1.578465409e-15,
	 1.446321044e+05,
	-5.811230650e+00
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	-1.562125524e+09,
	 1.161406778e+06,
	-3.302504720e+02,
	 4.710937520e-02,
	-3.354461380e-06,
	 1.167968599e-10,
	-1.589754791e-15,
	-8.857866270e+06,
	 2.852035602e+03
      }
   },
}
-- No transport properties in CEA, just set to O2
db['O2+'].ceaViscosity = db.O2.ceaViscosity
db['O2+'].ceaThermCond = db.O2.ceaThermCond

db.O3 = {}
db.O3.atomicConstituents = {O=3}
db.O3.M = {
   value = 47.9982e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.O3.charge = 0
db.O3.gamma = {
   value = 1.267,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.O3.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
	 -1.282314507e+04,
	  5.898216640e+02,
	 -2.547496763e+00,
	  2.690121526e-02,
	 -3.528258340e-05,
	  2.312290922e-08,
	 -6.044893270e-12,
	  1.348368701e+04,
	  3.852218580e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 -3.869662480e+07,
	  1.023344994e+05,
	 -8.961551600e+01,
	  3.706144970e-02,
	 -4.137638740e-06,
	 -2.725018591e-10,
	  5.248188110e-14,
	 -6.517918180e+05,
	  7.029109520e+02,
      }
   }
}
db.OH = {}
db.OH.atomicConstituents = {O=1,H=1,}
db.OH.charge = 0
db.OH.M = {
   value = 0.01700734,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.OH.gamma = {
   value = 1.38600000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.OH.sigma = {
   value = 2.750,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.OH.epsilon = {
   value = 80.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.OH.Lewis = {
   value = 0.736
}
db.OH.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         -1.998858990e+03,
          9.300136160e+01,
          3.050854229e+00,
          1.529529288e-03,
         -3.157890998e-06,
          3.315446180e-09,
         -1.138762683e-12,
          2.991214235e+03,
          4.674110790e+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
          1.017393379e+06,
         -2.509957276e+03,
          5.116547860e+00,
          1.305299930e-04,
         -8.284322260e-08,
          2.006475941e-11,
         -1.556993656e-15,
          2.019640206e+04,
         -1.101282337e+01,
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
          2.847234193e+08,
         -1.859532612e+05,
          5.008240900e+01,
         -5.142374980e-03,
          2.875536589e-07,
         -8.228817960e-12,
          9.567229020e-17,
          1.468393908e+06,
         -4.023555580e+02,
      }
   },
}
db.OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          3.99201543E+00,
         -2.40131752E-03,
          4.61793841E-06,
         -3.88113333E-09,
          1.36411470E-12,
          3.61508056E+03,
         -1.03925458E-01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          3.09288767E+00,
          5.48429716E-04,
          1.26505228E-07,
         -8.79461556E-11,
          1.17412376E-14,
          3.85865700E+03,
          4.47669610E+00,
      }
   }
}
db.OH.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  5.9711536e-01,
      B = -4.6100678e+02,
      C =  3.7606286e+04,
      D =  2.4041761e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4287721e-01,
      B = -1.8173747e+02,
      C = -8.8543767e+04,
      D =  1.9636057e+00
   },
}
db.OH.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.8627561e-01,
      B = -7.4033274e+02,
      C =  2.7559033e+04,
      D =  2.8308741e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.7918112e-01,
      B = -9.3769908e+03,
      C =  7.0509952e+06,
      D =  1.4203688e+01
   },
}
db.SF6 = {}--NONE OF SF6 DATA HAS BEEN UDPATED THIS IS MERELY A PLACED HOLDER
db.SF6.atomicConstituents = {}
db.SF6.charge = 0
db.SF6.M = {
   value = 0.146055,
   units = 'kg/mol',
}
db.SF6.gamma = {
   value = 1.3,
   note = "valid at low temperatures"
}
db.SF6.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.SF6.sutherlandVisc = {
   mu_ref = 14.8e-6, 
   T_ref = 293.15,
   S = 240.0,
   reference = "Crane Company (1988) - Flow of fluids through valves, fittings and pipes"
}
db.SF6.sutherlandThermCond = {
   T_ref = 273.0, --these have not been updated
   k_ref = 0.0241, --these have not been updated
   S = 194.0,--these have not been updated
   reference = "Table 1-3, White (2006)"
}
db.SF6.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper =  1000.0,
       coeffs = { 1.009950160e+04,--these have not been updated either
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
db.aC3H5 = {}
db.aC3H5.atomicConstituents = {C=3,H=5,}
db.aC3H5.charge = 0
db.aC3H5.M = {
   value = 41.071800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.aC3H5.gamma = {
   value = 1.1502e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.aC3H5.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.aC3H5.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.aC3H5.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.13631835E+01,
          0.19813821E-01,
          0.12497060E-04,
         -0.33355555E-07,
          0.15846571E-10,
          0.19245629E+05,
          0.17173214E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3000.0,
      coeffs = {
         0,
         0,
          0.65007877E+01,
          0.14324731E-01,
         -0.56781632E-05,
          0.11080801E-08,
         -0.90363887E-13,
          0.17482449E+05,
         -0.11243050E+02,
      }
   }
}
db.nC3H7 = {}
db.nC3H7.atomicConstituents = {C=3,H=7,}
db.nC3H7.charge = 0
db.nC3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.nC3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.nC3H7.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.nC3H7.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.nC3H7.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2, 
   segment0 ={
      T_lower = 300.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          0.10491173E+01,
          0.26008973E-01,
          0.23542516E-05,
         -0.19595132E-07,
          0.93720207E-11,
          0.10312346E+05,
          0.21136034E+02,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3000.0,
      coeffs = {
         0,
         0,
          0.77097479E+01,
          0.16031485E-01,
         -0.52720238E-05,
          0.75888352E-09,
         -0.38862719E-13,
          0.79762236E+04,
         -0.15515297E+02,
      }
   }
}
db.C2H3CHO = {}
db.C2H3CHO.atomicConstituents = {C=3,H=4,O=1,}
db.C2H3CHO.charge = 0
db.C2H3CHO.M = {   
	 value = 0.056063,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C2H3CHO.gamma = {   
	value = 1.1378,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C2H3CHO.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2H3CHO.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2H3CHO.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
	segment0 ={      
	T_lower = 298.1,      
	T_upper = 1000.0,      
	coeffs = {         
		 0.000000000000e+00,         
		 0.000000000000e+00,          
		 1.152958400000e+00,          
		 2.804021400000e-02,          
		-1.507215300000e-05,         
		 1.590584200000e-09,         
		 8.493037100000e-13,          
		-1.041769400000e+04,          
		 2.145327900000e+01,      
		}   
	},   
	segment1 = {      
	T_lower = 1000.0,      
	T_upper = 3000.0,      
	coeffs = {         
		 0.000000000000e+00,         
		 0.000000000000e+00,          
		 4.835318000000e+00,          
		 1.977260100000e-02,          
		-1.042662800000e-05,         
		 2.652580300000e-09,         
		-2.627820700000e-13,          
		-1.155783700000e+04,          
		 1.885314400000e+00,      
		}   
	}
}

