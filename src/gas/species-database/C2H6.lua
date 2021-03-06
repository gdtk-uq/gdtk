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
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.29142492E+00,
     -5.50154270E-03,
      5.99438288E-05,
     -7.08466285E-08,
      2.68685771E-11,
     -1.15222055E+04,
      2.66682316E+00,
   },
   segment1 = {
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
db.C2H6.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.862044e+05,
      3.406192e+03,
     -1.951705e+01,
      7.565836e-02,
     -8.204173e-05,
      5.061136e-08,
     -1.319282e-11,
     -2.702933e+04,
      1.298140e+02,
   },
   segment1 = {
      5.025782e+06,
     -2.033022e+04,
      3.322553e+01,
     -3.836703e-03,
      7.238406e-07,
     -7.319182e-11,
      3.065469e-15,
      1.115964e+05,
     -2.039411e+02
   }
}
db.C2H6.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.55619461E00,
      B = -0.15265690E03,
      C = 0.56050805E04,
      D = 0.18241467E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.65422199E00,
      B = 0.51041684E02,
      C = -0.51534435E05,
      D = 0.10006480E01,
    },
}
db.C2H6.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.87089937E00,
      B = -0.45633731E03,
      C = 0.31766620E05,
      D = 0.16351124E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.47062424E00,
      B = -0.96911156E03,
      C = 0.10907459E06,
      D = 0.48272647E01,
    },
}

db.C2H6.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.591309461762e+01,
      B = 4.747559388706e+00,
      C = -5.145597542839e-01,
      D = 2.154750958005e-02,
   }
}
db.C2H6.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.650116871377e+01,
      B = 2.218308188236e+00,
      C = 8.894351541917e-02,
      D = -1.530913269396e-02,
   }
}

db.C2H6.Hf = {
   value = -83851.544,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
