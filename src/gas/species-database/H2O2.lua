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
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -9.279533580e+04,
      1.564748385e+03,
     -5.976460140e+00,
      3.270744520e-02,
     -3.932193260e-05,
      2.509255235e-08,
     -6.465045290e-12,
     -2.494004728e+04,
      5.877174180e+01,
   },
   segment1 = {
      1.489428027e+06,
     -5.170821780e+03,
      1.128204970e+01,
     -8.042397790e-05,
     -1.818383769e-08,
      6.947265590e-12,
     -4.827831900e-16,
      1.418251038e+04,
     -4.650855660e+01,
   },
}
db.H2O2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.27611269E+00,
     -5.42822417E-04,
      1.67335701E-05,
     -2.15770813E-08,
      8.62454363E-12,
     -1.77025821E+04,
      3.43505074E+00,
   },
   segment1 = {
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
db.H2O2.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.99686871E00,
      B = -0.41461068E02,
      C = 0.87172900E04,
      D = -0.15770256E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.57419481E00,
      B = -0.50408983E03,
      C = 0.48898234E05,
      D = 0.17621537E01,
    },
}
db.H2O2.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.11075595E01,
      B = -0.20746382E03,
      C = 0.23930396E05,
      D = -0.12685243E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.46981213E00,
      B = -0.11937657E04,
      C = 0.22076993E06,
      D = 0.39203830E01,
    },
}

db.H2O2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.997002062318e+01,
      B = 2.898260266195e+00,
      C = -3.018510914973e-01,
      D = 1.346540112432e-02,
   }
}
db.H2O2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -9.226813248395e+00,
      B = 7.253630189137e-01,
      C = 1.013665976059e-01,
      D = -8.210973412477e-03,
   }
}

db.H2O2.Hf = {
   value = -135880.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
