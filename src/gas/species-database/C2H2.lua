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
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0,
      0,
      8.08681094E-01,
      2.33615629E-02,
     -3.55171815E-05,
      2.80152437E-08,
     -8.50072974E-12,
      2.64289807E+04,
      1.39397051E+01,
   },
   segment1 = {
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
db.C2H2.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      1.598112e+05,
     -2.216644e+03,
      1.265708e+01,
     -7.979651e-03,
      8.054993e-06,
     -2.433308e-09,
     -7.529233e-14,
      3.712619e+04,
     -5.244339e+01,
   },
   segment1 = {
      1.713847e+06,
     -5.929107e+03,
      1.236128e+01,
      1.314187e-04,
     -1.362764e-07,
      2.712656e-11,
     -1.302066e-15,
      6.266579e+04,
     -5.818961e+01
    }
}
db.C2H2.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.54922881E00,
      B = -0.17078109E03,
      C = 0.72130467E04,
      D = 0.19955795E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.65338952E00,
      B = 0.50419792E02,
      C = -0.56910493E05,
      D = 0.11190694E01,
    },
}
db.C2H2.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.72408606E00,
      B = -0.27145126E03,
      C = 0.11112107E05,
      D = 0.21630756E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.65646287E00,
      B = -0.43191905E03,
      C = 0.24326887E05,
      D = 0.27779508E01,
    },
}

db.C2H2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.530512903628e+01,
      B = 4.636835323190e+00,
      C = -5.132351164976e-01,
      D = 2.202694954932e-02,
   }
}
db.C2H2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.164258776526e+01,
      B = 5.580471585361e+00,
      C = -5.444293706233e-01,
      D = 2.047637535653e-02,
   }
}

db.C2H2.Hf = {
   value = 228200.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
