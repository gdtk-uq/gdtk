-- Helium ported from Rowan's He.lua file in the cfcfd3 collection
-- PJ, 2017-05-24
db.He = {}
db.He.type = "atom"
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
db.He.sigma = {
   value = 2.551,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
db.He.epsilon = {
   value = 10.22007017,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
db.He.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00, 
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00, 
      0.000000000e+00,
     -7.453750000e+02, 
      9.287239740e-01
   },
   segment1 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00, 
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00, 
      0.000000000e+00,
     -7.453750000e+02, 
      9.287239740e-01
   },
   segment2 = {
      3.396845420e+06,
     -2.194037652e+03,
      3.080231878e+00,
     -8.068957550e-05,
      6.252784910e-09,
     -2.574990067e-13,
      4.429960218e-18,
      1.650518960e+04,
     -4.048814390e+00
   },
   reference = 'cea2::thermo.inp'
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
db.He.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
