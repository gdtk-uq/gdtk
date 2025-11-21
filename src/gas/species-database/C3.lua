db.C3 = {}
db.C3.type = "molecule"
db.C3.molecule_type = "linear"
db.C3.atomicConstituents = {C=3,}
db.C3.charge = 0
db.C3.M = { 
   value = 36.0321000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C3.gamma = { 
   value = 1.2452406707995793,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'diatom -- assumed 1.4 at low temperatures'
}
db.C3.sigma = {
   value = 3.245,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Park, 2001'
}
db.C3.epsilon = {
   value = 535.3,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Park, 2001'
}
db.C3.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0,20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -4.354614480e+04,
      6.660183220e+02,
      1.451033157e+00,
      7.434513120e-03,
     -3.810152990e-06,
     -2.336961396e-11,
      4.407054530e-13,
      9.635170200e+04,
      2.025173297e+01
   },
   segment1 = {
      4.508098930e+06,
     -1.461033761e+04,
      2.281974644e+01,
     -8.544340610e-03,
      2.146069341e-06,
     -2.103867761e-10,
      6.351589060e-15,
      1.911976065e+05,
     -1.271869723e+02
   },
   segment2 = {
      1.539589859e+08,
     -2.089057498e+05,
      7.681111210e+01,
     -8.939056190e-03,
      5.594033240e-07,
     -1.743774353e-11,
      2.181541208e-16,
      1.650801763e+06,
     -6.081693320e+02
   },
   ref="from CEA2::thermo.inp"
}
db.C3.ceaViscosity = {
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
   ref="Copied from C2: No Data Available (NNG)"
}
db.C3.ceaThermCond = {
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
   ref="Copied from C2: No Data Available (NNG)"
}


db.C3.Hf = {
   value = 839948.646,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
