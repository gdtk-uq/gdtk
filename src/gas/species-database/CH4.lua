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
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      5.14987613E+00,
     -1.36709788E-02,
      4.91800599E-05,
     -4.84743026E-08,
      1.66693956E-11,
     -1.02466476E+04,
     -4.64130376E+00,
   },
   segment1 = {
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
db.CH4.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.766850998e+05,  
      2.786181020e+03, 
     -1.202577850e+01,
      3.917619290e-02, 
     -3.619054430e-05,  
      2.026853043e-08,
     -4.976705490e-12, 
     -2.331314360e+04,  
      8.904322750e+01
   },
   segment1 = { 
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



db.CH4.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.251804784118e+01,
      B = 3.644445849650e+00,
      C = -3.959557829050e-01,
      D = 1.742101397839e-02,
   }
}
db.CH4.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = 4.563206231014e+00,
      B = -6.029228026260e+00,
      C = 1.167469250228e+00,
      D = -6.185250865309e-02,
   }
}

db.CH4.Hf = {
   value = -74600.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
