db.Cl2 = {}
db.Cl2.type = "molecule"
db.Cl2.molecule_type = "linear"
db.Cl2.theta_v = {
   805.3,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Imported from Eilmer3'
}
db.Cl2.atomicConstituents = {Cl=2}
db.Cl2.charge = 0
db.Cl2.M = {
   value = 70.9060e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.Cl2.gamma = {
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.Cl2.sigma = {
   value = 4.217e-10,
   units = 'm',
   description = 'Lennard-Jones potential distance',
   reference = 'Svehla (1962), NASA Techncial Report R-132'
}
db.Cl2.epsilon = {
   value = 316.000,
   units = 'K',
   description = 'Lennard-Jones potential distance.',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
db.Cl2.r_eq = {
   value = 1.987e-10,
   units = 'm',
   description = 'equilibrium intermolecular distance',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
db.Cl2.SSH_mass_factor = {
   value = 1.0,
   units = "unitless",
   description = 'Mass factor = ( M ( Ma^2 + Mb^2 ) / ( 2 Ma Mb ( Ma + Mb ) )',
   reference = 'Thivet et al (1991) Phys. Fluids A 3 (11)'
}
db.Cl2.Hf = {
   value = 0.0,
   units = 'J/kg',
   description = 'heat of formation',
   reference = 'from CEA2::thermo.inp'
}
db.Cl2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
       3.462815170e+04,
      -5.547126520e+02,
       6.207589370e+00,
      -2.989632078e-03, 
       3.173027290e-06,
      -1.793629562e-09,
       4.260043590e-13,
       1.534069331e+03,
      -9.438331107e+00
   },
   segment1 = {
       6.092569420e+06,
      -1.949627662e+04,
       2.854535795e+01,
      -1.449968764e-02,
       4.463890770e-06,
      -6.358525860e-10,
       3.327360290e-14,
       1.212117724e+05,
      -1.690778824e+02
   }
}
db.Cl2.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A = 0.53516134e+00,
      B = -0.23624735e+03,
      C = 0.13738454e+05,
      D = 0.24970463e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.63348430e+00,
      B = -0.38786240e+02,
      C = -0.35830615e+05,
      D = 0.16699633e+01
   }
}
db.Cl2.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A = 0.34156262e+00,
      B = -0.46059166e+03,
      C = 0.34712872e+05,
      D = 0.37412367e+01
   },
   segment1 = {      
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A = 0.87392526e+00,
      B = 0.19876120e+03,
      C = -0.28784264e+05,
      D = -0.53204988e+00
   }
}
