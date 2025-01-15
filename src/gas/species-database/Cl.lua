db.Cl = {}
db.Cl.type = "atom"
db.Cl.atomicConstituents = {Cl=1}
db.Cl.charge = 0
db.Cl.M = {
   value = 0.035453,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.Cl.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Cl.sigma = {
   value = 3.613e-10,
   units = 'm',
   description = 'Lennard-Jones potential distance',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
db.Cl.epsilon = {
   value = 130.8,
   units = 'K',
   description = 'Lennard-Jones potential distance',
   reference = 'Svehla (1962), NASA Technical Report R-132'
}
db.Cl.H = {
   value = 121301.000,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.Cl.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
       2.276215854e+04,
      -2.168413293e+02,
       2.745185115e+00,
       2.451101694e-03,
      -5.458011990e-06,
       4.417986880e-09,
      -1.288134004e-12,
       1.501357068e+04,
       3.102963457e+00
   },
   segment1 = {
      -1.697329300e+05,
       6.081726460e+02,
       2.128664090e+00,
       1.307367034e-04,
      -2.644883596e-08,
       2.842504775e-12,
      -1.252911731e-16,
       9.934387400e+03,
       8.844772103e+00
   },
   segment2 = {
      -7.139687070e+07,
       4.499936330e+04,
      -9.264315350e+00,
       1.657437964e-03,
      -1.326219399e-07,
       5.533998870e-12,
      -8.390301878e-17,
      -3.405333030e+05,
       1.069111426e+02
   }
}

-- viscosity and thermal conductivity values curve fit from
-- Svehla 1952, NASA Technical report TR-132
-- Robert Watt Jun 18, 2025
db.Cl.ceaViscosity = {
   nsegments = 1,
   segment0 = {
      T_lower = 100.0,
      T_upper = 5000.0,
      A = 6.37186628e-01,
      B = -2.81707549e+01,
      C = -2.20693192e+03,
      D  = 1.73547094e+00
   },
}
db.Cl.ceaThermCond = {
   nsegments = 1,
   segment0 = {
      T_lower = 100.0,
      T_upper = 5000.0,
      A = 6.17917016e-01,
      B = -8.38073113e-01,
      C = -6.98395092e+03,
      D = 1.76653095e+00
   },
}
