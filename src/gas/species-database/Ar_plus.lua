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
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -5.731209170e+04,
      7.930791470e+02,
     -1.717121217e+00,
      1.044184018e-02,
     -1.180207501e-05,
      6.528134780e-09,
     -1.447558130e-12,
      1.790572230e+05,
      2.949150950e+01,
   },
   segment1 = {
     -3.835965400e+05,
      8.162019700e+02,
      2.301342628e+00,
     -4.952983770e-06,
      1.205108477e-08,
     -2.185050286e-12,
      1.265493898e-16,
      1.771811455e+05,
      7.947507480e+00,
   },
   segment2 = {
      1.006884827e+07,
     -6.624361280e+03,
      4.446908200e+00,
     -3.017567664e-04,
      2.612882069e-08,
     -1.201637769e-12,
      2.299206903e-17,
      2.349504137e+05,
     -1.032262257e+01,
   },
   ref="from CEA2::thermo.inp"
}

--No CEA transport data for Ar_plus, just use Ar
db['Ar+'].ceaViscosity = db.Ar.ceaViscosity
db['Ar+'].ceaThermCond = db.Ar.ceaThermCond

db['Ar+'].Hf = {
   value = 1526778.407,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
