-- Neon ported from Chris's Ne.lua file in the cfcfd3 collection
-- Yu Liu, 2018-06-08
-- LJ units fixed by NNG, 2021-02-15
db.Ne = {}
db.Ne.atomicConstituents = {Ne=1,}
db.Ne.charge = 0
db.Ne.M = {
   value = 20.1797000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.Ne.sigma = {
   value = 2.82,
   units = 'Angstrom',
   description = 'Hard sphere collision diameter',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
db.Ne.epsilon = {
   value = 32.7999876,
   units = 'K',
   description = 'Depth of the intermolecular potential minimum',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
db.Ne.gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Ne.ceaThermoCoeffs = {
   nsegments = 2,
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
      3.355322720e+00
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
      3.355322720e+00
   },
   segment2 = {
     -1.238252746e+07,
      6.958579580e+03, 
      1.016709287e+00,
      1.424664555e-04,
     -4.803933930e-09,
     -1.170213183e-13,
      8.415153652e-18,
     -5.663933630e+04,
      1.648438697e+01
   },
   reference = 'cea2::thermo.inp'
}
db.Ne.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A=0.68398511E+00, 
      B=0.18732366E+02, 
      C=-0.23663189E+04, 
      D=0.18284755E+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A=0.72333495E+00,
      B=0.10420872E+03, 
      C=-0.25429545E+05, 
      D=0.14942434E+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A=0.77549350E+00, 
      B=0.59414850E+03, 
      C=-0.69670786E+06, 
      D=0.97885712E+00
   },
   reference = "from CEA2::trans.inp which cites Bich et al. (1990)"
}
db.Ne.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A=0.68509965E+00, 
      B=0.19794924E+02, 
      C=-0.24525539E+04, 
      D=0.22586136E+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A=0.72278122E+00, 
      B=0.10528290E+03, 
      C=-0.26355706E+05, 
      D=0.19367337E+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A=0.77589413E+00, 
      B=0.61283778E+03, 
      C=-0.74015705E+06, 
      D=0.14114011E+01
   },
   reference = "from CEA2::trans.inp which cites Bich et al. (1990)"
}
db.Ne.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
