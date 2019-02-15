db.C2 = {}
db.C2.atomicConstituents = {C=2,}
db.C2.charge = 0
db.C2.M = { 
   value = 24.0214000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'diatom -- assumed 1.4 at low temperatures'
}
db.C2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = {  5.559634510e+05,
                 -9.980126440e+03,
                  6.681620370e+01,
                 -1.743432724e-01,
                  2.448523051e-04,
                 -1.703467580e-07,
                  4.684527730e-11,
                  1.445869634e+05,
                 -3.448229700e+02
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = { -9.689267930e+05,
                  3.561092990e+03,
                 -5.064138930e-01,
                  2.945154879e-03,
                 -7.139441190e-07,
                  8.670657250e-11,
                 -4.076906810e-15,
                  7.681796830e+04,
                  3.339985240e+01
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = {  6.315145920e+06,
                  1.365420661e+04,
                 -3.996903670e+00,
                  1.937561376e-03,
                 -1.584446580e-07,
                  5.520861660e-12,
                 -7.253735340e-17,
                  9.387024990e+03,
                  6.614329920e+01
      }
   },
   ref="from CEA2::thermo.inp"
}
db.C2.ceaViscosity = {
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
   ref="McBride et al (1983)"
}
db.C2.ceaThermCond = {
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
    ref="McBride et al (1983)"
}


