-- Curator: Rowan J. Gollan
-- Date: 08-Mar-2015
--
-- History:
--   08-Mar-2015 -- first cut.
--               -- Experiment to see if this form has longevity.
--   24-Nov-2015 -- Split into separate files.
--               -- Introduced a "defaults" table which is just
--               -- the properties for air.


db = {}

db.default = {}
db.default.atomicConstituents = {}
db.default.charge = 0
db.default.M = {
   value = 28.96e-3,
   units = 'kg/mol',
}
db.default.gamma = {
   value = 1.4,
   note = "valid at low temperatures for diatomic molecules"
}
db.default.sigma = {
   value = 3.621,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'taken from N2.'
}
db.default.epsilon = {
   value = 97.530,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'taken from N2.'
}
db.default.Lewis = {
   value = 1.152,
   reference = 'taken from N2'
}
db.default.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.default.sutherlandVisc = {
   mu_ref = 1.716e-5, 
   T_ref = 273.0,
   S = 111.0,
   reference = "Table 1-2, White (2006)"
}
db.default.sutherlandThermCond = {
   T_ref = 273.0, 
   k_ref = 0.0241, 
   S = 194.0,
   reference = "Table 1-3, White (2006)"
}
db.default.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper =  1000.0,
       coeffs = { 1.009950160e+04,
                 -1.968275610e+02,
                  5.009155110e+00,
	         -5.761013730e-03,
                  1.066859930e-05,
                 -7.940297970e-09,
                  2.185231910e-12,
                 -1.767967310e+02,
                 -3.921504225e+00 },
   },
   segment1 = { 
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {  2.415214430e+05,
                 -1.257874600e+03,
                  5.144558670e+00,
                 -2.138541790e-04,
                  7.065227840e-08,
                 -1.071483490e-11,
                  6.577800150e-16,
                  6.462263190e+03,
                 -8.147411905e+00 }
  }
}
-- CEA viscosity uses N2 values
db.default.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200,
      T_upper = 1000,
      A = 0.62526577,
      B = -31.779652,
      C = -1640.7983,
      D = 1.7454992,
    },
   segment1 = {
      T_lower = 1000,
      T_upper = 5000,
      A = 0.87395209,
      B = 561.52222,
      C = -173948.09,
      D = -0.39335958,
    },
   segment2 = {
      T_lower = 5000,
      T_upper = 15000,
      A = 0.88503551,
      B = 909.02171,
      C = -731290.61,
      D = -0.53503838,
    },
}
-- CEA thermal conductivity uses N2 values
db.default.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200,
      T_upper = 1000,
      A = 0.85439436,
      B = 105.73224,
      C = -12347.848,
      D = 0.47793128,
    },
    segment1 =  {
      T_lower = 1000,
      T_upper = 5000,
      A = 0.88407146,
      B = 133.57293,
      C = -11429.64,
      D = 0.24417019,
    },
    segment2 = {
      T_lower = 5000,
      T_upper = 15000,
      A = 2.4176185,
      B = 8047.7749,
      C = 3105580.2,
      D = -14.517761,
    },
}


db.A1 = {}
db.A1.atomicConstituents = {C=6,H=6,O=0,}
db.A1.charge = 0
db.A1.M = {   
	 value = 0.078112,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.A1.gamma = {   
	value = 1.1129,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.A1.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.A1.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.A1.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -4.899868000000e+00,          
           5.980693200000e-02,          
          -3.671008700000e-05,         
           3.274039900000e-09,         
           3.760088600000e-12,          
           9.182457000000e+03,          
           4.409564200000e+01,      
        },   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.724699400000e+01,          
           3.842016400000e-03,          
           8.277623200000e-06,         
          -4.896112000000e-09,         
           7.606454500000e-13,          
           2.664605500000e+03,          
          -7.194517500000e+01,      
	}
}

db["A1-"] = {}
db["A1-"].atomicConstituents = {C=6,H=5,O=0,}
db["A1-"].charge = 0
db["A1-"].M = {   
	value = 0.077104,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["A1-"].gamma = {   
	value = 1.1202,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["A1-"].sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["A1-"].epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["A1-"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -4.907614700000e+00,          
           5.979077100000e-02,          
          -4.563982700000e-05,         
           1.496499300000e-08,         
          -9.176782600000e-13,          
           3.873341000000e+04,          
           4.656778000000e+01,      
        },   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.449343900000e+01,          
           7.571268800000e-03,          
           3.789454200000e-06,         
          -3.076950000000e-09,         
           5.134782000000e-13,          
           3.318997700000e+04,          
          -5.428894000000e+01,      
	}
}

db.Ar = {}
db.Ar.atomicConstituents = {Ar=1,}
db.Ar.charge = 0
db.Ar.M = {
   value = 0.03994800,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.Ar.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Ar.sigma = {
   value = 3.330,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.Ar.epsilon = {
   value = 136.500,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.Ar.Lewis = {
   value = 1.173
}
db.Ar.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.Ar.ceaThermoCoeffs = {
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
      4.379674910e+00,
   },
   segment1 = {
      2.010538475e+01,
     -5.992661070e-02,
      2.500069401e+00,
     -3.992141160e-08,
      1.205272140e-11,
     -1.819015576e-15,
      1.078576636e-19,
     -7.449939610e+02,
      4.379180110e+00,
   },
   segment2 = {
     -9.951265080e+08,
      6.458887260e+05,
     -1.675894697e+02,
      2.319933363e-02,
     -1.721080911e-06,
      6.531938460e-11,
     -9.740147729e-16,
     -5.078300340e+06,
      1.465298484e+03,
   },
}
db.Ar.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0,
      0,
      0.02500000E+02,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
     -0.07453750E+04,
      0.04366000E+02,
   },
   segment1 = {
      0,
      0,
      0.02500000E+02,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
     -0.07453750E+04,
      0.04366000E+02,
   }
}
db.Ar.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.1205763e-01,
      B = -6.7714354e+01,
      C =  1.9040660e+02,
      D =  2.1588272e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.9357334e-01,
      B =  7.0953943e+01,
      C = -2.8386007e+04,
      D =  1.4856447e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.6608935e-01,
      B =  6.7867215e+02,
      C = -8.4991417e+05,
      D =  7.7935167e-01
   },
}
db.Ar.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0968928e-01,
      B = -7.0892249e+01,
      C =  5.8420624e+02,
      D =  1.9337152e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.9075463e-01,
      B =  6.2676058e+01,
      C = -2.5667413e+04,
      D =  1.2664189e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.6269502e-01,
      B =  6.2341752e+02,
      C = -7.1899552e+05,
      D =  5.6927918e-01
   },
}
db.Ar.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
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
db.C = {}
db.C.atomicConstituents = {C=1,}
db.C.type = "atom"
db.C.charge = 0
db.C.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.C.gamma = { 
   value = 1.664,
   units = 'non-dimensional',
    description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.C.sigma = {
   value = 3.298,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C.epsilon = {
   value = 71.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      6.495031470e+02,
     -9.649010860e-01,
      2.504675479e+00,
     -1.281448025e-05,
      1.980133654e-08,
     -1.606144025e-11,
      5.314483411e-15,
      8.545763110e+04,
      4.747924288e+00
   },
   segment1 = {
     -1.289136472e+05,
      1.719528572e+02,
      2.646044387e+00,
     -3.353068950e-04,
      1.742092740e-07,
     -2.902817829e-11,
      1.642182385e-15,
      8.410597850e+04,
      4.130047418e+00
   },
   segment2 = {
      4.432528010e+08,
     -2.886018412e+05,
      7.737108320e+01,
     -9.715281890e-03,
      6.649595330e-07,
     -2.230078776e-11,
      2.899388702e-16,
      2.355273444e+06,
     -6.405123160e+02
   }
}
db.C.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0,
      0,
      2.55423955E+00,
     -3.21537724E-04,
      7.33792245E-07,
     -7.32234889E-10,
      2.66521446E-13,
      8.54438832E+04,
      4.53130848E+00,
   },
   segment1 = {
      0,
      0,
      2.49266888E+00,
      4.79889284E-05,
     -7.24335020E-08,
      3.74291029E-11,
     -4.87277893E-15,
      8.54512953E+04,
      4.80150373E+00,
   }
}
db.C.Hf = {
   value = 716680.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['C+'] = {}
db['C+'].atomicConstituents = {C=1,}
db['C+'].charge = 1
db['C+'].M = { 
   value = 12.0101514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db['C+'].gamma = { 
   value = 1.664,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'copied from C.lua'
}
db['C+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      2.258535929e+03,
     -1.574575687e+00,
      2.503637730e+00,
     -5.202878370e-06,
      4.516908390e-09,
     -2.181431053e-12,
      4.495047033e-16,
      2.168951913e+05,
      4.345699505e+00
   },
   segment1 = {
      1.255112551e+04,
     -3.411874670e+01,
      2.543383218e+00,
     -2.805120849e-05,
      9.751641970e-09,
     -1.736855394e-12,
      1.246191931e-16,
      2.171001786e+05, 
      4.063913515e+00
   },
   segment2 = {
      5.618135320e+05,
     -6.047058900e+03,
      5.884541470e+00,
     -7.211894530e-04,
      6.823484110e-08,
     -2.599878590e-12,
      3.633868358e-17,
      2.581370458e+05,
     -2.280019759e+01
   }
}
db['C+'].Hf = {
   value = 1809444.482,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
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
db.C2.sigma = {
   value = 3.913,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
db.C2.epsilon = {
   value = 195.1013396,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Svehla (1962) NASA Technical Report R-132'
}
db.C2.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0,20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      5.559634510e+05,
     -9.980126440e+03,
      6.681620370e+01,
     -1.743432724e-01,
      2.448523051e-04,
     -1.703467580e-07,
      4.684527730e-11,
      1.445869634e+05,
     -3.448229700e+02
   },
   segment1 = {
     -9.689267930e+05,
      3.561092990e+03,
     -5.064138930e-01,
      2.945154879e-03,
     -7.139441190e-07,
      8.670657250e-11,
     -4.076906810e-15,
      7.681796830e+04,
      3.339985240e+01
   },
   segment2 = {
      6.315145920e+06,
      1.365420661e+04,
     -3.996903670e+00,
      1.937561376e-03,
     -1.584446580e-07,
      5.520861660e-12,
     -7.253735340e-17,
      9.387024990e+03,
      6.614329920e+01
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


db.C2.Hf = {
   value = 830457.322,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
--Auto-Generated File
db.C2H = {}
db.C2H.atomicConstituents = {C=2,H=1,O=0,Ar=0,He=0,N=0}
db.C2H.charge = 0
db.C2H.M = {
    value = 0.025029,
    units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H.gamma = {
   value = 1.246482,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = ' '
}
db.C2H.epsilon = {
   value = 209.000000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = ' '
}
db.C2H.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      1.343669487e+04,
     -5.067970720e+02,
      7.772107410e+00,
     -6.512339820e-03,
      1.030117855e-05,
     -5.880147670e-09,
      1.226901861e-12,
      6.892269990e+04,
     -1.871881626e+01
   },
   segment1 = {
      3.922334570e+06,
     -1.204751703e+04,
      1.756172920e+01,
     -3.655442940e-03,
      6.987685430e-07,
     -6.825162010e-11,
      2.719262793e-15,
      1.433266627e+05,
     -9.561634380e+01
   }
}
db.C2H.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0.000000000000e+00,
      0.000000000000e+00,
      2.889657330000e+00,
      1.340996110000e-02,
     -2.847695010000e-05,
      2.947910450000e-08,
     -1.093315110000e-11,
      6.683939320000e+04,
      6.222964380000e+00,
   },
   segment1 = {
      0.000000000000e+00,
      0.000000000000e+00,
      3.167806520000e+00,
      4.752219020000e-03,
     -1.837870770000e-06,
      3.041902520000e-10,
     -1.772327700000e-14,
      6.712106500000e+04,
      6.635894750000e+00,
   }
}
db.C2H.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.532486828957e+01,
      B = 4.636835182216e+00,
      C = -5.132350959604e-01,
      D = 2.202694856015e-02,
   }
}
db.C2H.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -4.671028652157e+00,
      B = -1.552634627462e+00,
      C = 4.455792117341e-01,
      D = -2.549920860076e-02,
   }
}
db.C2H.Hf = {
   value = 566200.482,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
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
db.C2H3 = {}
db.C2H3.atomicConstituents = {C=2,H=3,}
db.C2H3.charge = 0
db.C2H3.M = {
   value = 27.045220e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H3.gamma = {
   value = 1.2408e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H3.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H3.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0,
      0,
      3.21246645E+00,
      1.51479162E-03,
      2.59209412E-05,
     -3.57657847E-08,
      1.47150873E-11,
      3.48598468E+04,
      8.51054025E+00,
   },
   segment1 = {
      0,
      0,
      3.01672400E+00,
      1.03302292E-02,
     -4.68082349E-06,
      1.01763288E-09,
     -8.62607041E-14,
      3.46128739E+04,
      7.78732378E+00,
   }
}
db.C2H3.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -3.347897e+04,
      1.064104e+03,
     -6.403857e+00,
      3.934515e-02,
     -4.760046e-05,
      3.170071e-08,
     -8.633406e-12,
      3.039123e+04,
      5.809226e+01,
   },
   segment1 = {
      2.718080e+06,
     -1.030957e+04,
      1.836580e+01,
     -1.580131e-03,
      2.680595e-07,
     -2.439004e-11,
      9.209096e-16,
      9.765056e+04,
     -9.760087e+01
    }
}
db.C2H3.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.528613897362e+01,
      B = 4.636835198248e+00,
      C = -5.132350983080e-01,
      D = 2.202694867379e-02,
   }
}
db.C2H3.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.726898602161e+01,
      B = 3.149807052372e+00,
      C = -1.248614364768e-01,
      D = -2.248078033156e-03,
   }
}

db.C2H3.Hf = {
   value = 299686.817,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C2H3CHO = {}
db.C2H3CHO.atomicConstituents = {C=3,H=4,O=1,}
db.C2H3CHO.charge = 0
db.C2H3CHO.M = {   
	 value = 0.056063,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C2H3CHO.gamma = {   
	value = 1.1378,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C2H3CHO.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2H3CHO.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2H3CHO.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {298.15, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.152958400000e+00,          
           2.804021400000e-02,          
          -1.507215300000e-05,         
           1.590584200000e-09,         
           8.493037100000e-13,          
          -1.041769400000e+04,          
           2.145327900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           4.835318000000e+00,          
           1.977260100000e-02,          
          -1.042662800000e-05,         
           2.652580300000e-09,         
          -2.627820700000e-13,          
          -1.155783700000e+04,          
           1.885314400000e+00,      
	}
}

db.C2H4 = {}
db.C2H4.atomicConstituents = {C=2,H=4,}
db.C2H4.charge = 0
db.C2H4.M = {
   value = 28.053160e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H4.gamma = {
   value = 1.2393e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H4.sigma = {
   value = 3.971,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H4.epsilon = {
   value = 280.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H4.Lewis = {
   value = 1.402
}
db.C2H4.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.95920148E+00,
     -7.57052247E-03,
      5.70990292E-05,
     -6.91588753E-08,
      2.69884373E-11,
      5.08977593E+03,
      4.09733096E+00,
   },
   segment1 = {
      0,
      0,
      2.03611116E+00,
      1.46454151E-02,
     -6.71077915E-06,
      1.47222923E-09,
     -1.25706061E-13,
      4.93988614E+03,
      1.03053693E+01,
   }
}
db.C2H4.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.163606e+05,
      2.554852e+03,
     -1.609746e+01,
      6.625779e-02,
     -7.885082e-05,
      5.125225e-08,
     -1.370340e-11,
     -6.176191e+03,
      1.093338e+02,
   },
   segment1 = {
      3.408764e+06,
     -1.374848e+04,
      2.365898e+01,
     -2.423804e-03,
      4.431396e-07,
     -4.352683e-11,
      1.775411e-15,
      8.820429e+04,
     -1.371278e+02
   }
}
db.C2H4.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.55243600E00,
      B = -0.16260917E03,
      C = 0.64734038E04,
      D = 0.19463233E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.65385054E00,
      B = 0.51157317E02,
      C = -0.54731184E05,
      D = 0.10933538E01,
    },
}
db.C2H4.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.77957663E00,
      B = -0.47857623E03,
      C = 0.32147858E05,
      D = 0.21827872E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.48277394E00,
      B = -0.91773465E03,
      C = 0.11528060E06,
      D = 0.45824405E01,
    },
}

db.C2H4.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.580836606635e+01,
      B = 4.671473377835e+00,
      C = -4.951477412808e-01,
      D = 2.031861716368e-02,
   }
}
db.C2H4.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.876932653720e+01,
      B = 3.279081810443e+00,
      C = -7.600615678041e-02,
      D = -7.052525306715e-03,
   }
}

db.C2H4.Hf = {
   value = 52500.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C2H5 = {}
db.C2H5.atomicConstituents = {C=2,H=5,}
db.C2H5.charge = 0
db.C2H5.M = {
   value = 29.061100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H5.gamma = {
   value = 1.1963e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H5.sigma = {
   value = 4.302,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H5.epsilon = {
   value = 252.300,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H5.Lewis = {
   value = 1.551
}
db.C2H5.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.30646568E+00,
     -4.18658892E-03,
      4.97142807E-05,
     -5.99126606E-08,
      2.30509004E-11,
      1.28416265E+04,
      4.70720924E+00,
   },
   segment1 = {
      0,
      0,
      1.95465642E+00,
      1.73972722E-02,
     -7.98206668E-06,
      1.75217689E-09,
     -1.49641576E-13,
      1.28575200E+04,
      1.34624343E+01,
   }
}
db.C2H5.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.411313e+05,
      2.714285e+03,
     -1.534978e+01,
      6.451673e-02,
     -7.259144e-05,
      4.599116e-08,
     -1.218368e-11,
      5.981419e+02,
      1.090967e+02,
   },
   segment1 = {
      4.169220e+06,
     -1.662982e+04,
      2.795442e+01,
     -3.051716e-03,
      5.685160e-07,
     -5.682864e-11,
      2.355649e-15,
      1.137010e+05,
     -1.639358e+02
   }
}

db.C2H5.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.593014146501e+01,
      B = 4.747559219867e+00,
      C = -5.145597297061e-01,
      D = 2.154750839715e-02,
   }
}
db.C2H5.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.527203989367e+01,
      B = 1.873182769765e+00,
      C = 1.110291998328e-01,
      D = -1.533035825556e-02,
   }
}

db.C2H5.Hf = {
   value = 118658.24,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C2H6 = {}
db.C2H6.atomicConstituents = {C=2,H=6,}
db.C2H6.charge = 0
db.C2H6.M = {
   value = 30.069040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C2H6.gamma = {
   value = 1.1872e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C2H6.sigma = {
   value = 4.302,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H6.epsilon = {
   value = 252.300,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C2H6.Lewis = {
   value = 1.546
}
db.C2H6.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.29142492E+00,
     -5.50154270E-03,
      5.99438288E-05,
     -7.08466285E-08,
      2.68685771E-11,
     -1.15222055E+04,
      2.66682316E+00,
   },
   segment1 = {
      0,
      0,
      1.07188150E+00,
      2.16852677E-02,
     -1.00256067E-05,
      2.21412001E-09,
     -1.90002890E-13,
     -1.14263932E+04,
      1.51156107E+01,
   }
}
db.C2H6.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.862044e+05,
      3.406192e+03,
     -1.951705e+01,
      7.565836e-02,
     -8.204173e-05,
      5.061136e-08,
     -1.319282e-11,
     -2.702933e+04,
      1.298140e+02,
   },
   segment1 = {
      5.025782e+06,
     -2.033022e+04,
      3.322553e+01,
     -3.836703e-03,
      7.238406e-07,
     -7.319182e-11,
      3.065469e-15,
      1.115964e+05,
     -2.039411e+02
   }
}
db.C2H6.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.55619461E00,
      B = -0.15265690E03,
      C = 0.56050805E04,
      D = 0.18241467E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.65422199E00,
      B = 0.51041684E02,
      C = -0.51534435E05,
      D = 0.10006480E01,
    },
}
db.C2H6.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.87089937E00,
      B = -0.45633731E03,
      C = 0.31766620E05,
      D = 0.16351124E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.47062424E00,
      B = -0.96911156E03,
      C = 0.10907459E06,
      D = 0.48272647E01,
    },
}

db.C2H6.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.591309461762e+01,
      B = 4.747559388706e+00,
      C = -5.145597542839e-01,
      D = 2.154750958005e-02,
   }
}
db.C2H6.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.650116871377e+01,
      B = 2.218308188236e+00,
      C = 8.894351541917e-02,
      D = -1.530913269396e-02,
   }
}

db.C2H6.Hf = {
   value = -83851.544,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C2O = {}
db.C2O.atomicConstituents = {C=2,H=0,O=1,}
db.C2O.charge = 0
db.C2O.M = {   
	 value = 0.040021,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C2O.gamma = {   
	value = 1.2386,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C2O.sigma = {   
	value =  3.828 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2O.epsilon = {   
	value = 232.400 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           3.368851000000e+00,          
           8.241803000000e-03,          
          -8.765145000000e-06,         
           5.569262000000e-09,         
          -1.540009000000e-12,          
           3.317081000000e+04,          
           6.713314000000e+00,      
        },
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           4.849809000000e+00,          
           2.947585000000e-03,          
          -1.090729000000e-06,         
           1.792562000000e-10,         
          -1.115758000000e-14,          
           3.282055000000e+04,          
          -6.453226000000e-01,      
	}
}

db.C2O.Hf = {
   value = 291038.666,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C3H2 = {}
db.C3H2.atomicConstituents = {C=3,H=2,O=0,}
db.C3H2.charge = 0
db.C3H2.M = {   
	 value = 0.038048,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C3H2.gamma = {   
	value = 1.177,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C3H2.sigma = {   
	value =  4.100 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H2.epsilon = {   
	value = 209.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.691077000000e+00,          
           1.480366000000e-02,          
          -3.250551000000e-06,         
          -8.644363000000e-09,         
           5.284878000000e-12,          
           5.219072000000e+04,          
           8.757391000000e+00,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.530853000000e+00,          
           5.870316000000e-03,          
          -1.720777000000e-06,         
           2.127498000000e-10,         
          -8.291910000000e-15,          
           5.115214000000e+04,          
          -1.122728000000e+01,      
	}
}

db.C3H3 = {}
db.C3H3.atomicConstituents = {C=3,H=3,O=0,}
db.C3H3.charge = 0
db.C3H3.M = {   
	 value = 0.039056,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C3H3.gamma = {   
	value = 1.1464,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C3H3.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H3.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H3.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.351109270000e+00,          
           3.274112230000e-02,          
          -4.738271350000e-05,         
           3.763098080000e-08,         
          -1.185409230000e-11,          
           4.010577830000e+04,          
           1.520589240000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           7.142218800000e+00,          
           7.619020050000e-03,          
          -2.674599500000e-06,         
           4.249148010000e-10,         
          -2.514754150000e-14,          
           3.890874270000e+04,          
          -1.258484360000e+01,      
	}
}

db.C3H6 = {}
db.C3H6.atomicConstituents = {C=3,H=6,}
db.C3H6.charge = 0
db.C3H6.M = {
   value = 42.079740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H6.gamma = {
   value = 1.1474e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H6.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.C3H6.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.C3H6.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 ={
      0,
      0,
      0.01493307E+02,
      0.02092518E+00,
      0.04486794E-04,
     -0.01668912E-06,
      0.07158146E-10,
      0.01074826E+05,
      0.01614534E+03,
   },
   segment1 = {
      0,
      0,
      0.06732257E+02,
      0.01490834E+00,
     -0.04949899E-04,
      0.07212022E-08,
     -0.03766204E-12,
     -0.09235703E+04,
     -0.01331335E+03,
   }
}
db.C3H6.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.912462e+05,
      3.542074e+03,
     -2.114879e+01,
      8.901485e-02,
     -1.001429e-04,
      6.267959e-08,
     -1.637871e-11,
     -1.529962e+04,
      1.407641e+02,
   },
   segment1 = {
      5.017620e+06,
     -2.086084e+04,
      3.644156e+01,
     -3.881191e-03,
      7.278677e-07,
     -7.321204e-11,
      3.052176e-15,
      1.261245e+05,
     -2.195716e+02
   }
}
db.C3H6.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.52642893E00,
      B = -0.24304494E03,
      C = 0.14490001E05,
      D = 0.21036650E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.64243372E00,
      B = 0.83055174E01,
      C = -0.61290810E05,
      D = 0.11264132E01,
    },
}
db.C3H6.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.75434495E00,
      B = -0.56817108E03,
      C = 0.39706666E05,
      D = 0.23579094E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.46796950E00,
      B = -0.98032164E03,
      C = 0.12025017E06,
      D = 0.46607118E01,
    },
}

db.C3H6.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -2.622497446663e+01,
      B = 4.773760930715e+00,
      C = -5.124555504018e-01,
      D = 2.122189081184e-02,
   }
}
db.C3H6.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -3.193319260373e+01,
      B = 8.804340373459e+00,
      C = -8.583432800778e-01,
      D = 2.968532569582e-02,
   }
}

db.C3H7 = {}
db.C3H7.atomicConstituents = {C=3,H=7,}
db.C3H7.charge = 0
db.C3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H7.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H7.epsilon = {
   value = 266.8,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H7.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.10515518E+01,
      0.25991980E-01,
      0.23800540E-05,
     -0.19609569E-07,
      0.93732470E-11,
      0.10631863E+05,
      0.21122559E+02,
   },
   segment1 = {
      0,
      0,
      0.77026987E+01,
      0.16044203E-01,
     -0.52833220E-05,
      0.76298590E-09,
     -0.39392284E-13,
      0.82984336E+04,
     -0.15480180E+02,
   }
}
db.C3H8 = {}
db.C3H8.atomicConstituents = {C=3,H=8,}
db.C3H8.charge = 0
db.C3H8.M = {
   value = 44.095620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.C3H8.gamma = {
   value = 1.1267e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.C3H8.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H8.epsilon = {
   value = 266.80,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.C3H8.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.93355381E+00,
      0.26424579E-01,
      0.61059727E-05,
     -0.21977499E-07,
      0.95149253E-11,
     -0.13958520E+05,
      0.19201691E+02,
   },
   segment1 = {
      0,
      0,
      0.75341368E+01,
      0.18872239E-01,
     -0.62718491E-05,
      0.91475649E-09,
     -0.47838069E-13,
     -0.16467516E+05,
     -0.17892349E+02,
   }
}
db.C3H8.Hf = {
   value = -104680.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C4H = {}
db.C4H.atomicConstituents = {C=4,H=1,O=0,}
db.C4H.charge = 0
db.C4H.M = {   
	 value = 0.049051,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H.gamma = {   
	value = 1.1418,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.318629500000e+00,          
           3.858295600000e-02,          
          -7.138562300000e-05,         
           6.535635900000e-08,         
          -2.261766600000e-11,          
           9.545610600000e+04,          
           1.556758300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           7.769759300000e+00,          
           4.982997600000e-03,          
          -1.762854600000e-06,         
           2.814428400000e-10,         
          -1.668986900000e-14,          
           9.434590000000e+04,          
          -1.416527400000e+01,      
	}
}

db.C4H2 = {}
db.C4H2.atomicConstituents = {C=4,H=2,O=0,}
db.C4H2.charge = 0
db.C4H2.M = {   
	 value = 0.050059,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H2.gamma = {   
	value = 1.1268,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H2.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H2.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -3.920103000000e-01,          
           5.193756500000e-02,          
          -9.173734000000e-05,         
           8.047198600000e-08,         
          -2.689821800000e-11,          
           5.484526600000e+04,          
           2.095779400000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           8.663770800000e+00,          
           6.724718900000e-03,          
          -2.359339700000e-06,         
           3.750638000000e-10,         
          -2.223094000000e-14,          
           5.325227500000e+04,          
          -2.109350300000e+01,      
	}
}

db.C4H4 = {}
db.C4H4.atomicConstituents = {C=4,H=4,O=0,}
db.C4H4.charge = 0
db.C4H4.M = {   
	 value = 0.052074,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H4.gamma = {   
	value = 1.1281,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H4.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H4.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.915247900000e+00,          
           5.275087800000e-02,          
          -7.165594400000e-05,         
           5.507242300000e-08,         
          -1.728622800000e-11,          
           3.297850400000e+04,          
           3.141998300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.650709200000e+00,          
           1.612943400000e-02,          
          -7.193887500000e-06,         
           1.498178700000e-09,         
          -1.186411000000e-13,          
           3.119599200000e+04,          
          -9.795211800000e+00,      
	}
}

db.C4H6 = {}
db.C4H6.atomicConstituents = {C=4,H=6,O=0,}
db.C4H6.charge = 0
db.C4H6.M = {   
	 value = 0.05409,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H6.gamma = {   
	value = 1.1216,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H6.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H6.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H6.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.128446500000e-01,          
           3.436902200000e-02,          
          -1.110739200000e-05,         
          -9.210666000000e-09,         
           6.206517900000e-12,          
           1.180227000000e+04,          
           2.308999600000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           8.867313400000e+00,          
           1.491867000000e-02,          
          -3.154871600000e-06,         
          -4.184133000000e-10,         
           1.576125800000e-13,          
           9.133851600000e+03,          
          -2.332817100000e+01,      
	}
}

db.C4H612 = {}
db.C4H612.atomicConstituents = {C=4,H=6,O=0,}
db.C4H612.charge = 0
db.C4H612.M = {   
	 value = 0.05409,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H612.gamma = {   
	value = 1.1148,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H612.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H612.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H612.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.023467000000e+00,          
           3.495919000000e-02,          
          -2.200905000000e-05,         
           6.942272000000e-09,         
          -7.879187000000e-13,          
           1.811799000000e+04,          
           1.975066000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.781557000000e+01,          
          -4.257502000000e-03,          
           1.051185000000e-05,         
          -4.473844000000e-09,         
           5.848138000000e-13,          
           1.267342000000e+04,          
          -6.982662000000e+01,      
	}
}

db.C4H7 = {}
db.C4H7.atomicConstituents = {C=4,H=7,O=0,}
db.C4H7.charge = 0
db.C4H7.M = {   
	 value = 0.055098,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H7.gamma = {   
	value = 1.1079,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H7.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H7.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H7.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.869825400000e-01,          
           3.696449500000e-02,          
          -8.627744100000e-06,         
          -1.505182100000e-08,         
           8.989126300000e-12,          
           2.055130100000e+04,          
           2.448446700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.196339200000e+01,          
           1.142530500000e-02,          
           7.894890900000e-07,         
          -1.985887200000e-09,         
           3.687364500000e-13,          
           1.696297700000e+04,          
          -3.754290800000e+01,      
	}
}

db.C4H81 = {}
db.C4H81.atomicConstituents = {C=4,H=8,O=0,}
db.C4H81.charge = 0
db.C4H81.M = {   
	 value = 0.056106,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H81.gamma = {   
	value = 1.1073,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H81.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H81.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H81.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.181138000000e+00,          
           3.085338000000e-02,          
           5.086524700000e-06,         
          -2.465488800000e-08,         
           1.111019300000e-11,          
          -1.790400400000e+03,          
           2.106246900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.053584100000e+00,          
           3.435050700000e-02,          
          -1.588319700000e-05,         
           3.308966200000e-09,         
          -2.536104500000e-13,          
          -2.139723100000e+03,          
           1.554320100000e+01,      
	}
}

db.C5H4O = {}
db.C5H4O.atomicConstituents = {C=5,H=4,O=1,}
db.C5H4O.charge = 0
db.C5H4O.M = {   
	 value = 0.080084,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H4O.gamma = {   
	value = 1.1211,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H4O.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H4O.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H4O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -2.391535500000e+00,          
           4.736368000000e-02,          
          -3.072817100000e-05,         
           7.803155200000e-09,         
          -2.514572900000e-13,          
           4.374015200000e+03,          
           3.459433700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           4.792724200000e+00,          
           2.922168000000e-02,          
          -1.599748600000e-05,         
           4.206904900000e-09,         
          -4.281517900000e-13,          
           2.284928600000e+03,          
          -3.013189300000e+00,      
	}
}

db.C5H4OH = {}
db.C5H4OH.atomicConstituents = {C=5,H=5,O=1,}
db.C5H4OH.charge = 0
db.C5H4OH.M = {   
	 value = 0.081092,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H4OH.gamma = {   
	value = 1.0946,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H4OH.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H4OH.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H4OH.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.282223600000e+00,          
           4.904116000000e-02,          
          -1.368899700000e-05,         
          -2.913385800000e-08,         
           1.900696400000e-11,          
           8.008709800000e+03,          
           3.079835800000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.336791200000e+01,          
           1.520578500000e-02,          
          -5.459225800000e-06,         
           8.813486600000e-10,         
          -5.277445400000e-14,          
           3.841150600000e+03,          
          -4.592083900000e+01,      
	}
}

db.C5H5 = {}
db.C5H5.atomicConstituents = {C=5,H=5,O=0,}
db.C5H5.charge = 0
db.C5H5.M = {   
	 value = 0.065093,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H5.gamma = {   
	value = 1.1209,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H5.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -9.590371800000e-01,          
           3.139685900000e-02,          
           2.672379400000e-05,         
          -6.894187200000e-08,         
           3.330185600000e-11,          
           3.072912000000e+04,          
           2.907281600000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.084406600000e+01,          
           1.539283700000e-02,          
          -5.563042100000e-06,         
           9.018937100000e-10,         
          -5.415653100000e-14,          
           2.690056600000e+04,          
          -3.525494800000e+01,      
	}
}

db.C5H5O = {}
db.C5H5O.atomicConstituents = {C=5,H=5,O=1,}
db.C5H5O.charge = 0
db.C5H5O.M = {   
	 value = 0.081092,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H5O.gamma = {   
	value = 1.1011,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H5O.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5O.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.304283500000e-01,          
           3.232269100000e-02,          
           2.890044300000e-05,         
          -7.067997700000e-08,         
           3.340689100000e-11,          
           8.075308200000e+03,          
           2.533097400000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.260642200000e+01,          
           1.674726000000e-02,          
          -6.109857400000e-06,         
           9.967655700000e-10,         
          -6.011320100000e-14,          
           3.931345500000e+03,          
          -4.260427700000e+01,      
	}
}

db.C5H6 = {}
db.C5H6.atomicConstituents = {C=5,H=6,O=0,}
db.C5H6.charge = 0
db.C5H6.M = {   
	 value = 0.066101,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H6.gamma = {   
	value = 1.1228,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H6.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H6.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H6.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {298.15, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -2.897895800000e+00,          
           4.348477700000e-02,          
          -3.351100500000e-06,         
          -3.110375600000e-08,         
           1.691244400000e-11,          
           1.508474200000e+04,          
           3.689476000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.062432000000e+01,          
           1.773544800000e-02,          
          -6.233044600000e-06,         
           9.730831700000e-10,         
          -5.550013000000e-14,          
           1.077218800000e+04,          
          -3.577342200000e+01,      
	}
}

db.C6H2 = {}
db.C6H2.atomicConstituents = {C=6,H=2,O=0,}
db.C6H2.charge = 0
db.C6H2.M = {   
	 value = 0.07408,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H2.gamma = {   
	value = 1.0872,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H2.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H2.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.593262400000e+00,          
           8.053014500000e-02,          
          -1.480064900000e-04,         
           1.330003100000e-07,         
          -4.533231300000e-11,          
           8.327322700000e+04,          
           2.798087300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.322628100000e+01,          
           7.390430200000e-03,          
          -2.271538100000e-06,         
           2.587521700000e-10,         
          -5.535674100000e-15,          
           8.056525800000e+04,          
          -4.120117600000e+01,      
	}
}

db.C6H2.Hf = {
   value = 670000.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.C6H3 = {}
db.C6H3.atomicConstituents = {C=6,H=3,O=0,}
db.C6H3.charge = 0
db.C6H3.M = {   
	 value = 0.075088,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H3.gamma = {   
	value = 1.0866,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H3.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H3.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H3.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.179061900000e+00,          
           5.554736000000e-02,          
          -7.307616800000e-05,         
           5.207673600000e-08,         
          -1.504696400000e-11,          
           8.564731200000e+04,          
           1.917919900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.818834300000e+00,          
           2.793340800000e-02,          
          -1.782542700000e-05,         
           5.370253600000e-09,         
          -6.170762700000e-13,          
           8.518825000000e+04,          
          -9.214782700000e-01,      
	}
}

db.C6H5O = {}
db.C6H5O.atomicConstituents = {C=6,H=5,O=1,}
db.C6H5O.charge = 0
db.C6H5O.M = {   
	 value = 0.093103,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H5O.gamma = {   
	value = 1.0959,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H5O.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5O.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.821943300000e+00,          
           4.812251000000e-02,          
          -4.679230200000e-06,         
          -3.401859400000e-08,         
           1.864963700000e-11,          
           4.242918000000e+03,          
           3.352619900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.383398400000e+01,          
           1.761840300000e-02,          
          -6.069625700000e-06,         
           9.198817300000e-10,         
          -5.044918100000e-14,          
          -6.921254900000e+02,          
          -5.039299000000e+01,      
	}
}

db.C6H5OH = {}
db.C6H5OH.atomicConstituents = {C=6,H=6,O=1,}
db.C6H5OH.charge = 0
db.C6H5OH.M = {   
	 value = 0.094111,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H5OH.gamma = {   
	value = 1.0867,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H5OH.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5OH.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5OH.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.695653900000e+00,          
           5.227129900000e-02,          
          -7.202405000000e-06,         
          -3.585960300000e-08,         
           2.044907300000e-11,          
          -1.328412100000e+04,          
           3.254216000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.491207300000e+01,          
           1.837813500000e-02,          
          -6.198312800000e-06,         
           9.198322100000e-10,         
          -4.920956500000e-14,          
          -1.837519900000e+04,          
          -5.592410300000e+01,      
	}
}

db.CH = {}
db.CH.atomicConstituents = {C=1,H=1,}
db.CH.charge = 0
db.CH.M = {
   value = 13.0186400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH.sigma = {
   value = 2.750,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH.epsilon = {
   value = 80.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      2.220590133e+04,
     -3.405411530e+02,
      5.531452290e+00,
     -5.794964260e-03,
      7.969554880e-06,
     -4.465911590e-09,
      9.596338320e-13,
      7.240783270e+04,
     -9.107673050e+00
   },
   segment1 = {
      2.060763440e+06,
     -5.396206660e+03,
      7.856293850e+00,
     -7.965907450e-04,
      1.764308305e-07,
     -1.976386267e-11,
      5.030429510e-16,
      1.062236592e+05,
     -3.154757439e+01
   },
   segment2 = {
     -8.068368690e+08,
      4.575450540e+05,
     -9.843975080e+01,
      1.235244098e-02,
     -8.485608570e-07,
      3.040410624e-11,
     -4.400315170e-16,
     -3.595851590e+06,
      8.953477440e+02
   }
}
db.CH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.48981665E+00,
      3.23835541E-04,
     -1.68899065E-06,
      3.16217327E-09,
     -1.40609067E-12,
      7.07972934E+04,
      2.08401108E+00,
   },
   segment1 = {
      0,
      0,
      2.87846473E+00,
      9.70913681E-04,
      1.44445655E-07,
     -1.30687849E-10,
      1.76079383E-14,
      7.10124364E+04,
      5.48497999E+00,
   }
}
db.CH.ceaViscosity = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.63306702E00,
      B = -0.75491296E01,
      C = -0.21396736E04,
      D = 0.14781500E01,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.64807346E00,
      B = 0.31141665E01,
      C = -0.18292566E04,
      D = 0.13638041E01,
    },
}
db.CH.ceaThermCond = {
   notes = 'GORDON;NASATM86885,OCT1984',
   nsegments = 2,
   segment0 = {
      T_lower = 300.000,
      T_upper = 1000.000,
      A = 0.10884807E01,
      B = 0.27220319E03,
      C = -0.22480091E05,
      D = -0.64679028E00,
    },
   segment1 = {
      T_lower = 1000.000,
      T_upper = 5000.000,
      A = 0.44721563E00,
      B = -0.10513211E04,
      C = 0.28092259E06,
      D = 0.48038688E01,
    },
}

db.CH.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.795623025033e+01,
      B = 2.112205943159e+00,
      C = -1.984976383157e-01,
      D = 8.935785511923e-03,
   }
}
db.CH.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = 4.285614607933e+00,
      B = -4.158100803522e+00,
      C = 6.955243193806e-01,
      D = -3.195851658452e-02,
   }
}

db.CH.Hf = {
   value = 597370.604,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CH2 = {}
db.CH2.atomicConstituents = {C=1,H=2,}
db.CH2.charge = 0
db.CH2.M = {
   value = 14.0265800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH2.gamma = {
   value = 1.311,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH2.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2.Lewis = {
   value = 1.023
}
db.CH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.76267867E+00,
      9.68872143E-04,
      2.79489841E-06,
     -3.85091153E-09,
      1.68741719E-12,
      4.60040401E+04,
      1.56253185E+00,
   },
   segment1 = {
      0,
      0,
      2.87410113E+00,
      3.65639292E-03,
     -1.40894597E-06,
      2.60179549E-10,
     -1.87727567E-14,
      4.62636040E+04,
      6.17119324E+00,
   }
}
db.CH2.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      3.218922e+04,
     -2.877602e+02,
      4.203584e+00,
      3.455406e-03,
     -6.746193e-06,
      7.654572e-09,
     -2.870328e-12,
      4.733625e+04,
     -2.143629e+00,
   },
   segment1 = {
      2.550418e+06,
     -7.971625e+03,
      1.228924e+01,
     -1.699123e-03,
      2.991729e-07,
     -2.767007e-11,
      1.051342e-15,
      9.642217e+04,
     -6.094740e+01
   }
}

db.CH2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.274340688184e+01,
      B = 3.690876678686e+00,
      C = -4.016022859623e-01,
      D = 1.764923161171e-02,
   }
}
db.CH2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.125379834838e+00,
      B = -2.235463527624e+00,
      C = 5.038887351549e-01,
      D = -2.646402375981e-02,
   }
}

db.CH2.Hf = {
   value = 390364.517,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['CH2*'] = {}
db['CH2*'].atomicConstituents = {C=1,H=2,}
db['CH2*'].charge = 0
db['CH2*'].M = {
   value = 14.026580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db['CH2*'].gamma = {
   value = 1.3263e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db['CH2*'].sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db['CH2*'].epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}

db['CH2*'].grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.19860411E+00,
     -2.36661419E-03,
      8.23296220E-06,
     -6.68815981E-09,
      1.94314737E-12,
      5.04968163E+04,
     -7.69118967E-01,
   },
   segment1 = {
      0,
      0,
      2.29203842E+00,
      4.65588637E-03,
     -2.01191947E-06,
      4.17906000E-10,
     -3.39716365E-14,
      5.09259997E+04,
      8.62650169E+00,
   }
}
-- Same as GRIMech
db['CH2*'].ceaThermoCoeffs = db['CH2*'].grimechThermoCoeffs

db['CH2*'].chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.274340688184e+01,
      B = 3.690876678686e+00,
      C = -4.016022859623e-01,
      D = 1.764923161171e-02,
   }
}
db['CH2*'].chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = 3.563134631131e+00,
      B = -4.824793791738e+00,
      C = 8.889411555766e-01,
      D = -4.524261095664e-02,
   }
}

db.CH2CHO = {}
db.CH2CHO.atomicConstituents = {C=2,H=3,O=1,}
db.CH2CHO.charge = 0
db.CH2CHO.M = {
   value = 43.044620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CHO.gamma = {
   value = 1.1776e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CHO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CHO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.03409062E+02,
      0.10738574E-01,
      0.01891492E-04,
     -0.07158583E-07,
      0.02867385E-10,
      0.15214766E+04,
      0.09558290E+02,
   },
   segment1 = {
      0,
      0,
      0.05975670E+02,
      0.08130591E-01,
     -0.02743624E-04,
      0.04070304E-08,
     -0.02176017E-12,
      0.04903218E+04,
     -0.05045251E+02,
   }
}
-- Same as GRIMech
db.CH2CHO.ceaThermoCoeffs = db.CH2CHO.grimechThermoCoeffs

db.CH2CHO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -2.720702096187e+01,
      B = 4.986894409273e+00,
      C = -5.038737864236e-01,
      D = 1.949669378655e-02,
   }
}
db.CH2CHO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -2.577625433343e+01,
      B = 6.373198044217e+00,
      C = -5.427655805016e-01,
      D = 1.604699040543e-02,
   }
}

db.CH2CO = {}
db.CH2CO.atomicConstituents = {C=2,H=2,O=1,}
db.CH2CO.charge = 0
db.CH2CO.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CO.gamma = {
   value = 1.1908e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.13583630E+00,
      1.81188721E-02,
     -1.73947474E-05,
      9.34397568E-09,
     -2.01457615E-12,
     -7.04291804E+03,
      1.22156480E+01,
   },
   segment1 = {
      0,
      0,
      4.51129732E+00,
      9.00359745E-03,
     -4.16939635E-06,
      9.23345882E-10,
     -7.94838201E-14,
     -7.55105311E+03,
      6.32247205E-01,
   }
}
db.CH2CO.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      3.549598e+04,
     -4.063063e+02,
      3.718922e+00,
      1.583502e-02,
     -1.726196e-05,
      1.157377e-08,
     -3.305843e-12,
     -5.209993e+03,
      3.839604e+00,
   },
   segment1 = {
      2.013565e+06,
     -8.200887e+03,
      1.759694e+01,
     -1.464545e-03,
      2.695887e-07,
     -2.665675e-11,
      1.094205e-15,
      4.177777e+04,
     -8.725804e+01
   }
}

db.CH2CO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.193913326467e+01,
      B = 2.737948202373e+00,
      C = -1.862634718958e-01,
      D = 4.621015645441e-03,
   }
}
db.CH2CO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.799368137710e+01,
      B = 3.225542866669e+00,
      C = -1.223022955920e-01,
      D = -2.848863216417e-03,
   }
}

db.CH2O = {}
db.CH2O.atomicConstituents = {C=1,H=2,O=1}
db.CH2O.charge = 0
db.CH2O.M = {
   value = 30.025980e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2O.gamma = {
   value = 1.3065e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2O.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2O.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2O.Lewis = {
   value = 1.329
}
db.CH2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0.000000000e+00, 
      0.000000000e+00, 
      4.793723150e+00, 
     -9.908333690e-03, 
      3.732200080e-05, 
     -3.792852610e-08, 
      1.317726520e-11, 
     -1.430895670e+04, 
      6.028129000e-01, 
  },
  segment1 = {
     0.000000000e+00, 
     0.000000000e+00, 
     1.760690080e+00, 
     9.200000820e-03, 
    -4.422588130e-06, 
     1.006412120e-09, 
    -8.838556400e-14, 
    -1.399583230e+04, 
     1.365632300e+01, 
   }
}
-- Same as GRIMech
db.CH2O.ceaThermoCoeffs = db.CH2O.grimechThermoCoeffs

db.CH2O.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.944420616305e+01,
      B = 1.583460638757e+00,
      C = -1.418900703767e-02,
      D = -3.702993576843e-03,
   }
}
db.CH2O.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -3.458688239117e-01,
      B = -4.719241554972e+00,
      C = 1.048457032022e+00,
      D = -5.935855918886e-02,
   }
}

db.CH2O.Hf = {
   value = -108580.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CH2OH = {}
db.CH2OH.atomicConstituents = {C=1,H=3,O=1,}
db.CH2OH.charge = 0
db.CH2OH.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2OH.gamma = {
   value = 1.2070e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2OH.sigma = {
   value = 3.690,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2OH.epsilon = {
   value = 417.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.86388918E+00,
      5.59672304E-03,
      5.93271791E-06,
     -1.04532012E-08,
      4.36967278E-12,
     -3.19391367E+03,
      5.47302243E+00,
   },
   segment1 = {
      0,
      0,
      3.69266569E+00,
      8.64576797E-03,
     -3.75101120E-06,
      7.87234636E-10,
     -6.48554201E-14,
     -3.24250627E+03,
      5.81043215E+00,
   }
}
db.CH2OH.Hf = {
   value = -17800.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CH2_S = {}
db.CH2_S.atomicConstituents = {C=1,H=2,}
db.CH2_S.charge = 0
db.CH2_S.M = {
   value = 14.026580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2_S.gamma = {
   value = 1.3263e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2_S.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2_S.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2_S.Lewis = {
   value = 1.022
}
db.CH2_S.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.19860411E+00,
     -2.36661419E-03,
      8.23296220E-06,
     -6.68815981E-09,
      1.94314737E-12,
      5.04968163E+04,
     -7.69118967E-01,
   },
   segment1 = {
      0,
      0,
      2.29203842E+00,
      4.65588637E-03,
     -2.01191947E-06,
      4.17906000E-10,
     -3.39716365E-14,
      5.09259997E+04,
      8.62650169E+00,
   }
}
db.CH3 = {}
db.CH3.atomicConstituents = {C=1,H=3,}
db.CH3.charge = 0
db.CH3.M = {
   value = 15.0345200e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH3.gamma = {
   value = 1.276,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH3.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.Lewis = {
   value = 1.049
}
db.CH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.67359040E+00,
      2.01095175E-03,
      5.73021856E-06,
     -6.87117425E-09,
      2.54385734E-12,
      1.64449988E+04,
      1.60456433E+00,
   },
   segment1 = {
      0,
      0,
      2.28571772E+00,
      7.23990037E-03,
     -2.98714348E-06,
      5.95684644E-10,
     -4.67154394E-14,
      1.67755843E+04,
      8.48007179E+00,
   }
}
db.CH3.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -2.876189e+04,
      5.093269e+02,
      2.002144e-01,
      1.363606e-02,
     -1.433989e-05,
      1.013557e-08,
     -3.027332e-12,
      1.408272e+04,
      2.022773e+01,
   },
   segment1 = {
      2.760803e+06,
     -9.336531e+03,
      1.487730e+01,
     -1.439430e-03,
      2.444478e-07,
      -2.224556e-11,
      8.395066e-16,
      7.481809e+04,
     -7.919682e+01
   }
}

db.CH3.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.270870995825e+01,
      B = 3.690876383840e+00,
      C = -4.016022429639e-01,
      D = 1.764922953865e-02,
   }
}
db.CH3.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.243228188431e-01,
      B = -3.501395522712e+00,
      C = 7.418831266650e-01,
      D = -3.973086278935e-02,
   }
}

db.CH3.Hf = {
   value = 146658.04,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CH3CCH2 = {}
db.CH3CCH2.atomicConstituents = {C=3,H=5,O=0,}
db.CH3CCH2.charge = 0
db.CH3CCH2.M = {   
	 value = 0.041072,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.CH3CCH2.gamma = {   
	value = 1.1463,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.CH3CCH2.sigma = {   
	value =  4.982 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CCH2.epsilon = {   
	value = 266.800 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CCH2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.732920900000e+00,          
           2.239462000000e-02,          
          -5.149061100000e-06,         
          -6.759646600000e-09,         
           3.825321100000e-12,          
           2.904049800000e+04,          
           1.656887800000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.425552800000e+00,          
           1.551107200000e-02,          
          -5.667835000000e-06,         
           7.922438800000e-10,         
          -1.687803400000e-14,          
           2.784302700000e+04,          
          -3.352718400000e+00,      
	}
}

db.CH3CHCH = {}
db.CH3CHCH.atomicConstituents = {C=3,H=5,O=0,}
db.CH3CHCH.charge = 0
db.CH3CHCH.M = {   
	 value = 0.041072,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.CH3CHCH.gamma = {   
	value = 1.1482,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.CH3CHCH.sigma = {   
	value =  4.982 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CHCH.epsilon = {   
	value = 266.800 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CHCH.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           9.137293100000e-01,          
           2.643234300000e-02,          
          -1.175895000000e-05,         
          -2.303567800000e-09,         
           2.771548800000e-12,          
           3.091686700000e+04,          
           1.998926900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.372528100000e+00,          
           1.578050900000e-02,          
          -5.992285000000e-06,         
           9.308966400000e-10,         
          -3.655096600000e-14,          
           2.961476000000e+04,          
          -3.418647800000e+00,      
	}
}

db.CH3CHO = {}
db.CH3CHO.atomicConstituents = {C=2,H=4,O=1,}
db.CH3CHO.charge = 0
db.CH3CHO.M = {
   value = 44.052560e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3CHO.gamma = {
   value = 1.1762e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3CHO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3CHO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.47294595E+01,
     -0.31932858E-02,
      0.47534921E-04,
     -0.57458611E-07,
      0.21931112E-10,
     -0.21572878E+05,
      0.41030159E+01,
   },
   segment1 = {
      0,
      0,
      0.54041108E+01,
      0.11723059E-01,
     -0.42263137E-05,
      0.68372451E-09,
     -0.40984863E-13,
     -0.22593122E+05,
     -0.34807917E+01,
   }
}
db.CH3CHO.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.373904e+05,
      2.559938e+03,
     -1.340470e+01,
      5.922129e-02,
     -6.240006e-05,
      3.703324e-08,
     -9.342697e-12,
     -3.318731e+04,
      1.007418e+02,
   },
   segment1 = {
      3.321177e+06,
     -1.449720e+04,
      2.708421e+01,
     -2.879320e-03,
      5.556310e-07,
     -5.732675e-11,
      2.443965e-15,
      6.507756e+04,
     -1.536236e+02
    }
}
db.CH3CHO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 6000.000,
      A = -2.464943462372e+01,
      B = 3.963220206531e+00,
      C = -3.676625277378e-01,
      D = 1.349561598291e-02,
   }
}
db.CH3CHO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 6000.000,
      A = -2.472743688385e+01,
      B = 5.846139483528e+00,
      C = -4.559596622026e-01,
      D = 1.169834480741e-02,
   }
}

db.CH3CO = {}
db.CH3CO.atomicConstituents = {C=2,H=3,O=1,}
db.CH3CO.charge = 0
db.CH3CO.M = {   
	 value = 0.043044,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.CH3CO.gamma = {   
	value = 1.1951,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.CH3CO.sigma = {   
	value =  3.970 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CO.epsilon = {   
	value = 436.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CO.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           4.163425700000e+00,          
          -2.326161000000e-04,          
           3.426782000000e-05,         
          -4.410522700000e-08,         
           1.727561200000e-11,          
          -2.657452900000e+03,          
           7.346828000000e+00,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.944773100000e+00,          
           7.866720500000e-03,          
          -2.886588200000e-06,         
           4.727087500000e-10,         
          -2.859986100000e-14,          
          -3.787307500000e+03,          
          -5.013675100000e+00,      
	}
}

db.CH3O = {}
db.CH3O.atomicConstituents = {C=1,H=3,O=1,}
db.CH3O.charge = 0
db.CH3O.M = {
   value = 31.033920e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3O.gamma = {
   value = 1.2802e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3O.sigma = {
   value = 3.690,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3O.epsilon = {
   value = 417.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3O.Lewis = {
   value = 1.360
}
db.CH3O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 3000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.02106204E+02,
      0.07216595E-01,
      0.05338472E-04,
     -0.07377636E-07,
      0.02075610E-10,
      0.09786011E+04,
      0.13152177E+02,
   },
   segment1 = {
      0,
      0,
      0.03770799E+02,
      0.07871497E-01,
     -0.02656384E-04,
      0.03944431E-08,
     -0.02112616E-12,
      0.12783252E+03,
      0.02929575E+02,
   }
}
db.CH3O.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      8.657118e+04,
     -6.631685e+02,
      2.257456e+00,
      2.266284e-02,
     -2.970566e-05,
      2.199341e-08,
     -6.588043e-12,
      4.174102e+03,
      8.174778e+00,
   },
   segment1 = {
      2.101188e+06,
     -8.841969e+03,
      1.822646e+01,
     -1.743485e-03,
      3.340434e-07,
     -3.430673e-11,
      1.473898e-15,
      5.309582e+04,
     -9.422501e+01
   }
}

db.CH3O.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 3000.000,
      A = -2.596824781753e+01,
      B = 4.466550601397e+00,
      C = -4.304908052551e-01,
      D = 1.601059649628e-02,
   }
}
db.CH3O.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 3000.000,
      A = -2.187848840881e+01,
      B = 4.318601152219e+00,
      C = -1.996927776238e-01,
      D = -2.088356518496e-03,
   }
}

db.CH3O.Hf = {
   value = 13000.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CH3OH = {}
db.CH3OH.atomicConstituents = {C=1,H=4,O=1,}
db.CH3OH.charge = 0
db.CH3OH.M = {
   value = 32.041860e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH3OH.gamma = {
   value = 1.2320e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH3OH.sigma = {
   value = 3.626,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3OH.epsilon = {
   value = 481.8,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      5.71539582E+00,
     -1.52309129E-02,
      6.52441155E-05,
     -7.10806889E-08,
      2.61352698E-11,
     -2.56427656E+04,
     -1.50409823E+00,
   },
   segment1 = {
      0,
      0,
      1.78970791E+00,
      1.40938292E-02,
     -6.36500835E-06,
      1.38171085E-09,
     -1.17060220E-13,
     -2.53748747E+04,
      1.45023623E+01,
   }
}
db.CH3OH.Hf = {
   value = -200940.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
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
db.CHO = {}
db.CHO.atomicConstituents = {C=1,H=1,O=1,}
db.CHO.charge = 0
db.CHO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp',
}
db.CHO.gamma = {
   value = 1.316,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = "Gokel (2004), Dean's Handbook of Organic Chemistry",
}
db.CHO.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CHO.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CHO.Lewis = {
   value = 1.314
}
db.CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.22118584E+00,
     -3.24392532E-03,
      1.37799446E-05,
     -1.33144093E-08,
      4.33768865E-12,
      3.83956496E+03,
      3.39437243E+00,
   },
   segment1 = {
      0,
      0,
      2.77217438E+00,
      4.95695526E-03,
     -2.48445613E-06,
      5.89161778E-10,
     -5.33508711E-14,
      4.01191815E+03,
      9.79834492E+00,
   }
}
db.CN = {}
db.CN.atomicConstituents = {C=1,N=1,}
db.CN.charge = 0
db.CN.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CN.gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CN.sigma = {
   value = 3.856,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CN.epsilon = {
   value = 75.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CN.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      3.949148570e+03,
     -1.391590572e+02,
      4.930835320e+00,
     -6.304670510e-03,
      1.256836472e-05,
     -9.878300500e-09,
      2.843137221e-12,
      5.228455380e+04,
     -2.763115585e+00
   },
   segment1 = {
     -2.228006270e+06,
      5.040733390e+03,
     -2.121897722e-01,
      1.354901134e-03,
      1.325929798e-07,
     -6.937006370e-11,
      5.494952270e-15,
      1.784496132e+04,
      3.282563919e+01
   },
   segment2 = {
     -1.794798118e+08,
      1.054346069e+05,
     -1.729624170e+01,
      2.194895530e-03,
     -8.508938030e-08,
      9.318692990e-13,
      6.358139930e-18,
     -7.962594120e+05,
      1.913139639e+02
   }
}
db.CN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   nsegments = 2, 
   segment0 = {
      0,
      0,
      0.36129351E+01,
     -0.95551327E-03,
      0.21442977E-05,
     -0.31516323E-09,
     -0.46430356E-12,
      0.51708340E+05,
      0.39804995E+01,
   },
   segment1 = {
      0,
      0,
      0.37459805E+01,
      0.43450775E-04,
      0.29705984E-06,
     -0.68651806E-10,
      0.44134173E-14,
      0.51536188E+05,
      0.27867601E+01,
   }
}
db.CN.Hf = {
   value = 438683.552,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CO = {}
db.CO.atomicConstituents = {C=1,O=1,}
db.CO.type = "molecule"
db.CO.molecule_type = "linear"
db.CO.theta_v = {
   value = 3122,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Capitelli et al (2005), Table 11. omega_e in ground state converted to K'
}
db.CO.charge = 0
db.CO.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CO.gamma = {
   value = 1.3992e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CO.sigma = {
   value = 3.650,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO.epsilon = {
   value = 98.100,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO.Lewis = {
   value = 1.171
}
db.CO.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      1.489045326e+04,
     -2.922285939e+02,
      5.724527170e+00,
     -8.176235030e-03,
      1.456903469e-05,
     -1.087746302e-08,
      3.027941827e-12,
     -1.303131878e+04,
     -7.859241350e+00
   },
   segment1 = { 
      4.619197250e+05,
     -1.944704863e+03,
      5.916714180e+00,
     -5.664282830e-04,
      1.398814540e-07,
     -1.787680361e-11,
      9.620935570e-16,
     -2.466261084e+03,
     -1.387413108e+01
   },
   segment2 = {
      8.868662960e+08,
     -7.500377840e+05,
      2.495474979e+02,
     -3.956351100e-02,
      3.297772080e-06,
     -1.318409933e-10,
      1.998937948e-15,
      5.701421130e+06,
     -2.060704786e+03
   },
}

db.CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.57953347E+00,
     -6.10353680E-04,
      1.01681433E-06,
      9.07005884E-10,
     -9.04424499E-13,
     -1.43440860E+04,
      3.50840928E+00,
   },
   segment1 = {
      0,
      0,
      2.71518561E+00,
      2.06252743E-03,
     -9.98825771E-07,
      2.30053008E-10,
     -2.03647716E-14,
     -1.41518724E+04,
      7.81868772E+00,
   }
}


db.CO.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower =200.0,
      T_upper =1000.0,
      A= 0.62526577e+00, 
      B=-0.31779652e+02, 
      C=-0.16407983e+04, 
      D=0.17454992e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A= 0.87395209e+00,
      B= 0.56152222e+03,
      C= -0.17394809e+06,
      D= -0.39335958e+00
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A= 0.88503551e+00,
      B= 0.90902171e+03,
      C= -0.73129061e+06,
      D= -0.53503838e+00
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
}

db.CO.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper =1000.0, 
      A = 0.85439436e+00,
      B = 0.10573224e+03,
      C = -0.12347848e+05,
      D = 0.47793128e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A = 0.88407146e+00,
      B = 0.13357293e+03,
      C = -0.11429640e+05,
      D = 0.24417019E+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A = 0.24175411e+01,
      B = 0.80462671e+04,
      C = 0.31090740e+07,
      D = -0.14516932e+02
   },
      --ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
}


db.CO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.957998515812e+01,
      B = 2.674022664222e+00,
      C = -2.730698167015e-01,
      D = 1.223479766828e-02,
   }
}
db.CO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.707713918561e+00,
      B = -1.810436476889e+00,
      C = 4.069171041651e-01,
      D = -2.093217501136e-02,
   }
}

db.CO.Hf = {
   value = -110535.196,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['CO+'] = {}
db['CO+'].atomicConstituents = {C=1,O=1,}
db['CO+'].charge = 1
db['CO+'].M = {
   value = 28.095514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db['CO+'].gamma = {
   value = 1.3992e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db['CO+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -2.178786658e+04,
      1.288857032e+02, 
      3.769057550e+00,
     -3.431730130e-03, 
      8.193945750e-06,
     -6.463814690e-09,
      1.803727574e-12,
      1.482345898e+05,
      3.990547070e+00
   },
   segment1 = { 
      2.316847506e+05,
     -1.057646148e+03,
      4.554257780e+00,
      4.495520320e-04,
     -2.489507047e-07,
      5.267566420e-11,
     -3.289510270e-15,
      1.555050724e+05,
     -3.873462640e+00
   },
   segment2 = {
     -3.035604054e+08,
      2.393118392e+05,
     -7.034999240e+01,
      1.139551440e-02,
     -8.315173100e-07,
      2.863705515e-11,
     -3.803269410e-16,
     -1.688617704e+06,
      6.291980420e+02
   },
}
-- No CEA transport data for CO+, just use CO
db['CO+'].ceaViscosity = db.CO.ceaViscosity 
db['CO+'].ceaThermCond = db.CO.ceaThermCond 


db['CO+'].Hf = {
   value = 1247789.204,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.CO2 = {}
db.CO2.atomicConstituents = {C=1, O=2}
db.CO2.type = "molecule"
db.CO2.molecule_type = "linear"
db.CO2.theta_v = {
   -- CO2 has 4 vibrational modes, with characteristic vibrational temperatures: 960 (degeneracy 2), 1992, 3380
   -- For theta_v, use lowest value (this is what Park et al [1993] seem to do)
   value = 960,
   units = 'K',
   description = 'charateristic vibrational temperature',
   reference = 'Park et al (1993)'
}
db.CO2.charge = 0
db.CO2.M = {
   value = 0.04401,
   units = 'kg/mol',
}
db.CO2.gamma = {
   value = 1.3,
   note = "valid at low temperatures"
}
db.CO2.sigma = {
   value = 3.763,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO2.epsilon = {
   value = 244.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CO2.Lewis = {
   value = 1.404
}
db.CO2.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.CO2.sutherlandVisc = {
   mu_ref = 14.8e-6, 
   T_ref = 293.15,
   S = 240.0,
   reference = "Crane Company (1988) - Flow of fluids through valves, fittings and pipes"
}
db.CO2.sutherlandThermCond = {
   T_ref = 273.0, --these have not been updated
   k_ref = 0.0241, --these have not been updated
   S = 194.0,--these have not been updated
   reference = "Table 1-3, White (2006)"
}

db.CO2.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      4.943650540e+04,
     -6.264116010e+02,
      5.301725240e+00,
      2.503813816e-03,
     -2.127308728e-07,
     -7.689988780e-10,
      2.849677801e-13,
     -4.528198460e+04,
     -7.048279440e+00
   },
   segment1 = { 
      1.176962419e+05,
     -1.788791477e+03,
      8.291523190e+00,
     -9.223156780e-05,
      4.863676880e-09,
     -1.891053312e-12,
      6.330036590e-16,
     -3.908350590e+04,
     -2.652669281e+01
   },
   segment2 = { 
     -1.544423287e+09,
      1.016847056e+06,
     -2.561405230e+02,
      3.369401080e-02,
     -2.181184337e-06,
      6.991420840e-11,
     -8.842351500e-16,
     -8.043214510e+06,
      2.254177493e+03
   } -- from thermo.inp Gurvich, 1991 pt1 p211 pt2 p200
}

db.CO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.35677352E+00,
      8.98459677E-03,
     -7.12356269E-06,
      2.45919022E-09,
     -1.43699548E-13,
     -4.83719697E+04,
      9.90105222E+00,
   },
   segment1 = {
      0,
      0,
      3.85746029E+00,
      4.41437026E-03,
     -2.21481404E-06,
      5.23490188E-10,
     -4.72084164E-14,
     -4.87591660E+04,
      2.27163806E+00,
   }
}

db.CO2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower =200.0,
      T_upper =1000.0,
      A= 0.51137258e+00,
      B= -0.22951321e+03,
      C= 0.13710678e+05,
      D= 0.27075538e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A= 0.63978285e+00,
      B= -0.42637076e+02,
      C= -0.15522605e+05,
      D= 0.16628843e+01
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A= 0.72150912e+00,
      B= 0.75012895e+03,
      C= -0.11825507e+07,
      D= 0.85493645e+00
     -- ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
   }
}

db.CO2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper =1000.0, 
      A = 0.48056568e+00,
      B = -0.50786720e+03,
      C = 0.35088811e+05,
      D = 0.36747794e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0, 
      A = 0.69857277e+00,
      B = -0.11830477e+03,
      C = -0.50688859e+05,
      D =  0.18650551e+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0, 
      A = 0.10518358e+01,
      B = -0.42555944e+04,
      C = 0.14288688e+08,
      D = -0.88950473e+00
   }
      --ref = 'from CEA2::trans.inp which cites Bousheri et al. (1987) and Svehla (1994)'
}


db.CO2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.541411738166e+01,
      B = 4.756161515092e+00,
      C = -5.184638085888e-01,
      D = 2.183133880431e-02,
   }
}
db.CO2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.119099884033e+01,
      B = 5.186927559697e+00,
      C = -4.741229077145e-01,
      D = 1.610702319175e-02,
   }
}

db.CO2.Hf = {
   value = -393510.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.H = {}
db.H.atomicConstituents = {H=1,}
db.H.charge = 0
db.H.M = {
   value = 0.00100794,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.H.sigma = {
   value = 2.050,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H.epsilon = {
   value = 145.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H.Lewis = {
   value = 0.189
}
db.H.ceaThermoCoeffs = {
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
      2.547370801e+04,
     -4.466828530e-01,
   },
   segment1 = {
      6.078774250e+01,
     -1.819354417e-01,
      2.500211817e+00,
     -1.226512864e-07,
      3.732876330e-11,
     -5.687744560e-15,
      3.410210197e-19,
      2.547486398e+04,
     -4.481917770e-01,
   },
   segment2 = {
      2.173757694e+08,
     -1.312035403e+05,
      3.399174200e+01,
     -3.813999680e-03,
      2.432854837e-07,
     -7.694275540e-12,
      9.644105630e-17,
      1.067638086e+06,
     -2.742301051e+02,
   },
}
db.H.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.4226149e-01,
      B = -4.0132865e+02,
      C =  1.8554165e+05,
      D =  4.6741844e-02
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7486623e-01,
      B = -2.5022902e+03,
      C =  7.0955048e+06,
      D = -9.3888455e-01
   },
}
db.H.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.4166119e-01,
      B = -4.0487203e+02,
      C =  1.8775642e+05,
      D =  3.4843121e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7447639e-01,
      B = -2.5089452e+03,
      C =  7.1081294e+06,
      D =  2.4970991e+00
   },
}
db.H.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.50000000E+00,
      7.05332819E-13,
     -1.99591964E-15,
      2.30081632E-18,
     -9.27732332E-22,
      2.54736599E+04,
     -4.46682853E-01,
   },
   segment1 = {
      0,
      0,
      2.50000001E+00,
     -2.30842973E-11,
      1.61561948E-14,
     -4.73515235E-18,
      4.98197357E-22,
      2.54736599E+04,
     -4.46682914E-01,
   }
}

db.H.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.287685911498e+01,
      B = 3.709357281377e+00,
      C = -4.038665036507e-01,
      D = 1.774159302331e-02,
   }
}
db.H.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.253722034102e+01,
      B = 3.709357074341e+00,
      C = -4.038664735031e-01,
      D = 1.774159157190e-02,
   }
}

db.H.Hf = {
   value = 217998.828,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['H+'] = {}
db['H+'].atomicConstituents = {H=1}
db['H+'].charge = 1
db['H+'].M = {
   value = 1.0073914e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
--Use Gamma value for H.
db['H+'].gamma = { 
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

db['H+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      1.840214877e+05,
     -1.140646644e+00
   },
   segment1 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      1.840214877e+05,
     -1.140646644e+00
   },
   segment2 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      1.840214877e+05,
     -1.140646644e+00
   },
   ref="from CEA2::thermo.inp"
}

-- No CEA transport data for H_plus, just use H
db['H+'].ceaViscosity = db.H.ceaViscosity
db['H+'].ceaThermCond = db.H.ceaThermCond


db['H+'].Hf = {
   value = 1536245.928,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.H2 = {}
db.H2.atomicConstituents = {H=2,}
db.H2.charge = 0
db.H2.M = {
   value = 0.00201588,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.H2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.H2.sigma = {
   value = 2.920,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2.epsilon = {
   value = 38.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2.Lewis = {
   value = 0.317
}
db.H2.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      4.078323210e+04,
     -8.009186040e+02,
      8.214702010e+00,
     -1.269714457e-02,
      1.753605076e-05,
     -1.202860270e-08,
      3.368093490e-12,
      2.682484665e+03,
     -3.043788844e+01,
   },
   segment1 = {
      5.608128010e+05,
     -8.371504740e+02,
      2.975364532e+00,
      1.252249124e-03,
     -3.740716190e-07,
      5.936625200e-11,
     -3.606994100e-15,
      5.339824410e+03,
     -2.202774769e+00,
   },
   segment2 = {
      4.966884120e+08,
     -3.147547149e+05,
      7.984121880e+01,
     -8.414789210e-03,
      4.753248350e-07,
     -1.371873492e-11,
      1.605461756e-16,
      2.488433516e+06,
     -6.695728110e+02,
   },
}
db.H2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.4553182e-01,
      B =  4.3555109e+01,
      C = -3.2579340e+03,
      D =  1.3556243e-01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.6730605e-01,
      B =  6.7931897e+02,
      C = -2.1025179e+05,
      D = -1.8251697e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  1.0126129e+00,
      B =  1.4973739e+03,
      C = -1.4428484e+06,
      D = -2.3254928e+00
   },
}
db.H2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  1.0059461e+00,
      B =  2.7951262e+02,
      C = -2.9792018e+04,
      D =  1.1996252e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  1.0582450e+00,
      B =  2.4875372e+02,
      C =  1.1736907e+04,
      D =  8.2758695e-01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -2.2364420e-01,
      B = -6.9650442e+03,
      C = -7.7771313e+04,
      D =  1.3189369e+01
   },
}
db.H2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.34433112E+00,
      7.98052075E-03,
     -1.94781510E-05,
      2.01572094E-08,
     -7.37611761E-12,
     -9.17935173E+02,
      6.83010238E-01,
   },
   segment1 = {
      0,
      0,
      3.33727920E+00,
     -4.94024731E-05,
      4.99456778E-07,
     -1.79566394E-10,
      2.00255376E-14,
     -9.50158922E+02,
     -3.20502331E+00,
   }
}

db.H2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.625410533148e+01,
      B = 1.046163975892e+00,
      C = -5.543497377779e-02,
      D = 2.537276639956e-03,
   }
}
db.H2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -5.540632470662e+00,
      B = 8.558861129169e-01,
      C = -6.135027230087e-02,
      D = 5.213698814113e-03,
   }
}

db.H2.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.H2C4O = {}
db.H2C4O.atomicConstituents = {C=4,H=2,O=1,}
db.H2C4O.charge = 0
db.H2C4O.M = {   
	 value = 0.066058,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.H2C4O.gamma = {   
	value = 1.13,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.H2C4O.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.H2C4O.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.H2C4O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 4000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           4.810971000000e+00,          
           1.313999000000e-02,          
           9.865073000000e-07,         
          -6.120720000000e-09,         
           1.640003000000e-12,          
           2.545803000000e+04,          
           2.113424000000e+00,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.026888000000e+01,          
           4.896164000000e-03,          
          -4.885081000000e-07,         
          -2.708566000000e-10,         
           5.107013000000e-14,          
           2.346903000000e+04,          
          -2.815985000000e+01,      
	}
}

db.H2CC = {}
db.H2CC.atomicConstituents = {H=2,C=2,}
db.H2CC.charge = 0
db.H2CC.M = {
   value = 26.0372800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2CC.gamma = {
   value = 1.242,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2CC.sigma = {
   value = 4.100,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file'
}
db.H2CC.epsilon = {
   value = 209.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file'
}
db.H2CC.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.32815483E+01,
      0.69764791E-02,
     -0.23855244E-05,
     -0.12104432E-08,
      0.98189545E-12,
      0.48621794E+05,
      0.59203910E+01,
   },
   segment1 = {
      0,
      0,
      0.42780340E+01,
      0.47562804E-02,
     -0.16301009E-05,
      0.25462806E-09,
     -0.14886379E-13,
      0.48316688E+05,
      0.64023701E+00,
   }
}
db.H2CC.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
         -1.466042e+04,
      2.789476e+02,
      1.276230e+00,
      1.395015e-02,
     -1.475703e-05,
      9.476298e-09,
     -2.567602e-12,
      4.736110e+04,
      1.658226e+01,
   },
   segment1 = {
      1.940839e+06,
     -6.892718e+03,
      1.339582e+01,
     -9.368969e-04,
      1.470804e-07,
     -1.220040e-11,
      4.122392e-16,
      9.107113e+04,
     -6.337503e+01
    }
}

db.H2CC.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 6000.000,
      A = -2.396866153965e+01,
      B = 4.040169611213e+00,
      C = -4.252689499838e-01,
      D = 1.774234492713e-02,
   }
}
db.H2CC.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 6000.000,
      A = -1.534295035569e+01,
      B = 2.865727279816e+00,
      C = -1.600083918257e-01,
      D = 2.376301694764e-03,
   }
}

db.H2CN = {}
db.H2CN.atomicConstituents = {C=1,H=2,N=1,}
db.H2CN.charge = 0
db.H2CN.M = {
   value = 28.033280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.H2CN.gamma = {
   value = 1.2769e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.H2CN.sigma = {
   value = 3.63,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2CN.epsilon = {
   value = 569.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2CN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   T_break_points = {300.0, 1000.0, 4000.0},
   T_blend_ranges = {400.0},
   nsegments = 2, 
   segment0 ={
      0,
      0,
      0.28516610E+01,
      0.56952331E-02,
      0.10711400E-05,
     -0.16226120E-08,
     -0.23511081E-12,
      0.28637820E+05,
      0.89927511E+01,
   },
   segment1 = {
      0,
      0,
      0.52097030E+01,
      0.29692911E-02,
     -0.28555891E-06,
     -0.16355500E-09,
      0.30432589E-13,
      0.27677109E+05,
     -0.44444780E+01,
   }
}
db.H2O = {}
db.H2O.atomicConstituents = {O=1,H=2,}
db.H2O.charge = 0
db.H2O.M = {
   value = 0.01801528,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.H2O.gamma = {
   value = 1.32900000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.H2O.sigma = {
   value = 2.605,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O.epsilon = {
   value = 572.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.H2O.Lewis = {
   value = 0.854
}
db.H2O.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -3.947960830e+04,
      5.755731020e+02,
      9.317826530e-01,
      7.222712860e-03,
     -7.342557370e-06,
      4.955043490e-09,
     -1.336933246e-12,
     -3.303974310e+04,
      1.724205775e+01,
   },
   segment1 = {
      1.034972096e+06,
     -2.412698562e+03,
      4.646110780e+00,
      2.291998307e-03,
     -6.836830480e-07,
      9.426468930e-11,
     -4.822380530e-15,
     -1.384286509e+04,
     -7.978148510e+00,
   },
}
db.H2O.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  5.0019557e-01,
      B = -6.9712796e+02,
      C =  8.8163892e+04,
      D =  3.0836508e+00
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  5.8988538e-01,
      B = -5.3769814e+02,
      C =  5.4263513e+04,
      D =  2.3386375e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4330087e-01,
      B = -9.5668913e+01,
      C = -3.7742283e+05,
      D =  1.8125190e+00
   },
}
db.H2O.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  1.0966389e+00,
      B = -5.5513429e+02,
      C =  1.0623408e+05,
      D = -2.4664550e-01
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  3.9367933e-01,
      B = -2.2524226e+03,
      C =  6.1217458e+05,
      D =  5.8011317e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.1858737e-01,
      B = -1.4096649e+04,
      C =  1.9179190e+07,
      D =  1.4345613e+01
   },
}
db.H2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.19864056E+00,
     -2.03643410E-03,
      6.52040211E-06,
     -5.48797062E-09,
      1.77197817E-12,
     -3.02937267E+04,
     -8.49032208E-01,
   },
   segment1 = {
      0,
      0,
      3.03399249E+00,
      2.17691804E-03,
     -1.64072518E-07,
     -9.70419870E-11,
      1.68200992E-14,
     -3.00042971E+04,
      4.96677010E+00,
   }
}

db.H2O.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.596254204075e+01,
      B = 1.678935618571e-01,
      C = 1.929065548334e-01,
      D = -1.356763840380e-02,
   }
}
db.H2O.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = 5.416062110310e+00,
      B = -5.968633913268e+00,
      C = 1.091216246422e+00,
      D = -5.536357184347e-02,
   }
}

db.H2O.Hf = {
   value = -241826.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
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
db.HCCO = {}
db.HCCO.atomicConstituents = {C=2,H=1,O=1,}
db.HCCO.charge = 0
db.HCCO.M = {
   value = 41.028740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCCO.gamma = {
   value = 1.2067e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCCO.sigma = {
   value = 2.500,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCO.epsilon = {
   value = 150.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 4000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.22517214E+01,
      0.17655021E-01,
     -0.23729101E-04,
      0.17275759E-07,
     -0.50664811E-11,
      0.20059449E+05,
      0.12490417E+02,
   },
   segment1 = {
      0,
      0,
      0.56282058E+01,
      0.40853401E-02,
     -0.15934547E-05,
      0.28626052E-09,
     -0.19407832E-13,
      0.19327215E+05,
     -0.39302595E+01,
   }
}
db.HCCO.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      6.959613e+04,
     -1.164594e+03,
      9.456616e+00,
     -2.331241e-03,
      5.161874e-06,
     -3.526170e-09,
      8.599143e-13,
      2.535004e+04,
     -2.726355e+01,
   },
   segment1 = {
      1.093922e+06,
     -4.498228e+03,
      1.246446e+01,
     -6.343317e-04,
      1.108549e-07,
     -1.125489e-11,
      5.689152e-16,
      4.652280e+04,
     -5.099070e+01
    }
}

db.HCCO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 4000.000,
      A = -2.047911174161e+01,
      B = 3.289579534067e+00,
      C = -3.425514016760e-01,
      D = 1.478685784020e-02,
   }
}
db.HCCO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 4000.000,
      A = -1.709974427710e+01,
      B = 4.224884718432e+00,
      C = -3.726422242450e-01,
      D = 1.275222078039e-02,
   }
}

db.HCCO.Hf = {
   value = 176568.1,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HCCOH = {}
db.HCCOH.atomicConstituents = {C=2,H=2,O=1,}
db.HCCOH.charge = 0
db.HCCOH.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCCOH.gamma = {
   value = 1.1656e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCCOH.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCOH.epsilon = {
   value = 436.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCCOH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.12423733E+01,
      0.31072201E-01,
     -0.50866864E-04,
      0.43137131E-07,
     -0.14014594E-10,
      0.80316143E+04,
      0.13874319E+02,
   },
   segment1 = {
      0,
      0,
      0.59238291E+01,
      0.67923600E-02,
     -0.25658564E-05,
      0.44987841E-09,
     -0.29940101E-13,
      0.72646260E+04,
     -0.76017742E+01,
   }
}
db.HCN = {}
db.HCN.atomicConstituents = {C=1,H=1,N=1,}
db.HCN.charge = 0
db.HCN.M = {
   value = 27.0253400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HCN.gamma = {
   value = 1.301,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HCN.sigma = {
   value = 3.63,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCN.epsilon = {
   value = 569.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.22589886E+01,
      0.10051170E-01,
     -0.13351763E-04,
      0.10092349E-07,
     -0.30089028E-11,
      0.14712633E+05,
      0.89164419E+01,
   },
   segment1 = {
      0,
      0,
      0.38022392E+01,
      0.31464228E-02,
     -0.10632185E-05,
      0.16619757E-09,
     -0.97997570E-14,
      0.14407292E+05,
      0.15754601E+01,
   }
}
db.HCN.Hf = {
   value = 133082.46,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HCNN = {}
db.HCNN.atomicConstituents = {C=1,H=1,N=2,}
db.HCNN.charge = 0
db.HCNN.M = {
   value = 41.032040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCNN.gamma = {
   value = 1.2032e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCNN.sigma = {
   value = 2.50,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCNN.epsilon = {
   value = 150.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCNN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.25243194E+01,
      0.15960619E-01,
     -0.18816354E-04,
      0.12125540E-07,
     -0.32357378E-11,
      0.54261984E+05,
      0.11675870E+02,
   },
   segment1 = {
      0,
      0,
      0.58946362E+01,
      0.39895959E-02,
     -0.15982380E-05,
      0.29249395E-09,
     -0.20094686E-13,
      0.53452941E+05,
     -0.51030502E+01,
   }
}
db.HCNO = {}
db.HCNO.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HCNO.charge = 0
db.HCNO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HCNO.gamma = {
   value = 1.2154e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HCNO.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCNO.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCNO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1382.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.64727989E+00,
      1.27505342E-02,
     -1.04794236E-05,
      4.41432836E-09,
     -7.57521466E-13,
      1.92990252E+04,
      1.07332972E+01,
   },
   segment1 = {
      0,
      0,
      6.59860456E+00,
      3.02778626E-03,
     -1.07704346E-06,
      1.71666528E-10,
     -1.01439391E-14,
      1.79661339E+04,
     -1.03306599E+01,
   }
}
db.HCO = {}
db.HCO.atomicConstituents = {C=1,H=1,O=1,}
db.HCO.charge = 0
db.HCO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp',
}
db.HCO.gamma = {
   value = 1.316,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = "Gokel (2004), Dean's Handbook of Organic Chemistry",
}
db.HCO.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCO.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.22118584E+00,
     -3.24392532E-03,
      1.37799446E-05,
     -1.33144093E-08,
      4.33768865E-12,
      3.83956496E+03,
      3.39437243E+00,
   },
   segment1 = {
      0,
      0,
      2.77217438E+00,
      4.95695526E-03,
     -2.48445613E-06,
      5.89161778E-10,
     -5.33508711E-14,
      4.01191815E+03,
      9.79834492E+00,
   }
}
db.HCO.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.189852e+04,
      2.151536e+02,
      2.730224e+00,
      1.806516e-03,
      4.984301e-06,
     -5.814568e-09,
      1.869690e-12,
      2.905756e+03,
      1.136773e+01,
   },
   segment1 = {
      6.949606e+05,
     -3.656223e+03,
      9.604731e+00,
     -1.117129e-03,
      2.875328e-07,
     -3.626248e-11,
      1.808330e-15,
      2.543704e+04,
     -3.582474e+01
    }
}

db.HCO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.946125270029e+01,
      B = 1.583449384081e+00,
      C = -1.418736403518e-02,
      D = -3.703072867374e-03,
   }
}
db.HCO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.790317887628e+00,
      B = -3.149564742925e+00,
      C = 7.513430765797e-01,
      D = -4.266290638275e-02,
   }
}

db.HCO.Hf = {
   value = 42397.85,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HI = {}
db.HI.atomicConstituents = {I=1,H=1,}
db.HI.charge = 0
db.HI.M = {
   value = 0.12791241,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.HI.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.HI.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      1.872881730e+04,
     -3.431788840e+02,
      5.956712430e+00,
     -8.543439600e-03,
      1.454780274e-05,
     -1.049104164e-08,
      2.839734003e-12,
      3.682950720e+03,
     -8.149756090e+00,
   },
   segment1 = {
      4.724921450e+05,
     -1.923465741e+03,
      5.758048970e+00,
     -4.066266380e-04,
      9.474332050e-08,
     -1.033534431e-11,
      4.611614790e-16,
      1.394857037e+04,
     -1.182487652e+01,
   },
}
db.HI.Hf = {
   value = 26359.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HNCO = {}
db.HNCO.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HNCO.charge = 0
db.HNCO.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HNCO.gamma = {
   value = 1.2173e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HNCO.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNCO.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {300.0, 1478.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.63096317E+00,
      7.30282357E-03,
     -2.28050003E-06,
     -6.61271298E-10,
      3.62235752E-13,
     -1.55873636E+04,
      6.19457727E+00,
   },
   segment1 = {
      0,
      0,
      6.22395134E+00,
      3.17864004E-03,
     -1.09378755E-06,
      1.70735163E-10,
     -9.95021955E-15,
     -1.66599344E+04,
     -8.38224741E+00,
   }
}
db.HNCO.Hf = {
   value = -118056.529,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HNO = {}
db.HNO.atomicConstituents = {H=1,N=1,O=1}
db.HNO.M = {
   value = 31.01404e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO.charge = 0
db.HNO.gamma = {
   value = 1.325,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO.sigma = {
   value = 3.492,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNO.epsilon = {
   value = 116.7,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HNO.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -6.854764860e+04,
      9.551627200e+02,
     -6.000720210e-01,
      7.995176750e-03,
     -6.547079160e-07,
     -3.670513400e-09,
      1.783392519e-12,
      6.435351260e+03,
      3.048166179e+01,
   },
   segment1 = {
     -5.795614980e+06,
      1.945457427e+04,
     -2.152568374e+01,
      1.797428992e-02,
     -4.976040670e-06,
      6.397924170e-10,
     -3.142619368e-14,
     -1.104192372e+05,
      1.818650338e+02,
   }
}
db.HNO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.45334916E+01,
     -0.56696171E-02,
      0.18473207E-04,
     -0.17137094E-07,
      0.55454573E-11,
      0.11548297E+05,
      0.17498417E+01,
   },
   segment1 = {
      0,
      0,
      0.29792509E+01,
      0.34944059E-02,
     -0.78549778E-06,
      0.57479594E-10,
     -0.19335916E-15,
      0.11750582E+05,
      0.86063728E+01,
   }
}
db.HNO.Hf = {
   value = 102032.725,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HNO2 = {}
db.HNO2.atomicConstituents = {H=1,N=1,O=2}
db.HNO2.charge = 0
db.HNO2.M = {
   value = 0.04701344,
   unit = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO2.gamma = {
   value = 1.218,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      8.591985060e+03,
      1.203644046e+02,
      9.412979120e-01,
      1.942891839e-02,
     -2.253174194e-05,
      1.384587594e-08,
     -3.473550460e-12,
     -1.106337202e+04,
      2.073967331e+01,
   },
   segment1 = {
      8.787904130e+05,
     -3.990455030e+03,
      1.187349269e+01,
     -4.881900610e-04,
      7.133636790e-08,
     -5.376303340e-12,
      1.581778986e-16,
      1.246343241e+04,
     -4.608874688e+01,
   }
}
db.HNO2.Hf = {
   value = -78451.922,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HNO3 = {}
db.HNO3.atomicConstituents = {H=1,N=1,O=3}
db.HNO3.M = {
   value = 63.01284e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HNO3.charge = 0
db.HNO3.gamma = {
   value = 1.181,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HNO3.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      9.202869010e+03,
      1.093774496e+02,
     -4.521042450e-01,
      2.984914503e-02,
     -3.190635500e-05,
      1.720931528e-08,
     -3.782649830e-12,
     -1.764048507e+04,
      2.746644879e+01,
   },
   segment1 = {
     -9.497809640e+04,
     -2.733024468e+03,
      1.449426995e+01,
     -7.821868050e-04,
      1.702693665e-07,
     -1.930543961e-11,
      8.870455120e-16,
     -4.882517780e+03,
     -5.928392985e+01,
   }
}
db.HNO3.Hf = {
   value = -133912.869,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HO2 = {}
db.HO2.atomicConstituents = {O=2,H=1,}
db.HO2.charge = 0
db.HO2.M = {
   value = 0.03300674,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.HO2.gamma = {
   value = 1.31200000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.HO2.sigma = {
   value = 3.458,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HO2.epsilon = {
   value = 107.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HO2.Lewis = {
   value = 1.079
}
db.HO2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -7.598882540e+04,
      1.329383918e+03,
     -4.677388240e+00,
      2.508308202e-02,
     -3.006551588e-05,
      1.895600056e-08,
     -4.828567390e-12,
     -5.873350960e+03,
      5.193602140e+01,
   },
   segment1 = {
     -1.810669724e+06,
      4.963192030e+03,
     -1.039498992e+00,
      4.560148530e-03,
     -1.061859447e-06,
      1.144567878e-10,
     -4.763064160e-15,
     -3.200817190e+04,
      4.066850920e+01,
   },
}
db.HO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      4.30179801E+00,
     -4.74912051E-03,
      2.11582891E-05,
     -2.42763894E-08,
      9.29225124E-12,
      2.94808040E+02,
      3.71666245E+00,
   },
   segment1 = {
      0,
      0,
      4.01721090E+00,
      2.23982013E-03,
     -6.33658150E-07,
      1.14246370E-10,
     -1.07908535E-14,
      1.11856713E+02,
      3.78510215E+00,
   }
}

db.HO2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.998505987333e+01,
      B = 2.898260095296e+00,
      C = -3.018510665245e-01,
      D = 1.346539991802e-02,
   }
}
db.HO2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -9.026815991823e+00,
      B = 7.723598729587e-01,
      C = 6.217081555068e-02,
      D = -4.947480321696e-03,
   }
}

db.HO2.Hf = {
   value = 12020.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.HOCN = {}
db.HOCN.atomicConstituents = {C=1,H=1,N=1,O=1,}
db.HOCN.charge = 0
db.HOCN.M = {
   value = 43.024740e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.HOCN.gamma = {
   value = 1.2185e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.HOCN.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HOCN.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HOCN.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1368.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.78604952E+00,
      6.88667922E-03,
     -3.21487864E-06,
      5.17195767E-10,
      1.19360788E-14,
     -2.82698400E+03,
      5.63292162E+00,
   },
   segment1 = {
      0,
      0,
      5.89784885E+00,
      3.16789393E-03,
     -1.11801064E-06,
      1.77243144E-10,
     -1.04339177E-14,
     -3.70653331E+03,
     -6.18167825E+00,
   }
}
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
-- Helium_plus ported from Daniel's He_plus.lua file in the cfcfd3 collection
-- Yu Liu, 2018-06-27
db['He+'] = {}
db['He+'].atomicConstituents = {He=1,}
db['He+'].charge = 1
db['He+'].M = {
   value = 4.0020534e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db['He+'].gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['He+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      0.000000000e+00,  
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,  
      0.000000000e+00,
      0.000000000e+00, 
      2.853233739e+05, 
      1.621665557e+00
   },
   segment1 = {
      0.000000000e+00,  
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,  
      0.000000000e+00,
      0.000000000e+00, 
      2.853233739e+05, 
      1.621665557e+00
   },
   segment2 = {
      0.000000000e+00,  
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,  
      0.000000000e+00,
      0.000000000e+00, 
      2.853233739e+05, 
      1.621665557e+00
   },
   reference = 'cea2::thermo.inp'
}
-- No CEA transport data for He+, just use He
db['He+'].ceaViscosity = db.He.ceaViscosity 
db['He+'].ceaThermCond = db.He.ceaThermCond 
db['He+'].Hf = {
   value = 2378521.473,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
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
-- Neon_plus 
-- Yu Liu, 2018-06-16
db['Ne+'] = {}
db['Ne+'].atomicConstituents = {Ne=1,}
db['Ne+'].charge = 1
db['Ne+'].M = {
   value = 20.1797000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db['Ne+'].gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['Ne+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      7.281551480E+04,
     -8.695697990E+02, 
      6.108646970E+00,
     -5.841356930E-03, 
      5.041044170E-06,
     -2.293759207E-09, 
      4.339065680E-13,                 
      2.545996890E+05,
     -1.673449355E+01
   },
   segment1 = {
     -1.112742658E+05, 
      4.765697970E+02, 
      2.196650531E+00, 
      1.102593151E-04,
     -2.287564425E-08,
      2.510218183E-12,
     -1.126646096E-16,                
      2.472536944E+05, 
      7.466140540E+00
   },
   segment2 = {
     -5.615474110E+04, 
      1.418980160E+02,
      2.475716842E+00, 
      1.944430992E-06,
     -6.323099200E-11,
     -1.313313446E-16, 
      3.534699010E-20,               
      2.494452217E+05,
      5.366882220E+00
   },
   reference = 'cea2::thermo.inp'
}
-- No CEA transport data for Ne+, just use Ne
db['Ne+'].ceaViscosity = db.Ne.ceaViscosity 
db['Ne+'].ceaThermCond = db.Ne.ceaThermCond 
db['Ne+'].Hf = {
   value = 2086965.946,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.I2 = {}
db.I2.atomicConstituents = {I=2,}
db.I2.charge = 0
db.I2.M = {
   value = 0.25380894,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
db.I2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.I2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -5.087968770e+03,
     -1.249585210e+01,
      4.504219090e+00,
      1.370962533e-04,
     -1.390523014e-07,
      1.174813853e-10,
     -2.337541043e-14,
      6.213469810e+03,
      5.583836940e+00,
   },
   segment1 = {
     -5.632594160e+06,
      1.793961560e+04,
     -1.723055169e+01,
      1.244214080e-02,
     -3.332768580e-06,
      4.125477940e-10,
     -1.960461713e-14,
     -1.068505292e+05,
      1.600531883e+02,
   },
}
db.I2.Hf = {
   value = 62420.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.Kr = {}
db.Kr.atomicConstituents = {Kr=1,}
db.Kr.charge = 0
db.Kr.M = {
   value = 0.0838,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db.Kr.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.Kr.ceaThermoCoeffs = {
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
      5.490956510e+00
   },
   segment1 = {
      2.643639057e+02,
     -7.910050820e-01,
      2.500920585e+00,
     -5.328164110e-07,
      1.620730161e-10,
     -2.467898017e-14,
      1.478585040e-18,
     -7.403488940e+02,
      5.484398150e+00
   },
   segment2 = {
     -1.375531087e+09,
      9.064030530e+05,
     -2.403481435e+02,
      3.378312030e-02,
     -2.563103877e-06,
      9.969787790e-11,
     -1.521249677e-15,
     -7.111667370e+06,
      2.086866326e+03
   },
}
db.Kr.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58597795e+00,
      B = -0.12924832e+03,
      C =  0.47495759e+04,
      D =  0.25793650e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.68985845e+00,
      B =  0.56296306e+02,
      C = -0.36082600e+05,
      D =  0.17170715e+01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  0.76582939e+00,
      B =  0.68610377e+03,
      C = -0.82537190e+06,
      D =  0.97565128e+00
   },
}
db.Kr.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  0.58008139e+00,
      B = -0.13792556e+03,
      C =  0.60771460e+04,
      D =  0.16420039e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.68859431e+00,
      B =  0.51765647e+02,
      C = -0.34512131e+05,
      D =  0.74332130e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  0.76365443e+00,
      B =  0.65175847e+03,
      C = -0.73589800e+06,
      D =  0.12112126e-01
   },
}
db.Kr.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.N = {}
db.N.type = "atom"
db.N.atomicConstituents = {N=1}
db.N.charge = 0
db.N.M = {
   value = 14.0067e-3,
   units = 'kg/mol',
   reference = 'CEA2::thermo.inp'
}
db.N.gamma =  {
   value = 5./3.,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.N.sigma = {
   value = 3.298,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Svehla (1962)'
}
db.N.epsilon = {
   value = 71.4, 
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Svehla (1962)'
}
db.N.electronic_levels = {
  Te = {
    value = {0.0, 19227.95, 28839.18, 83335.6, 86192.79, 88132.45, 93581.55},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {4, 10, 6, 12, 6, 12, 2},
    units = 'NA',
    description = 'degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db.N.ionisation_energy = {
  value = 14.53413,
  units = "eV",
  description = "Ionisation energy",
  reference = "https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html"
}
db.N.ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      5.61061063e+04,
      4.19425139e+00
   },
   segment1 = {
     -2.27073277e+05,
      8.14052944e+02,
      1.32705137e+00,
      8.62721731e-04,
     -3.35747089e-07,
      6.29010687e-11,
     -3.90674587e-15,
      5.10943141e+04,
      1.22823713e+01
   },
   segment2 = {
     -2.04738994e+09,
      1.45842847e+06,
     -4.18833824e+02,
      6.25994407e-02,
     -4.96542822e-06,
      1.98247052e-10,
     -3.05470194e-15,
     -1.12727730e+07,
      3.58487417e+03
   },
   segment3 = {
      5.74291902e+11,
     -1.29039294e+08,
      1.15381467e+04,
     -5.25078568e-01,
      1.29219090e-05,
     -1.63974231e-10,
      8.41878585e-16,
      1.15261836e+09,
     -1.11649232e+05
   },
   reference="Johnston et al. AIAA Paper 2015-3110"
}
db.N.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.25000000E+01,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
      0.00000000E+00,
      0.56104637E+05,
      0.41939087E+01,
   },
   segment1 = {
      0,
      0,
      0.24159429E+01,
      0.17489065E-03,
     -0.11902369E-06,
      0.30226245E-10,
     -0.20360982E-14,
      0.56133773E+05,
      0.46496096E+01,
   }
}
db.N.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83724737e+00,
      B = 0.43997150e+03,
      C = -0.17450753e+06,
      D = 0.10365689e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.89986588e+00,
      B = 0.14112801e+04,
      C = -0.18200478e+07,
      D = -0.55811716e+00
   },
   reference = 'from CEA2::trans.inp which cites Levin et al (1990)'
}
db.N.ceaThermCond = {
   nsegments = 2, 
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A = 0.83771661e+00,
      B = 0.44243270e+03,
      C = -0.17578446e+06,
      D = 0.89942915e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = 0.90001710e+00,
      B = 0.14141175e+04,
      C = -0.18262403e+07,
      D = 0.24048513e+00
   },
   reference = 'from CEA2::trans.inp which cites Levin et al (1990)'
}

db.N.Hf = {
   value = 472680.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['N+'] = {}
db['N+'].atomicConstituents = {N=1}
db['N+'].type = "atom"
db['N+'].charge = 1
db['N+'].M = {
   value = 14.0061514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['N+'].gamma =  {
   value = 1.641, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
db['N+'].electronic_levels = {
  Te = {
    value = {0.0, 70.07, 188.19, 22036.47, 47031.63, 67312.23},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {1, 3, 5, 5, 1, 5},
    units = 'NA',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db['N+'].ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      5.23707921e+03,
      2.29995832e+00,
      2.48237856e+00,
      2.73749076e-05,
     -3.13444758e-08,
      1.85011133e-11,
     -4.44735098e-15,
      2.25626950e+05,
      5.04771460e+00
   },
   segment1 = {
      3.51764922e+05,
     -1.04468568e+03,
      3.68498994e+00,
     -6.41320819e-04,
      1.66764737e-07,
     -1.82509101e-11,
      7.37893962e-16,
      2.32355850e+05,
     -3.48540450e+00
   },
   segment2 = {
      2.08931831e+07,
     -1.36213868e+04,
      5.54565345e+00,
     -2.68374643e-04,
      1.46878484e-08,
     -4.22829994e-13,
      5.70569781e-18,
      3.33626095e+05,
     -2.20294391e+01
   },
   segment3 = {
      1.48351337e+09,
     -2.69100020e+05,
      1.91340133e+01,
     -2.13380850e-04,
     -1.51589694e-08,
      5.96512853e-13,
     -5.46224064e-18,
      2.70008612e+06,
     -1.64949269e+02
   },
   reference="Johnston et al. AIAA Paper 2015-3110"
}
-- No CEA transport data for N+, just use N
db['N+'].ceaViscosity = db.N.ceaViscosity 
db['N+'].ceaThermCond = db.N.ceaThermCond 


db['N+'].Hf = {
   value = 1882127.624,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.N2 = {}
db.N2.type = "molecule"
db.N2.molecule_type = "linear"
db.N2.theta_v = {
   value = 3393.44,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Capitelli et al (2005), Table 15. omega_e in ground state converted to K'
}
db.N2.theta_D = {
  value = 78740.00,
  units = 'K',
  description = 'characteristic dissociation tempereature',
  reference = '???'
}
db.N2.electronic_levels = {
  g  = {
   value = {1, 3, 6, 6, 3, 1, 2},
   units = 'NA',
   description = 'degeneracy of electronic levels',
   reference = 'NIST'
  },
  Te = {
    value = {0.0, 50203.66, 59619.09, 59808.00, 66272.5, 68152.66, 69283.06},
    units = 'cm^(-1)',
    description = 'Electronic excitation energy',
    reference = 'NIST',
  }
}
db.N2.ionisation_energy = {
  value = 15.581,
  units = 'eV',
  descritpion = 'Ionisation energy',
  reference = 'https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=20'
}
db.N2.atomicConstituents = {N=2}
db.N2.charge = 0
db.N2.M = {
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.N2.gamma = {
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.N2.sigma = {
   value = 3.621,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2.epsilon = {
   value = 97.530,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2.Lewis = {
   value = 1.152
}
db.N2.entropyRefValues = {
   s1 = 6839.91,
   T1 = 298.15,
   p1 = 101.325e3,
   description = 'Standard state entropy at 1 bar',
   reference = 'NIST Chemistry WebBook: http://webbook.nist.gov/chemistry/'
}
db.N2.sutherlandVisc = {
   mu_ref = 1.663e-05,
   T_ref = 273.0,
   S = 107.0,
   reference = "Table 1-2, White (2006)"
}
db.N2.sutherlandThermCond = {
   k_ref = 0.0242,
   T_ref = 273.0,
   S = 150.0,
   reference = "Table 1-3, White (2006)"
}
db.N2.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      2.210371497e+04,
     -3.818461820e+02,
      6.082738360e+00,
     -8.530914410e-03,
      1.384646189e-05,
     -9.625793620e-09,
      2.519705809e-12,
      7.108460860e+02,
     -1.076003744e+01
   },
   segment1 = {
      5.877124060e+05,
     -2.239249073e+03,
      6.066949220e+00,
     -6.139685500e-04,
      1.491806679e-07,
     -1.923105485e-11,
      1.061954386e-15,
      1.283210415e+04,
     -1.586640027e+01
   },
   segment2 = {
      8.310139160e+08,
     -6.420733540e+05,
      2.020264635e+02,
     -3.065092046e-02,
      2.486903333e-06,
     -9.705954110e-11,
      1.437538881e-15,
      4.938707040e+06,
     -1.672099740e+03
   },
   ref="from CEA2::thermo.inp"
}
db.N2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.03298677E+02,
      0.14082404E-02,
     -0.03963222E-04,
      0.05641515E-07,
     -0.02444854E-10,
     -0.10208999E+04,
      0.03950372E+02,
   },
   segment1 = {
      0,
      0,
      0.02926640E+02,
      0.14879768E-02,
     -0.05684760E-05,
      0.10097038E-09,
     -0.06753351E-13,
     -0.09227977E+04,
      0.05980528E+02,
   }
}
db.N2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200,
      T_upper = 1000,
      A = 0.62526577,
      B = -31.779652,
      C = -1640.7983,
      D = 1.7454992,
    },
   segment1 = {
      T_lower = 1000,
      T_upper = 5000,
      A = 0.87395209,
      B = 561.52222,
      C = -173948.09,
      D = -0.39335958,
    },
   segment2 = {
      T_lower = 5000,
      T_upper = 15000,
      A = 0.88503551,
      B = 909.02171,
      C = -731290.61,
      D = -0.53503838,
    },
}
db.N2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200,
      T_upper = 1000,
      A = 0.85439436,
      B = 105.73224,
      C = -12347.848,
      D = 0.47793128,
    },
    segment1 =  {
      T_lower = 1000,
      T_upper = 5000,
      A = 0.88407146,
      B = 133.57293,
      C = -11429.64,
      D = 0.24417019,
    },
    segment2 = {
      T_lower = 5000,
      T_upper = 15000,
      A = 2.4176185,
      B = 8047.7749,
      C = 3105580.2,
      D = -14.517761,
    },
}

db.N2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -1.756744424273e+01,
      B = 1.819353013613e+00,
      C = -1.517111911440e-01,
      D = 6.520135228694e-03,
   }
}
db.N2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -4.129764724745e+00,
      B = -1.096196318523e+00,
      C = 2.922595492645e-01,
      D = -1.493636112632e-02,
   }
}

db.N2.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['N2+'] = {}
db['N2+'].type = "molecule"
db['N2+'].molecule_type = "linear"
db['N2+'].theta_v = {
   value = 3393.44,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Copied from N2'
}
db['N2+'].atomicConstituents = {N=2}
db['N2+'].charge = 1
db['N2+'].M = {
   value = 28.0128514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db['N2+'].gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db['N2+'].electronic_levels = {
  Te = {
    value = {0.0, 9167.46, 25461.11, 51663.2, 64609.03},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {2, 4, 2, 4, 2},
    units = 'cm^(-1)',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST',
  }
}
db['N2+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -3.474047470e+04,
      2.696222703e+02,
      3.164916370e+00,
     -2.132239781e-03,
      6.730476400e-06,
     -5.637304970e-09,
      1.621756000e-12,
      1.790004424e+05,
      6.832974166e+00
   },
   segment1 = {
     -2.845599002e+06,
      7.058893030e+03,
     -2.884886385e+00,
      3.068677059e-03,
     -4.361652310e-07,
      2.102514545e-11, 
      5.411996470e-16,
      1.340388483e+05,
      5.090897022e+01
   },
   segment2 = {
     -3.712829770e+08,
      3.139287234e+05,
     -9.603518050e+01, 
      1.571193286e-02,
     -1.175065525e-06,
      4.144441230e-11,
     -5.621893090e-16,
     -2.217361867e+06,
      8.436270947e+02
   },
   ref="from CEA2::thermo.inp"
}
-- No CEA transport data for N2+, just use N2
db['N2+'].ceaViscosity = db.N2.ceaViscosity
db['N2+'].ceaThermCond = db.N2.ceaThermCond



db['N2+'].Hf = {
   value = 1509508.424,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.N2O = {}
db.N2O.atomicConstituents = {N=2,O=1,}
db.N2O.charge = 0
db.N2O.M = {
   value = 44.012800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.N2O.gamma = {
   value = 1.2735e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.N2O.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2O.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.N2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.22571502E+01,
      0.11304728E-01,
     -0.13671319E-04,
      0.96819806E-08,
     -0.29307182E-11,
      0.87417744E+04,
      0.10757992E+02,
   },
   segment1 = {
      0,
      0,
      0.48230729E+01,
      0.26270251E-02,
     -0.95850874E-06,
      0.16000712E-09,
     -0.97752303E-14,
      0.80734048E+04,
     -0.22017207E+01,
   }
}
db.N2O.Hf = {
   value = 81600.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.NCO = {}
db.NCO.atomicConstituents = {C=1,N=1,O=1,}
db.NCO.charge = 0
db.NCO.M = {
   value = 42.016800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NCO.gamma = {
   value = 1.2609e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NCO.sigma = {
   value = 3.828,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NCO.epsilon = {
   value = 232.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.28269308E+01,
      0.88051688E-02,
     -0.83866134E-05,
      0.48016964E-08,
     -0.13313595E-11,
      0.14682477E+05,
      0.95504646E+01,
   },
   segment1 = {
      0,
      0,
      0.51521845E+01,
      0.23051761E-02,
     -0.88033153E-06,
      0.14789098E-09,
     -0.90977996E-14,
      0.14004123E+05,
     -0.25442660E+01,
   }
}
db.NCO.Hf = {
   value = 131847.241,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.NH = {}
db.NH.atomicConstituents = {N=1,H=1,}
db.NH.charge = 0
db.NH.M = {
   value = 15.0146400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NH.gamma = {
   value = 1.398,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NH.sigma = {
   value = 2.65,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH.epsilon = {
   value = 80.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.34929085E+01,
      0.31179198E-03,
     -0.14890484E-05,
      0.24816442E-08,
     -0.10356967E-11,
      0.41880629E+05,
      0.18483278E+01,
   },
   segment1 = {
      0,
      0,
      0.27836928E+01,
      0.13298430E-02,
     -0.42478047E-06,
      0.78348501E-10,
     -0.55044470E-14,
      0.42120848E+05,
      0.57407799E+01,
   }
}
db.NH.Hf = {
   value = 357032.001,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.NH2 = {}
db.NH2.atomicConstituents = {N=1,H=2,}
db.NH2.charge = 0
db.NH2.M = {
   value = 16.022580e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH2.gamma = {
   value = 1.3254e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH2.sigma = {
   value = 2.65,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH2.epsilon = {
   value = 80.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.42040029E+01,
     -0.21061385E-02,
      0.71068348E-05,
     -0.56115197E-08,
      0.16440717E-11,
      0.21885910E+05,
     -0.14184248E+00,
   },
   segment1 = {
      0,
      0,
      0.28347421E+01,
      0.32073082E-02,
     -0.93390804E-06,
      0.13702953E-09,
     -0.79206144E-14,
      0.22171957E+05,
      0.65204163E+01,
   }
}
db.NH2.Hf = {
   value = 189134.713,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.NH3 = {}
db.NH3.atomicConstituents = {N=1,H=3,}
db.NH3.charge = 0
db.NH3.M = {
   value = 17.030520e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NH3.gamma = {
   value = 1.3036e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NH3.sigma = {
   value = 2.92,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH3.epsilon = {
   value = 481.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.42860274E+01,
     -0.46605230E-02,
      0.21718513E-04,
     -0.22808887E-07,
      0.82638046E-11,
     -0.67417285E+04,
     -0.62537277E+00,
   },
   segment1 = {
      0,
      0,
      0.26344521E+01,
      0.56662560E-02,
     -0.17278676E-05,
      0.23867161E-09,
     -0.12578786E-13,
     -0.65446958E+04,
      0.65662928E+01,
   }
}
db.NH3.Hf = {
   value = -45940.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.NNH = {}
db.NNH.atomicConstituents = {N=2,H=1,}
db.NNH.charge = 0
db.NNH.M = {
   value = 29.021340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.NNH.gamma = {
   value = 1.3152e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.NNH.sigma = {
   value = 3.798,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NNH.epsilon = {
   value = 71.4,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NNH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.43446927E+01,
     -0.48497072E-02,
      0.20059459E-04,
     -0.21726464E-07,
      0.79469539E-11,
      0.28791973E+05,
      0.29779410E+01,
   },
   segment1 = {
      0,
      0,
      0.37667544E+01,
      0.28915082E-02,
     -0.10416620E-05,
      0.16842594E-09,
     -0.10091896E-13,
      0.28650697E+05,
      0.44705067E+01,
   }
}
db.NO = {}
db.NO.type = "molecule"
db.NO.molecule_type = "linear"
db.NO.theta_v = {
   value = 2739.70, 
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Capitelli et al (2005), Table 19. omega_e in ground state converted to K'
}
db.NO.electronic_levels = {
  Te = {
    --value = {60.55, 38440.00, 43965.7, 45932.3, 48680.0, 52175.7},
    value = {0.0, 38440.00, 43965.7, 45932.3, 48680.0, 52175.7},
    --value = {0.0, 120.94, 38807, 43966, 45932,47950,52186,53085,53637},
    units = 'cm^(-1)',
    description = 'Electronic energies',
    reference = 'NIST'
  },
  g = {
    value = {4, 8, 2, 4, 4, 4},
    --value={2,2,8,2,4,4,4,2,4},
    units = 'NA',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST',
  }
}
db.NO.atomicConstituents = {O=1,N=1,}
db.NO.charge = 0
db.NO.M = {
   value = 0.03000610,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NO.gamma = {
   value = 1.38600000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NO.sigma = {
   value = 3.621,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NO.epsilon = {
   value = 97.53,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NO.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -1.143916503e+04,
      1.536467592e+02,
      3.431468730e+00,
     -2.668592368e-03,
      8.481399120e-06,
     -7.685111050e-09,
      2.386797655e-12,
      9.098214410e+03,
      6.728725490e+00,
   },
   segment1 = {
      2.239018716e+05,
     -1.289651623e+03,
      5.433936030e+00,
     -3.656034900e-04,
      9.880966450e-08,
     -1.416076856e-11,
      9.380184620e-16,
      1.750317656e+04,
     -8.501669090e+00,
   },
   segment2 = {
     -9.575303540e+08,
      5.912434480e+05,
     -1.384566826e+02,
      1.694339403e-02,
     -1.007351096e-06,
      2.912584076e-11,
     -3.295109350e-16,
     -4.677501240e+06,
      1.242081216e+03,
   },
}
db.NO.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0262029e-01,
      B = -6.2017783e+01,
      C = -1.3954524e+02,
      D =  2.0268332e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.8009050e-01,
      B =  3.0486891e+02,
      C = -9.4847722e+04,
      D =  5.2873381e-01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.0580582e-01,
      B =  6.2427878e+02,
      C = -5.7879210e+05,
      D =  2.6516450e-01
   },
}
db.NO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.42184763E+01,
     -0.46389760E-02,
      0.11041022E-04,
     -0.93361354E-08,
      0.28035770E-11,
      0.98446230E+04,
      0.22808464E+01,
   },
   segment1 = {
      0,
      0,
      0.32606056E+01,
      0.11911043E-02,
     -0.42917048E-06,
      0.69457669E-10,
     -0.40336099E-14,
      0.99209746E+04,
      0.63693027E+01,
   }
}
db.NO.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  9.5028758e-01,
      B =  7.6667058e+01,
      C = -9.9894764e+03,
      D = -6.2776717e-03
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.6215238e-01,
      B =  4.4568223e+02,
      C = -2.3856466e+05,
      D =  4.6209876e-01
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.0377865e+00,
      B = -3.4486864e+04,
      C =  6.7451187e+07,
      D =  2.0928749e+01
   },
}
db.NO.Hf = {
   value = 91271.31,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['NO+'] = {}
db['NO+'].type = "molecule"
db['NO+'].molecule_type = "linear"
db['NO+'].theta_v = {
   value = 2739.70, 
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Copied from NO'
}
db['NO+'].electronic_levels = {
  Te = {
    value={0.0, 52190.00, 59240.00, 61880.0, 67720.0, 71450.0, 73471.72},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {1, 3, 6, 6, 3, 1, 2, 2},
    units = 'cm^(-1)',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db['NO+'].atomicConstituents = {O=1,N=1,}
db['NO+'].charge = 1
db['NO+'].M = {
   value = 30.0055514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['NO+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      1.398106635e+03,
     -1.590446941e+02,
      5.122895400e+00,
     -6.394388620e-03,
      1.123918342e-05,
     -7.988581260e-09,
      2.107383677e-12,
      1.187495132e+05,
     -4.398433810e+00
   },
   segment1 = {
      6.069876900e+05,
     -2.278395427e+03,
      6.080324670e+00,
     -6.066847580e-04,
      1.432002611e-07,
     -1.747990522e-11,
      8.935014060e-16,
      1.322709615e+05,
     -1.519880037e+01
   },
   segment2 = {
      2.676400347e+09,
     -1.832948690e+06,
      5.099249390e+02, 
     -7.113819280e-02,
      5.317659880e-06,
     -1.963208212e-10,
      2.805268230e-15,
      1.443308939e+07,
     -4.324044462e+03 
   },
}
-- No CEA transport data for NO+, just use NO
db['NO+'].ceaViscosity = db.NO.ceaViscosity  
db['NO+'].ceaThermCond = db.NO.ceaThermCond  
db['NO+'].Hf = {
   value = 990809.704,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.NO2 = {}
db.NO2.atomicConstituents = {O=2,N=1}
db.NO2.charge = 0
db.NO2.M = {
   value = 0.0460055000,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.NO2.gamma = {
   value = 1.287,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.NO2.sigma = {
   value = 3.5,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NO2.epsilon = {
   value = 200.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.NO2.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -5.642038780e+04,
      9.633085720e+02,
     -2.434510974e+00,
      1.927760886e-02,
     -1.874559328e-05,
      9.145497730e-09,
     -1.777647635e-12,
     -1.547925037e+03,
      4.067851210e+01,
   },
   segment1 = {
      7.213001570e+05,
     -3.832615200e+03,
      1.113963285e+01,
     -2.238062246e-03,
      6.547723430e-07,
     -7.611335900e-11,
      3.328361050e-15,
      2.502497403e+04,
     -4.305130040e+01,
   }
}
db.NO2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.39440312E+01,
     -0.15854290E-02,
      0.16657812E-04,
     -0.20475426E-07,
      0.78350564E-11,
      0.28966179E+04,
      0.63119917E+01,
   },
   segment1 = {
      0,
      0,
      0.48847542E+01,
      0.21723956E-02,
     -0.82806906E-06,
      0.15747510E-09,
     -0.10510895E-13,
      0.23164983E+04,
     -0.11741695E+00,
   }
}
db.NO2.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A =  0.57379100e+00,
      B = -0.12636034e+03,
      C =  0.21566823e+04,
      D =  0.22287492e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.64239645e+00,
      B =  0.60012144e+00,
      C = -0.27020876e+05,
      D = 0.16570566e+01
   }
}
db.NO2.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 300.0,
      T_upper = 1000.0,
      A =  0.48574998e+00,
      B = -0.50702110e+03,
      C =  0.46605820e+05,
      D =  0.36444556e+01
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  0.97660465e+00,
      B =  0.72760751e+03,
      C = -0.32527989e+06,
      D = -0.60899123e+00
   }
}
db.NO2.Hf = {
   value = 34193.019,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.O = {}
db.O.type = "atom"
db.O.atomicConstituents = {O=1,}
db.O.charge = 0
db.O.M = {
   value = 0.01599940,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.O.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.O.sigma = {
   value = 2.750,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O.epsilon = {
   value = 80.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O.Lewis = {
   value = 0.712
}
db.O.electronic_levels = {
  Te = {
    --value = {77.97, 15867.86, 33792.58, 73768.2, 76794.98, 86629.09},
    value = {0.0, 15867.86, 33792.58, 73768.2, 76794.98, 86629.09},
    --value = {0.0, 158.3, 227.0, 15867.9, 33792.6},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {9, 5, 1, 5, 3, 15},
    --value = {5, 3, 1, 5, 1},
    units = 'cm^(-1)',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db.O.ionisation_energy = {
  value = 13.618055,
  units = "eV",
  description = "Ionisation energy",
  reference = "https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html"
}
db.O.ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      0.00000000e+00,
      0.00000000e+00,
      2.50000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      2.92727578e+04,
      5.20440739e+00
   },
   segment1 = {
     -4.59493000e+05,
      1.36004299e+03,
      9.30779889e-01,
      9.08157296e-04,
     -2.81746654e-07,
      4.49451675e-11,
     -2.68564490e-15,
      2.06167786e+04,
      1.63928984e+01
   },
   segment2 = {
     -1.36180488e+09,
      9.17683684e+05,
     -2.46884942e+02,
      3.48520165e-02,
     -2.61713430e-06,
      9.91544159e-11,
     -1.44192000e-15,
     -7.15328977e+06,
      2.14015002e+03
   },
   segment3 = {
      4.38121309e+11,
     -8.98310623e+07,
      7.33626546e+03,
     -3.05369185e-01,
      6.90318054e-06,
     -8.09257822e-11,
      3.86322459e-16,
      8.10576891e+08,
     -7.16401144e+04
   },
   reference = "Johnston et al. AIAA Paper 2015-3110"
}
db.O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.16826710E+00,
     -3.27931884E-03,
      6.64306396E-06,
     -6.12806624E-09,
      2.11265971E-12,
      2.91222592E+04,
      2.05193346E+00,
   },
   segment1 = {
      0,
      0,
      2.56942078E+00,
     -8.59741137E-05,
      4.19484589E-08,
     -1.00177799E-11,
      1.22833691E-15,
      2.92175791E+04,
      4.78433864E+00,
   }
}
db.O.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7269241e-01,
      B =  8.3842977e+01,
      C = -5.8502098e+04,
      D =  8.5100827e-01
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7669586e-01,
      B =  1.0158420e+03,
      C = -1.0884566e+06,
      D = -1.8001077e-01
   },
}
db.O.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7271664e-01,
      B =  8.3989100e+01,
      C = -5.8580966e+04,
      D =  1.5179900e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7676666e-01,
      B =  1.0170744e+03,
      C = -1.0906690e+06,
      D =  4.8644232e-01
   },
}

db.O.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.785314441989e+01,
      B = 2.112206043926e+00,
      C = -1.984976530122e-01,
      D = 8.935786220530e-03,
   }
}
db.O.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.027818773205e+01,
      B = 2.112205758221e+00,
      C = -1.984976113175e-01,
      D = 8.935784208992e-03,
   }
}

db.O.Hf = {
   value = 249175.003,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['O+'] = {}
db['O+'].type = "atom"
db['O+'].atomicConstituents = {O=1,}
db['O+'].charge = 1
db['O+'].M = {
   value = 15.9988514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['O+'].gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['O+'].electronic_levels = {
  Te = {
    value = {0.0, 91.25, 61871.81, 61903.47, 61944.19, 107807.11},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {2, 4, 2, 4, 6, 6,},
    units = 'NA',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db['O+'].ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      0.00000000e+00,
      0.00000000e+00,
      2.50000016e+00,
      0.00000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      0.00000000e+00,
      1.87935284e+05,
      4.39337767e+00
   },
   segment1 = {
     -3.14817544e+05,
      9.52751354e+02,
      1.37960003e+00,
      6.50225275e-04,
     -1.93781966e-07,
      2.73210815e-11,
     -1.29823418e-15,
      1.81897151e+05,
      1.23661521e+01
   },
   segment2 = {
     -1.89172344e+08,
      1.32258091e+05,
     -3.33573784e+01,
      4.61791665e-03,
     -2.81231154e-07,
      8.25139893e-12,
     -9.46296238e-17,
     -8.44265776e+05,
      3.12572886e+02
   },
   segment3 = {
     -1.25053902e+09,
      1.95862039e+05,
     -7.04252592e+00,
      2.94936789e-04,
     -7.50877176e-09,
      1.46522532e-13,
     -1.07948513e-18,
     -1.65453653e+06,
      1.03404000e+02
   },
   reference="Johnston et al. AIAA Paper 2015-3110"
}
-- No CEA transport data, just use O
db['O+'].ceaViscosity = db.O.ceaViscosity
db['O+'].ceaThermCond = db.O.ceaThermCond 

db['O+'].Hf = {
   value = 1568787.228,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.O2 = {}
db.O2.type = "molecule"
db.O2.molecule_type = "linear"
db.O2.theta_v = {
   value = 2273.53,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Capitelli et al (2005), Table 25. omega_e in ground state converted to K'
}
db.O2.theta_D = {
  value = 78740.00,
  units = 'K',
  description = 'characteristic dissociation temperature',
  reference = '???'
}
db.O2.electronic_levels = {
  g = {
    value = {3, 2, 1, 1, 6, 3},
    units = 'NA',
    description = 'degeneracy of electronic levels',
    reference = 'NIST'
  },
  Te = {
    value = {0.0, 7918.04, 13195.10, 33057.3, 34690.3, 35396.6},
    units = 'cm^(-1)',
    description = 'electronic excitation energy',
    reference = 'NIST'
  }
}
db.O2.atomicConstituents = {O=2,}
db.O2.charge = 0
db.O2.M = {
   value = 0.03199880,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.O2.gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.O2.sigma = {
   value = 3.458,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O2.epsilon = {
   value = 107.400,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.O2.Lewis = {
   value = 1.086
}
db.O2.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -3.425563420e+04,
      4.847000970e+02,
      1.119010961e+00,
      4.293889240e-03,
     -6.836300520e-07,
     -2.023372700e-09,
      1.039040018e-12,
     -3.391454870e+03,
      1.849699470e+01,
   },
   segment1 = {
     -1.037939022e+06,
      2.344830282e+03,
      1.819732036e+00,
      1.267847582e-03,
     -2.188067988e-07,
      2.053719572e-11,
     -8.193467050e-16,
     -1.689010929e+04,
      1.738716506e+01,
   },
   segment2 = {
      4.975294300e+08,
     -2.866106874e+05,
      6.690352250e+01,
     -6.169959020e-03,
      3.016396027e-07,
     -7.421416600e-12,
      7.278175770e-17,
      2.293554027e+06,
     -5.530621610e+02,
   },
}
db.O2.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.78245636E+00,
     -2.99673416E-03,
      9.84730201E-06,
     -9.68129509E-09,
      3.24372837E-12,
     -1.06394356E+03,
      3.65767573E+00,
   },
   segment1 = {
      0,
      0,
      3.28253784E+00,
      1.48308754E-03,
     -7.57966669E-07,
      2.09470555E-10,
     -2.16717794E-14,
     -1.08845772E+03,
      5.45323129E+00,
   }
}
db.O2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0916180e-01,
      B = -5.2244847e+01,
      C = -5.9974009e+02,
      D =  2.0410801e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.2216486e-01,
      B =  1.7550839e+02,
      C = -5.7974816e+04,
      D =  1.0901044e+00
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.3981127e-01,
      B =  3.9194906e+02,
      C = -3.7833168e+05,
      D =  9.0931780e-01
   },
}
db.O2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.7229167e-01,
      B =  6.8463210e+00,
      C = -5.8933377e+03,
      D =  1.2210365e+00
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.0917351e-01,
      B =  2.9124182e+02,
      C = -7.9650171e+04,
      D =  6.4851631e-02
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.1218262e+00,
      B = -1.9286378e+04,
      C =  2.3295011e+07,
      D =  2.0342043e+01
   },
}

db.O2.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.000056657134e+01,
      B = 2.898260367046e+00,
      C = -3.018511062902e-01,
      D = 1.346540184145e-02,
   }
}
db.O2.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.172264668036e+01,
      B = 2.184456084794e+00,
      C = -1.756136715791e-01,
      D = 7.271800113398e-03,
   }
}

db.O2.Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db['O2+'] = {}
db['O2+'].type = "molecule"
db['O2+'].molecule_type = "linear"
db['O2+'].theta_v = {
   value = 2273.53,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Copied from O2'
}
db['O2+'].atomicConstituents = {O=2,}
db['O2+'].charge = 1
db['O2+'].M = {
   value = 31.9982514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['O2+'].gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db['O2+'].electronic_levels = {
  Te = {
    value = {0.0, 32964.00, 40669.00, 49552.00},
    units = 'cm^(-1)',
    description = 'Electronic energy levels',
    reference = 'NIST'
  },
  g = {
    value = {4, 8, 4, 4},
    units = 'NA',
    description = 'Degeneracy of electronic energy levels',
    reference = 'NIST'
  }
}
db['O2+'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -8.607205450e+04,
      1.051875934e+03,
     -5.432380470e-01,
      6.571166540e-03,
     -3.274263750e-06,
      5.940645340e-11,
      3.238784790e-13,
      1.345544668e+05,
      2.902709750e+01
   },
   segment1 = {
      7.384654880e+04,
     -8.459559540e+02,
      4.985164160e+00,
     -1.611010890e-04,
      6.427083990e-08,
     -1.504939874e-11,
      1.578465409e-15,
      1.446321044e+05,
     -5.811230650e+00
   },
   segment2 = {
     -1.562125524e+09,
      1.161406778e+06,
     -3.302504720e+02,
      4.710937520e-02,
     -3.354461380e-06,
      1.167968599e-10,
     -1.589754791e-15,
     -8.857866270e+06,
      2.852035602e+03
   },
}
-- No transport properties in CEA, just set to O2
db['O2+'].ceaViscosity = db.O2.ceaViscosity
db['O2+'].ceaThermCond = db.O2.ceaThermCond

db['O2+'].Hf = {
   value = 1171828.436,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.O3 = {}
db.O3.atomicConstituents = {O=3}
db.O3.M = {
   value = 47.9982e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.O3.charge = 0
db.O3.gamma = {
   value = 1.267,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.O3.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.282314507e+04,
      5.898216640e+02,
     -2.547496763e+00,
      2.690121526e-02,
     -3.528258340e-05,
      2.312290922e-08,
     -6.044893270e-12,
      1.348368701e+04,
      3.852218580e+01,
   },
   segment1 = {
     -3.869662480e+07,
      1.023344994e+05,
     -8.961551600e+01,
      3.706144970e-02,
     -4.137638740e-06,
     -2.725018591e-10,
      5.248188110e-14,
     -6.517918180e+05,
      7.029109520e+02,
   }
}
db.O3.Hf = {
   value = 141800.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.OH = {}
db.OH.atomicConstituents = {O=1,H=1,}
db.OH.charge = 0
db.OH.M = {
   value = 0.01700734,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.OH.gamma = {
   value = 1.38600000,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.OH.sigma = {
   value = 2.750,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.OH.epsilon = {
   value = 80.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.OH.Lewis = {
   value = 0.736
}
db.OH.ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -1.998858990e+03,
      9.300136160e+01,
      3.050854229e+00,
      1.529529288e-03,
     -3.157890998e-06,
      3.315446180e-09,
     -1.138762683e-12,
      2.991214235e+03,
      4.674110790e+00,
   },
   segment1 = {
      1.017393379e+06,
     -2.509957276e+03,
      5.116547860e+00,
      1.305299930e-04,
     -8.284322260e-08,
      2.006475941e-11,
     -1.556993656e-15,
      2.019640206e+04,
     -1.101282337e+01,
   },
   segment2 = {
      2.847234193e+08,
     -1.859532612e+05,
      5.008240900e+01,
     -5.142374980e-03,
      2.875536589e-07,
     -8.228817960e-12,
      9.567229020e-17,
      1.468393908e+06,
     -4.023555580e+02,
   },
}
db.OH.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.99201543E+00,
     -2.40131752E-03,
      4.61793841E-06,
     -3.88113333E-09,
      1.36411470E-12,
      3.61508056E+03,
     -1.03925458E-01,
   },
   segment1 = {
      0,
      0,
      3.09288767E+00,
      5.48429716E-04,
      1.26505228E-07,
     -8.79461556E-11,
      1.17412376E-14,
      3.85865700E+03,
      4.47669610E+00,
   }
}
db.OH.ceaViscosity = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  5.9711536e-01,
      B = -4.6100678e+02,
      C =  3.7606286e+04,
      D =  2.4041761e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4287721e-01,
      B = -1.8173747e+02,
      C = -8.8543767e+04,
      D =  1.9636057e+00
   },
}
db.OH.ceaThermCond = {
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.8627561e-01,
      B = -7.4033274e+02,
      C =  2.7559033e+04,
      D =  2.8308741e+00
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.7918112e-01,
      B = -9.3769908e+03,
      C =  7.0509952e+06,
      D =  1.4203688e+01
   },
}

db.OH.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.782259788034e+01,
      B = 2.112205682453e+00,
      C = -1.984976002697e-01,
      D = 8.935783676410e-03,
   }
}
db.OH.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -9.245117757609e-01,
      B = -1.760258809859e+00,
      C = 3.295703318949e-01,
      D = -1.399005438687e-02,
   }
}

db.OH.Hf = {
   value = 37278.206,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.SF6 = {}--NONE OF SF6 DATA HAS BEEN UDPATED THIS IS MERELY A PLACED HOLDER
db.SF6.atomicConstituents = {}
db.SF6.charge = 0
db.SF6.M = {
   value = 0.146055,
   units = 'kg/mol',
}
db.SF6.gamma = {
   value = 1.3,
   note = "valid at low temperatures"
}
db.SF6.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.SF6.sutherlandVisc = {
   mu_ref = 14.8e-6, 
   T_ref = 293.15,
   S = 240.0,
   reference = "Crane Company (1988) - Flow of fluids through valves, fittings and pipes"
}
db.SF6.sutherlandThermCond = {
   T_ref = 273.0, --these have not been updated
   k_ref = 0.0241, --these have not been updated
   S = 194.0,--these have not been updated
   reference = "Table 1-3, White (2006)"
}

db.SF6.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      3.309526740e+05,
     -4.737685050e+03,
      2.247738068e+01,
      1.046954309e-02,
     -2.560641961e-05,
      2.153716967e-08,
     -6.516098960e-12,
     -1.255360583e+05,
     -1.091760145e+02
   },
   segment1 = {
     -7.306726500e+05,
     -6.367056550e+02,
      1.947442853e+01,
     -1.894325671e-04,
      4.178722830e-08,
     -4.783744950e-12,
      2.213516129e-16,
     -1.510609837e+05,
     -8.147574587e+01
   }
}
db.SF6.Hf = {
   value = -1219400.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db.aC3H4 = {}
db.aC3H4.atomicConstituents = {C=3,H=4,O=0,}
db.aC3H4.charge = 0
db.aC3H4.M = {   
	 value = 0.040064,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.aC3H4.gamma = {   
	value = 1.1637,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.aC3H4.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.aC3H4.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.aC3H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.613044500000e+00,          
           1.212257500000e-02,          
           1.853988000000e-05,         
          -3.452514900000e-08,         
           1.533507900000e-11,          
           2.154156700000e+04,          
           1.022613900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.316872200000e+00,          
           1.113372800000e-02,          
          -3.962937800000e-06,         
           6.356423800000e-10,         
          -3.787554000000e-14,          
           2.011749500000e+04,          
          -1.099576600000e+01,      
	}
}

db.aC3H5 = {}
db.aC3H5.atomicConstituents = {C=3,H=5,}
db.aC3H5.charge = 0
db.aC3H5.M = {
   value = 41.071800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.aC3H5.gamma = {
   value = 1.1502e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.aC3H5.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.aC3H5.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.aC3H5.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 3000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.13631835E+01,
      0.19813821E-01,
      0.12497060E-04,
     -0.33355555E-07,
      0.15846571E-10,
      0.19245629E+05,
      0.17173214E+02,
   },
   segment1 = {
      0,
      0,
      0.65007877E+01,
      0.14324731E-01,
     -0.56781632E-05,
      0.11080801E-08,
     -0.90363887E-13,
      0.17482449E+05,
     -0.11243050E+02,
   }
}
db.aC3H5.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -4.315996e+04,
      1.441601e+03,
     -1.197014e+01,
      7.319796e-02,
     -9.066358e-05,
      6.077059e-08,
     -1.658826e-11,
      1.232157e+04,
      8.563173e+01,
   },
   segment1 = {
      4.094571e+06,
     -1.676676e+04,
      3.123006e+01,
     -2.885450e-03,
      5.211344e-07,
     -5.058284e-11,
      2.039933e-15,
      1.185720e+05,
     -1.823070e+02
    }
}
db.aC3H5.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 3000.000,
      A = -2.737260561578e+01,
      B = 5.270213437016e+00,
      C = -5.843750571241e-01,
      D = 2.467473070508e-02,
   }
}
db.aC3H5.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 3000.000,
      A = -3.326474626507e+01,
      B = 9.549029292733e+00,
      C = -9.860839719784e-01,
      D = 3.644127580075e-02,
   }
}

db.air = {}
db.air.atomicConstituents = {}
db.air.charge = 0
db.air.M = {
   value = 28.96e-3,
   units = 'kg/mol',
}
db.air.gamma = {
   value = 1.4,
   note = "valid at low temperatures"
}
db.air.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.air.sutherlandVisc = {
   mu_ref = 1.716e-5, 
   T_ref = 273.0,
   S = 111.0,
   reference = "Table 1-2, White (2006)"
}
db.air.sutherlandThermCond = {
   T_ref = 273.0, 
   k_ref = 0.0241, 
   S = 194.0,
   reference = "Table 1-3, White (2006)"
}
db.air.ceaThermoCoeffs = {
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      1.009950160e+04,
     -1.968275610e+02,
      5.009155110e+00,
     -5.761013730e-03,
      1.066859930e-05,
     -7.940297970e-09,
      2.185231910e-12,
     -1.767967310e+02,
     -3.921504225e+00
   },
   segment1 = { 
      2.415214430e+05,
     -1.257874600e+03,
      5.144558670e+00,
     -2.138541790e-04,
      7.065227840e-08,
     -1.071483490e-11,
      6.577800150e-16,
      6.462263190e+03,
     -8.147411905e+00
   }
}

db["c-C6H4"] = {}
db["c-C6H4"].atomicConstituents = {C=6,H=4,O=0,}
db["c-C6H4"].charge = 0
db["c-C6H4"].M = {   
	value = 0.076096,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["c-C6H4"].gamma = {   
	value = 1.1145,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["c-C6H4"].sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["c-C6H4"].epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["c-C6H4"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -3.099126800000e+00,          
           5.403056400000e-02,          
          -4.083900400000e-05,         
           1.073883700000e-08,         
           9.807849000000e-13,          
           5.220571100000e+04,          
           3.741520700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.384920900000e+01,          
           7.880792000000e-03,          
           1.824383600000e-06,         
          -2.116916600000e-09,         
           3.745997700000e-13,          
           4.744634000000e+04,          
          -5.040495300000e+01,      
	}
}

db.cC3H4 = {}
db.cC3H4.atomicConstituents = {C=3,H=4,O=0,}
db.cC3H4.charge = 0
db.cC3H4.M = {   
	 value = 0.040064,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.cC3H4.gamma = {   
	value = 1.1849,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.cC3H4.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.cC3H4.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.cC3H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -2.462104700000e-02,          
           2.319721500000e-02,          
          -1.847435700000e-06,         
          -1.592759300000e-08,         
           8.684615500000e-12,          
           3.233413700000e+04,          
           2.272976200000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.699993100000e+00,          
           1.035737200000e-02,          
          -3.455116700000e-06,         
           5.065294900000e-10,         
          -2.668227600000e-14,          
           3.019905100000e+04,          
          -1.337877000000e+01,      
	}
}

db['e-'] = {}
db['e-'].type = "electron"
db['e-'].atomicConstituents = {}
db['e-'].charge = -1
db['e-'].M = {
   value = 0.000548579903e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['e-'].gamma = {
   value = 1.667, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}

db['e-'].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0},
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
     -1.172081224e+01
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
     -1.172081224e+01
   },
   segment2 = {
      0.000000000e+00,
      0.000000000e+00,
      2.500000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
      0.000000000e+00,
     -7.453750000e+02,
     -1.172081224e+01
   },
}


db['e-'].Hf = {
   value = 0.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
db["i-C4H3"] = {}
db["i-C4H3"].atomicConstituents = {C=4,H=3,O=0,}
db["i-C4H3"].charge = 0
db["i-C4H3"].M = {   
	value = 0.051066,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["i-C4H3"].gamma = {   
	value = 1.1112,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["i-C4H3"].sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["i-C4H3"].epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["i-C4H3"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,  
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.083041200000e+00,          
           4.083427400000e-02,          
          -6.215968500000e-05,         
           5.167935800000e-08,         
          -1.702918400000e-11,          
           5.800512900000e+04,          
           1.361746200000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           9.097816500000e+00,          
           9.220711900000e-03,          
          -3.387844100000e-06,         
           4.916049800000e-10,         
          -1.452978000000e-14,          
           5.660057400000e+04,          
          -1.980259700000e+01,      
	}
}

db["i-C4H5"] = {}
db["i-C4H5"].atomicConstituents = {C=4,H=5,O=0,}
db["i-C4H5"].charge = 0
db["i-C4H5"].M = {   
	value = 0.053082,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["i-C4H5"].gamma = {   
	value = 1.1232,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["i-C4H5"].sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["i-C4H5"].epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["i-C4H5"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.993290000000e-02,          
           3.800567200000e-02,          
          -2.755945000000e-05,         
           7.783555100000e-09,         
           4.020938300000e-13,          
           3.749622300000e+04,          
           2.439424100000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.022909200000e+01,          
           9.485013800000e-03,          
          -9.040644500000e-08,         
          -1.259610000000e-09,         
           2.478146800000e-13,          
           3.464281200000e+04,          
          -2.856452900000e+01,      
	}
}

db.iC3H7 = {}
db.iC3H7.atomicConstituents = {C=3,H=7,O=0,}
db.iC3H7.charge = 0
db.iC3H7.M = {   
	 value = 0.043087,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.iC3H7.gamma = {   
	value = 1.1429,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.iC3H7.sigma = {   
	value =  4.982 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.iC3H7.epsilon = {   
	value = 266.800 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.iC3H7.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.444919900000e+00,          
           2.099911200000e-02,          
           7.703622200000e-06,         
          -1.847625300000e-08,         
           7.128296200000e-12,          
           9.422372400000e+03,          
           2.011631700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.519274100000e+00,          
           1.722010400000e-02,          
          -5.736421700000e-06,         
           8.413073200000e-10,         
          -4.456591300000e-14,          
           7.322719300000e+03,          
          -9.083021500000e+00,      
	}
}

db["l-C6H4"] = {}
db["l-C6H4"].atomicConstituents = {C=6,H=4,O=0,}
db["l-C6H4"].charge = 0
db["l-C6H4"].M = {   
	value = 0.076096,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["l-C6H4"].gamma = {   
	value = 1.0855,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["l-C6H4"].sigma = {   
	value =  5.349 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["l-C6H4"].epsilon = {   
	value = 412.300 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["l-C6H4"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.959022500000e-01,          
           5.805331800000e-02,          
          -6.776675600000e-05,         
           4.337676200000e-08,         
          -1.141886400000e-11,          
           6.000137100000e+04,          
           2.231897000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.271518200000e+01,          
           1.383966200000e-02,          
          -4.376544000000e-06,         
           3.154163600000e-10,         
           4.661902600000e-14,          
           5.703114800000e+04,          
          -3.946460000000e+01,      
	}
}

db["n-C4H3"] = {}
db["n-C4H3"].atomicConstituents = {C=4,H=3,O=0,}
db["n-C4H3"].charge = 0
db["n-C4H3"].M = {   
	value = 0.051066,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["n-C4H3"].gamma = {   
	value = 1.1261,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["n-C4H3"].sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["n-C4H3"].epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["n-C4H3"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -3.168411300000e-01,          
           4.691210000000e-02,          
          -6.809381000000e-05,         
           5.317992100000e-08,         
          -1.652300500000e-11,          
           6.247619900000e+04,          
           2.462255900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.432827900000e+00,          
           1.686098100000e-02,          
          -9.431310900000e-06,         
           2.570389500000e-09,         
          -2.745630900000e-13,          
           6.160068000000e+04,          
          -1.567398100000e+00,      
	}
}

db["n-C4H5"] = {}
db["n-C4H5"].atomicConstituents = {C=4,H=5,O=0,}
db["n-C4H5"].charge = 0
db["n-C4H5"].M = {   
	value = 0.053082,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["n-C4H5"].gamma = {   
	value = 1.1185,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["n-C4H5"].sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["n-C4H5"].epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["n-C4H5"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.630532100000e-01,          
           3.983013700000e-02,          
          -3.400012800000e-05,         
           1.514723300000e-08,         
          -2.466582500000e-12,          
           4.142976600000e+04,          
           2.353616300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           9.850197800000e+00,          
           1.077900800000e-02,          
          -1.367212500000e-06,         
          -7.720053500000e-10,         
           1.836631400000e-13,          
           3.884030100000e+04,          
          -2.600184600000e+01,      
	}
}

db.nC3H7 = {}
db.nC3H7.atomicConstituents = {C=3,H=7,}
db.nC3H7.charge = 0
db.nC3H7.M = {
   value = 43.087680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.nC3H7.gamma = {
   value = 1.1314e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.nC3H7.sigma = {
   value = 4.982,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'USC-Mech II transport file.'
}
db.nC3H7.epsilon = {
   value = 266.800,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'USC-Mech II transport file.'
}
db.nC3H7.grimechThermoCoeffs = {
   notes = 'data from USC-Mech II',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 3000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.10491173E+01,
      0.26008973E-01,
      0.23542516E-05,
     -0.19595132E-07,
      0.93720207E-11,
      0.10312346E+05,
      0.21136034E+02,
   },
   segment1 = {
      0,
      0,
      0.77097479E+01,
      0.16031485E-01,
     -0.52720238E-05,
      0.75888352E-09,
     -0.38862719E-13,
      0.79762236E+04,
     -0.15515297E+02,
   }
}
db.nC3H7.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -1.895337e+05,
      3.949517e+03,
     -2.606216e+01,
      1.121920e-01,
     -1.365292e-04,
      9.023663e-08,
     -2.441057e-11,
     -7.227877e+03,
      1.673706e+02,
   },
   segment1 = {
      5.646513e+06,
     -2.291087e+04,
      3.987275e+01,
     -4.106233e-03,
      7.562558e-07,
     -7.478263e-11,
      3.068984e-15,
      1.483007e+05,
     -2.403781e+02
   }
}

db.nC3H7.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 3000.000,
      A = -2.734864796957e+01,
      B = 5.270213096852e+00,
      C = -5.843750081534e-01,
      D = 2.467472836750e-02,
   }
}
db.nC3H7.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 3000.000,
      A = -3.138426993014e+01,
      B = 8.645992725716e+00,
      C = -8.407403308989e-01,
      D = 2.905510973744e-02,
   }
}

db.pC3H4 = {}
db.pC3H4.atomicConstituents = {C=3,H=4,O=0,}
db.pC3H4.charge = 0
db.pC3H4.M = {   
	 value = 0.040064,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.pC3H4.gamma = {   
	value = 1.158,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.pC3H4.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.pC3H4.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.pC3H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.680386900000e+00,          
           1.579965100000e-02,          
           2.507059600000e-06,         
          -1.365762300000e-08,         
           6.615428500000e-12,          
           2.080237400000e+04,          
           9.876935100000e+00,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.025240000000e+00,          
           1.133654200000e-02,          
          -4.022339100000e-06,         
           6.437606300000e-10,         
          -3.829963500000e-14,          
           1.962094200000e+04,          
          -8.604378500000e+00,      
	}
}

