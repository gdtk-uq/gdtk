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
