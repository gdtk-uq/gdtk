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

