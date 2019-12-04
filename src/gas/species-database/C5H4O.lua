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

