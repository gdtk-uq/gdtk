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

