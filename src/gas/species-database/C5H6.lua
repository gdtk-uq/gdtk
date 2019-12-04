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

