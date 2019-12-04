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

