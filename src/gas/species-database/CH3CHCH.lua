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

