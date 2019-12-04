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

