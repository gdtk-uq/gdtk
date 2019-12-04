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

