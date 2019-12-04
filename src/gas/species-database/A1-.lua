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

