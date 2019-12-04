db["c-C6H4"] = {}
db["c-C6H4"].atomicConstituents = {C=6,H=4,O=0,}
db["c-C6H4"].charge = 0
db["c-C6H4"].M = {   
	value = 0.076096,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["c-C6H4"].gamma = {   
	value = 1.1145,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["c-C6H4"].sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["c-C6H4"].epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["c-C6H4"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -3.099126800000e+00,          
           5.403056400000e-02,          
          -4.083900400000e-05,         
           1.073883700000e-08,         
           9.807849000000e-13,          
           5.220571100000e+04,          
           3.741520700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.384920900000e+01,          
           7.880792000000e-03,          
           1.824383600000e-06,         
          -2.116916600000e-09,         
           3.745997700000e-13,          
           4.744634000000e+04,          
          -5.040495300000e+01,      
	}
}

