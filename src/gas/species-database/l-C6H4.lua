db["l-C6H4"] = {}
db["l-C6H4"].atomicConstituents = {C=6,H=4,O=0,}
db["l-C6H4"].charge = 0
db["l-C6H4"].M = {   
	value = 0.076096,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["l-C6H4"].gamma = {   
	value = 1.0855,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["l-C6H4"].sigma = {   
	value =  5.349 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["l-C6H4"].epsilon = {   
	value = 412.300 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["l-C6H4"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.959022500000e-01,          
           5.805331800000e-02,          
          -6.776675600000e-05,         
           4.337676200000e-08,         
          -1.141886400000e-11,          
           6.000137100000e+04,          
           2.231897000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.271518200000e+01,          
           1.383966200000e-02,          
          -4.376544000000e-06,         
           3.154163600000e-10,         
           4.661902600000e-14,          
           5.703114800000e+04,          
          -3.946460000000e+01,      
	}
}

