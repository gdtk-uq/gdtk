db["i-C4H5"] = {}
db["i-C4H5"].atomicConstituents = {C=4,H=5,O=0,}
db["i-C4H5"].charge = 0
db["i-C4H5"].M = {   
	value = 0.053082,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["i-C4H5"].gamma = {   
	value = 1.1232,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["i-C4H5"].sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["i-C4H5"].epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["i-C4H5"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.993290000000e-02,          
           3.800567200000e-02,          
          -2.755945000000e-05,         
           7.783555100000e-09,         
           4.020938300000e-13,          
           3.749622300000e+04,          
           2.439424100000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.022909200000e+01,          
           9.485013800000e-03,          
          -9.040644500000e-08,         
          -1.259610000000e-09,         
           2.478146800000e-13,          
           3.464281200000e+04,          
          -2.856452900000e+01,      
	}
}

