db["n-C4H5"] = {}
db["n-C4H5"].atomicConstituents = {C=4,H=5,O=0,}
db["n-C4H5"].charge = 0
db["n-C4H5"].M = {   
	value = 0.053082,   
	units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db["n-C4H5"].gamma = {   
	value = 1.1185,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db["n-C4H5"].sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["n-C4H5"].epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db["n-C4H5"].grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.630532100000e-01,          
           3.983013700000e-02,          
          -3.400012800000e-05,         
           1.514723300000e-08,         
          -2.466582500000e-12,          
           4.142976600000e+04,          
           2.353616300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           9.850197800000e+00,          
           1.077900800000e-02,          
          -1.367212500000e-06,         
          -7.720053500000e-10,         
           1.836631400000e-13,          
           3.884030100000e+04,          
          -2.600184600000e+01,      
	}
}

