db.C5H5 = {}
db.C5H5.atomicConstituents = {C=5,H=5,O=0,}
db.C5H5.charge = 0
db.C5H5.M = {   
	 value = 0.065093,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H5.gamma = {   
	value = 1.1209,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H5.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -9.590371800000e-01,          
           3.139685900000e-02,          
           2.672379400000e-05,         
          -6.894187200000e-08,         
           3.330185600000e-11,          
           3.072912000000e+04,          
           2.907281600000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.084406600000e+01,          
           1.539283700000e-02,          
          -5.563042100000e-06,         
           9.018937100000e-10,         
          -5.415653100000e-14,          
           2.690056600000e+04,          
          -3.525494800000e+01,      
	}
}

