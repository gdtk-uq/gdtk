db.C5H4OH = {}
db.C5H4OH.atomicConstituents = {C=5,H=5,O=1,}
db.C5H4OH.charge = 0
db.C5H4OH.M = {   
	 value = 0.081092,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H4OH.gamma = {   
	value = 1.0946,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H4OH.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H4OH.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H4OH.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.282223600000e+00,          
           4.904116000000e-02,          
          -1.368899700000e-05,         
          -2.913385800000e-08,         
           1.900696400000e-11,          
           8.008709800000e+03,          
           3.079835800000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.336791200000e+01,          
           1.520578500000e-02,          
          -5.459225800000e-06,         
           8.813486600000e-10,         
          -5.277445400000e-14,          
           3.841150600000e+03,          
          -4.592083900000e+01,      
	}
}

