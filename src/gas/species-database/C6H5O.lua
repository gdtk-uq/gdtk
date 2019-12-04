db.C6H5O = {}
db.C6H5O.atomicConstituents = {C=6,H=5,O=1,}
db.C6H5O.charge = 0
db.C6H5O.M = {   
	 value = 0.093103,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H5O.gamma = {   
	value = 1.0959,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H5O.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5O.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.821943300000e+00,          
           4.812251000000e-02,          
          -4.679230200000e-06,         
          -3.401859400000e-08,         
           1.864963700000e-11,          
           4.242918000000e+03,          
           3.352619900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.383398400000e+01,          
           1.761840300000e-02,          
          -6.069625700000e-06,         
           9.198817300000e-10,         
          -5.044918100000e-14,          
          -6.921254900000e+02,          
          -5.039299000000e+01,      
	}
}

