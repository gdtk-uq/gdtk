db.C6H5OH = {}
db.C6H5OH.atomicConstituents = {C=6,H=6,O=1,}
db.C6H5OH.charge = 0
db.C6H5OH.M = {   
	 value = 0.094111,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H5OH.gamma = {   
	value = 1.0867,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H5OH.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5OH.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H5OH.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.695653900000e+00,          
           5.227129900000e-02,          
          -7.202405000000e-06,         
          -3.585960300000e-08,         
           2.044907300000e-11,          
          -1.328412100000e+04,          
           3.254216000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.491207300000e+01,          
           1.837813500000e-02,          
          -6.198312800000e-06,         
           9.198322100000e-10,         
          -4.920956500000e-14,          
          -1.837519900000e+04,          
          -5.592410300000e+01,      
	}
}

