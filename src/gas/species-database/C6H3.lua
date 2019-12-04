db.C6H3 = {}
db.C6H3.atomicConstituents = {C=6,H=3,O=0,}
db.C6H3.charge = 0
db.C6H3.M = {   
	 value = 0.075088,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H3.gamma = {   
	value = 1.0866,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H3.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H3.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H3.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.179061900000e+00,          
           5.554736000000e-02,          
          -7.307616800000e-05,         
           5.207673600000e-08,         
          -1.504696400000e-11,          
           8.564731200000e+04,          
           1.917919900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.818834300000e+00,          
           2.793340800000e-02,          
          -1.782542700000e-05,         
           5.370253600000e-09,         
          -6.170762700000e-13,          
           8.518825000000e+04,          
          -9.214782700000e-01,      
	}
}

