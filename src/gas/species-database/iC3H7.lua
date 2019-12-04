db.iC3H7 = {}
db.iC3H7.atomicConstituents = {C=3,H=7,O=0,}
db.iC3H7.charge = 0
db.iC3H7.M = {   
	 value = 0.043087,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.iC3H7.gamma = {   
	value = 1.1429,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.iC3H7.sigma = {   
	value =  4.982 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.iC3H7.epsilon = {   
	value = 266.800 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.iC3H7.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.444919900000e+00,          
           2.099911200000e-02,          
           7.703622200000e-06,         
          -1.847625300000e-08,         
           7.128296200000e-12,          
           9.422372400000e+03,          
           2.011631700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.519274100000e+00,          
           1.722010400000e-02,          
          -5.736421700000e-06,         
           8.413073200000e-10,         
          -4.456591300000e-14,          
           7.322719300000e+03,          
          -9.083021500000e+00,      
	}
}

