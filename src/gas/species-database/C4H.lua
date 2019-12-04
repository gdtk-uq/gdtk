db.C4H = {}
db.C4H.atomicConstituents = {C=4,H=1,O=0,}
db.C4H.charge = 0
db.C4H.M = {   
	 value = 0.049051,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H.gamma = {   
	value = 1.1418,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.318629500000e+00,          
           3.858295600000e-02,          
          -7.138562300000e-05,         
           6.535635900000e-08,         
          -2.261766600000e-11,          
           9.545610600000e+04,          
           1.556758300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           7.769759300000e+00,          
           4.982997600000e-03,          
          -1.762854600000e-06,         
           2.814428400000e-10,         
          -1.668986900000e-14,          
           9.434590000000e+04,          
          -1.416527400000e+01,      
	}
}

