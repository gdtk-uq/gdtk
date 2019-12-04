db.pC3H4 = {}
db.pC3H4.atomicConstituents = {C=3,H=4,O=0,}
db.pC3H4.charge = 0
db.pC3H4.M = {   
	 value = 0.040064,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.pC3H4.gamma = {   
	value = 1.158,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.pC3H4.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.pC3H4.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.pC3H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.680386900000e+00,          
           1.579965100000e-02,          
           2.507059600000e-06,         
          -1.365762300000e-08,         
           6.615428500000e-12,          
           2.080237400000e+04,          
           9.876935100000e+00,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.025240000000e+00,          
           1.133654200000e-02,          
          -4.022339100000e-06,         
           6.437606300000e-10,         
          -3.829963500000e-14,          
           1.962094200000e+04,          
          -8.604378500000e+00,      
	}
}

