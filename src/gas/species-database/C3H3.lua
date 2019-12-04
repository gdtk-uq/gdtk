db.C3H3 = {}
db.C3H3.atomicConstituents = {C=3,H=3,O=0,}
db.C3H3.charge = 0
db.C3H3.M = {   
	 value = 0.039056,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C3H3.gamma = {   
	value = 1.1464,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C3H3.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H3.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H3.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.351109270000e+00,          
           3.274112230000e-02,          
          -4.738271350000e-05,         
           3.763098080000e-08,         
          -1.185409230000e-11,          
           4.010577830000e+04,          
           1.520589240000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           7.142218800000e+00,          
           7.619020050000e-03,          
          -2.674599500000e-06,         
           4.249148010000e-10,         
          -2.514754150000e-14,          
           3.890874270000e+04,          
          -1.258484360000e+01,      
	}
}

