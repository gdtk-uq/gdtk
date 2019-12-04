db.aC3H4 = {}
db.aC3H4.atomicConstituents = {C=3,H=4,O=0,}
db.aC3H4.charge = 0
db.aC3H4.M = {   
	 value = 0.040064,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.aC3H4.gamma = {   
	value = 1.1637,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.aC3H4.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.aC3H4.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.aC3H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.613044500000e+00,          
           1.212257500000e-02,          
           1.853988000000e-05,         
          -3.452514900000e-08,         
           1.533507900000e-11,          
           2.154156700000e+04,          
           1.022613900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.316872200000e+00,          
           1.113372800000e-02,          
          -3.962937800000e-06,         
           6.356423800000e-10,         
          -3.787554000000e-14,          
           2.011749500000e+04,          
          -1.099576600000e+01,      
	}
}

