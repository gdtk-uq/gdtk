db.C3H2 = {}
db.C3H2.atomicConstituents = {C=3,H=2,O=0,}
db.C3H2.charge = 0
db.C3H2.M = {   
	 value = 0.038048,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C3H2.gamma = {   
	value = 1.177,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C3H2.sigma = {   
	value =  4.100 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H2.epsilon = {   
	value = 209.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C3H2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.691077000000e+00,          
           1.480366000000e-02,          
          -3.250551000000e-06,         
          -8.644363000000e-09,         
           5.284878000000e-12,          
           5.219072000000e+04,          
           8.757391000000e+00,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.530853000000e+00,          
           5.870316000000e-03,          
          -1.720777000000e-06,         
           2.127498000000e-10,         
          -8.291910000000e-15,          
           5.115214000000e+04,          
          -1.122728000000e+01,      
	}
}

