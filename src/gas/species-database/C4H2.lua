db.C4H2 = {}
db.C4H2.atomicConstituents = {C=4,H=2,O=0,}
db.C4H2.charge = 0
db.C4H2.M = {   
	 value = 0.050059,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H2.gamma = {   
	value = 1.1268,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H2.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H2.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -3.920103000000e-01,          
           5.193756500000e-02,          
          -9.173734000000e-05,         
           8.047198600000e-08,         
          -2.689821800000e-11,          
           5.484526600000e+04,          
           2.095779400000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           8.663770800000e+00,          
           6.724718900000e-03,          
          -2.359339700000e-06,         
           3.750638000000e-10,         
          -2.223094000000e-14,          
           5.325227500000e+04,          
          -2.109350300000e+01,      
	}
}

