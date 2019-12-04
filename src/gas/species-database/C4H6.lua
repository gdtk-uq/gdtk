db.C4H6 = {}
db.C4H6.atomicConstituents = {C=4,H=6,O=0,}
db.C4H6.charge = 0
db.C4H6.M = {   
	 value = 0.05409,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H6.gamma = {   
	value = 1.1216,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H6.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H6.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H6.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.128446500000e-01,          
           3.436902200000e-02,          
          -1.110739200000e-05,         
          -9.210666000000e-09,         
           6.206517900000e-12,          
           1.180227000000e+04,          
           2.308999600000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           8.867313400000e+00,          
           1.491867000000e-02,          
          -3.154871600000e-06,         
          -4.184133000000e-10,         
           1.576125800000e-13,          
           9.133851600000e+03,          
          -2.332817100000e+01,      
	}
}

