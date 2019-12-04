db.C4H612 = {}
db.C4H612.atomicConstituents = {C=4,H=6,O=0,}
db.C4H612.charge = 0
db.C4H612.M = {   
	 value = 0.05409,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H612.gamma = {   
	value = 1.1148,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H612.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H612.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H612.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.023467000000e+00,          
           3.495919000000e-02,          
          -2.200905000000e-05,         
           6.942272000000e-09,         
          -7.879187000000e-13,          
           1.811799000000e+04,          
           1.975066000000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.781557000000e+01,          
          -4.257502000000e-03,          
           1.051185000000e-05,         
          -4.473844000000e-09,         
           5.848138000000e-13,          
           1.267342000000e+04,          
          -6.982662000000e+01,      
	}
}

