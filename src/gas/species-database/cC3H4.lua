db.cC3H4 = {}
db.cC3H4.atomicConstituents = {C=3,H=4,O=0,}
db.cC3H4.charge = 0
db.cC3H4.M = {   
	 value = 0.040064,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.cC3H4.gamma = {   
	value = 1.1849,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.cC3H4.sigma = {   
	value =  4.760 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.cC3H4.epsilon = {   
	value = 252.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.cC3H4.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -2.462104700000e-02,          
           2.319721500000e-02,          
          -1.847435700000e-06,         
          -1.592759300000e-08,         
           8.684615500000e-12,          
           3.233413700000e+04,          
           2.272976200000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           6.699993100000e+00,          
           1.035737200000e-02,          
          -3.455116700000e-06,         
           5.065294900000e-10,         
          -2.668227600000e-14,          
           3.019905100000e+04,          
          -1.337877000000e+01,      
	}
}

