db.CH3CCH2 = {}
db.CH3CCH2.atomicConstituents = {C=3,H=5,O=0,}
db.CH3CCH2.charge = 0
db.CH3CCH2.M = {   
	 value = 0.041072,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.CH3CCH2.gamma = {   
	value = 1.1463,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.CH3CCH2.sigma = {   
	value =  4.982 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CCH2.epsilon = {   
	value = 266.800 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.CH3CCH2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.732920900000e+00,          
           2.239462000000e-02,          
          -5.149061100000e-06,         
          -6.759646600000e-09,         
           3.825321100000e-12,          
           2.904049800000e+04,          
           1.656887800000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           5.425552800000e+00,          
           1.551107200000e-02,          
          -5.667835000000e-06,         
           7.922438800000e-10,         
          -1.687803400000e-14,          
           2.784302700000e+04,          
          -3.352718400000e+00,      
	}
}

