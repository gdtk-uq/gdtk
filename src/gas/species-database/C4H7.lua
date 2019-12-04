db.C4H7 = {}
db.C4H7.atomicConstituents = {C=4,H=7,O=0,}
db.C4H7.charge = 0
db.C4H7.M = {   
	 value = 0.055098,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H7.gamma = {   
	value = 1.1079,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H7.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H7.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H7.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.869825400000e-01,          
           3.696449500000e-02,          
          -8.627744100000e-06,         
          -1.505182100000e-08,         
           8.989126300000e-12,          
           2.055130100000e+04,          
           2.448446700000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.196339200000e+01,          
           1.142530500000e-02,          
           7.894890900000e-07,         
          -1.985887200000e-09,         
           3.687364500000e-13,          
           1.696297700000e+04,          
          -3.754290800000e+01,      
	}
}

