db.C4H81 = {}
db.C4H81.atomicConstituents = {C=4,H=8,O=0,}
db.C4H81.charge = 0
db.C4H81.M = {   
	 value = 0.056106,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C4H81.gamma = {   
	value = 1.1073,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C4H81.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H81.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C4H81.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 5000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.181138000000e+00,          
           3.085338000000e-02,          
           5.086524700000e-06,         
          -2.465488800000e-08,         
           1.111019300000e-11,          
          -1.790400400000e+03,          
           2.106246900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.053584100000e+00,          
           3.435050700000e-02,          
          -1.588319700000e-05,         
           3.308966200000e-09,         
          -2.536104500000e-13,          
          -2.139723100000e+03,          
           1.554320100000e+01,      
	}
}

