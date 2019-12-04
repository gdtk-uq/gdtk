db.C2H3CHO = {}
db.C2H3CHO.atomicConstituents = {C=3,H=4,O=1,}
db.C2H3CHO.charge = 0
db.C2H3CHO.M = {   
	 value = 0.056063,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C2H3CHO.gamma = {   
	value = 1.1378,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C2H3CHO.sigma = {   
	value =  5.176 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2H3CHO.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C2H3CHO.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {298.15, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.152958400000e+00,          
           2.804021400000e-02,          
          -1.507215300000e-05,         
           1.590584200000e-09,         
           8.493037100000e-13,          
          -1.041769400000e+04,          
           2.145327900000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           4.835318000000e+00,          
           1.977260100000e-02,          
          -1.042662800000e-05,         
           2.652580300000e-09,         
          -2.627820700000e-13,          
          -1.155783700000e+04,          
           1.885314400000e+00,      
	}
}

