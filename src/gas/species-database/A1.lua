db.A1 = {}
db.A1.atomicConstituents = {C=6,H=6,O=0,}
db.A1.charge = 0
db.A1.M = {   
	 value = 0.078112,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.A1.gamma = {   
	value = 1.1129,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.A1.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.A1.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.A1.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -4.899868000000e+00,          
           5.980693200000e-02,          
          -3.671008700000e-05,         
           3.274039900000e-09,         
           3.760088600000e-12,          
           9.182457000000e+03,          
           4.409564200000e+01,      
        },   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.724699400000e+01,          
           3.842016400000e-03,          
           8.277623200000e-06,         
          -4.896112000000e-09,         
           7.606454500000e-13,          
           2.664605500000e+03,          
          -7.194517500000e+01,      
	}
}

