db.C5H5O = {}
db.C5H5O.atomicConstituents = {C=5,H=5,O=1,}
db.C5H5O.charge = 0
db.C5H5O.M = {   
	 value = 0.081092,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C5H5O.gamma = {   
	value = 1.1011,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C5H5O.sigma = {   
	value =  5.29  ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5O.epsilon = {   
	value = 464.8   ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C5H5O.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {200.0, 1000.0, 6000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           2.304283500000e-01,          
           3.232269100000e-02,          
           2.890044300000e-05,         
          -7.067997700000e-08,         
           3.340689100000e-11,          
           8.075308200000e+03,          
           2.533097400000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.260642200000e+01,          
           1.674726000000e-02,          
          -6.109857400000e-06,         
           9.967655700000e-10,         
          -6.011320100000e-14,          
           3.931345500000e+03,          
          -4.260427700000e+01,      
	}
}

