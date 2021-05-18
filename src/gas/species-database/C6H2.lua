db.C6H2 = {}
db.C6H2.atomicConstituents = {C=6,H=2,O=0,}
db.C6H2.charge = 0
db.C6H2.M = {   
	 value = 0.07408,   
	 units = 'kg/mol',   
	description = 'molecular mass',   
	reference = 'Periodic table'
	}
db.C6H2.gamma = {   
	value = 1.0872,   
	units = 'non-dimensional',   
	description = 'ratio of specific heats at 300.0K',   
	reference = 'evaluated using Cp/R from Chemkin-II coefficients'
	}
db.C6H2.sigma = {   
	value =  5.180 ,   
	units = 'Angstrom',   
	description = 'Lennard-Jones potential distance',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H2.epsilon = {   
	value = 357.000 ,   
	units = 'K',   
	description = 'Lennard-Jones potential well depth.',   
	reference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'
	}
db.C6H2.grimechThermoCoeffs = {   
	notes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',   
	nsegments = 2,    
        T_break_points = {300.0, 1000.0, 3000.0},
        T_blend_ranges = {400.0},
	segment0 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
          -1.593262400000e+00,          
           8.053014500000e-02,          
          -1.480064900000e-04,         
           1.330003100000e-07,         
          -4.533231300000e-11,          
           8.327322700000e+04,          
           2.798087300000e+01,      
	},   
	segment1 = {      
           0.000000000000e+00,         
           0.000000000000e+00,          
           1.322628100000e+01,          
           7.390430200000e-03,          
          -2.271538100000e-06,         
           2.587521700000e-10,         
          -5.535674100000e-15,          
           8.056525800000e+04,          
          -4.120117600000e+01,      
	}
}

db.C6H2.Hf = {
   value = 670000.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
