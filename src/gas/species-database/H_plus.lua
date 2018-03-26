db['H+'] = {}
db['H+'].atomicConstituents = {H=1}
db['H+'].charge = 1
db['H+'].M = {
   value = 1.0073914e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
--Use Gamma value for H.
db['H+'].gamma = { 
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

db['H+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 298.15,
      T_upper = 1000.0,
      coeffs = { 
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.840214877e+05,
	-1.140646644e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.840214877e+05,
	-1.140646644e+00
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = { 
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.840214877e+05,
	-1.140646644e+00
      }
   },
   ref="from CEA2::thermo.inp"
}

-- No CEA transport data for H_plus, just use H
db['H+'].ceaViscosity = db.H.ceaViscosity
db['H+'].ceaThermCond = db.H.ceaThermCond


