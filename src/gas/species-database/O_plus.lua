db.O_plus = {}
db.O_plus.atomicConstituents = {O=1,}
db.O_plus.charge = 1
db.O_plus.M = {
   value = 15.9988514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db.O_plus.gamma = {
   value = 1.66666667,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db.O_plus.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.15,
      T_upper = 1000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 1.879352842e+05,
	 4.393376760e+00
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	-2.166513208e+05,
	 6.665456150e+02,
	 1.702064364e+00,  
	 4.714992810e-04,
	-1.427131823e-07,
	 2.016595903e-11,
	-9.107157762e-16,
	 1.837191966e+05,
	 1.005690382e+01
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	-2.143835383e+08,
	 1.469518523e+05,
	-3.680864540e+01,
	 5.036164540e-03,
	-3.087873854e-07,
	 9.186834870e-12,
	-1.074163268e-16,
	-9.614208960e+05,
	 3.426193080e+02
      }
   },
}
-- No CEA transport data, just use O
db.O_plus.ceaViscosity = db.O.ceaViscosity
db.O_plus.ceaThermCond = db.O.ceaThermCond 

