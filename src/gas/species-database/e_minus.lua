db['e-'] = {}
db['e-'].atomicConstituents = {}
db['e-'].charge = -1
db['e-'].M = {
   value = 0.000548579903e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::therm.inp'
}
db['e-'].gamma = {
   value = 1.667, 
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'Cp/Cv from CEA2 at room temperature'
}
-- The coefficients are the same over all temperature ranges
-- as given by CEA (which makes sense because there is no
-- internal energy storage in an electron), so we can just
-- use one range.
db['e-'].ceaThermoCoeffs = {
   nsegments = 1,
   segment0 = {
      T_lower = 298.15,
      T_upper = 20000.0,
      coeffs = {
	 0.000000000e+00,
	 0.000000000e+00,
	 2.500000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	 0.000000000e+00,
	-7.453750000e+02,
	-1.172081224e+01
      }
   }
}


