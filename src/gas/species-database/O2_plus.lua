db['O2+'] = {}
db['O2+'].atomicConstituents = {O=2,}
db['O2+'].charge = 1
db['O2+'].M = {
   value = 31.9982514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db['O2+'].gamma = {
   value = 1.40000000,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db['O2+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.15,
      T_upper = 1000.0,
      coeffs = {
	-8.607205450e+04,
	 1.051875934e+03,
	-5.432380470e-01,
	 6.571166540e-03,
	-3.274263750e-06,
	 5.940645340e-11,
	 3.238784790e-13,
	 1.345544668e+05,
	 2.902709750e+01
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 7.384654880e+04,
	-8.459559540e+02,
	 4.985164160e+00,
	-1.611010890e-04,
	 6.427083990e-08,
	-1.504939874e-11,
	 1.578465409e-15,
	 1.446321044e+05,
	-5.811230650e+00
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
	-1.562125524e+09,
	 1.161406778e+06,
	-3.302504720e+02,
	 4.710937520e-02,
	-3.354461380e-06,
	 1.167968599e-10,
	-1.589754791e-15,
	-8.857866270e+06,
	 2.852035602e+03
      }
   },
}
-- No transport properties in CEA, just set to O2
db['O2+'].ceaViscosity = db.O2.ceaViscosity
db['O2+'].ceaThermCond = db.O2.ceaThermCond

