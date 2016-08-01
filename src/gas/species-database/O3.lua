db.O3 = {}
db.O3.atomicConstituents = {O=3}
db.O3.M = {
   value = 47.9982e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.O3.charge = 0
db.O3.gamma = {
   value = 1.267,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.O3.ceaThermoCoeffs = {
   nsegments = 2,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
	 -1.282314507e+04,
	  5.898216640e+02,
	 -2.547496763e+00,
	  2.690121526e-02,
	 -3.528258340e-05,
	  2.312290922e-08,
	 -6.044893270e-12,
	  1.348368701e+04,
	  3.852218580e+01,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
	 -3.869662480e+07,
	  1.023344994e+05,
	 -8.961551600e+01,
	  3.706144970e-02,
	 -4.137638740e-06,
	 -2.725018591e-10,
	  5.248188110e-14,
	 -6.517918180e+05,
	  7.029109520e+02,
      }
   }
}
