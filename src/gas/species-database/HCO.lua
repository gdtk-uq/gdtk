db.HCO = {}
db.HCO.atomicConstituents = {C=1,H=1,O=1,}
db.HCO.charge = 0
db.HCO.M = {
   value = 29.018040e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp',
}
db.HCO.gamma = {
   value = 1.316,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = "Gokel (2004), Dean's Handbook of Organic Chemistry",
}
db.HCO.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCO.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.HCO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   segment0 ={
      T_lower = 200.0,
      T_upper = 1000.0,
      coeffs = {
         0,
         0,
          4.22118584E+00,
         -3.24392532E-03,
          1.37799446E-05,
         -1.33144093E-08,
          4.33768865E-12,
          3.83956496E+03,
          3.39437243E+00,
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 3500.0,
      coeffs = {
         0,
         0,
          2.77217438E+00,
          4.95695526E-03,
         -2.48445613E-06,
          5.89161778E-10,
         -5.33508711E-14,
          4.01191815E+03,
          9.79834492E+00,
      }
   }
}
