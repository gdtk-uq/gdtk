-- Neon_plus 
-- Yu Liu, 2018-06-16
db['Ne+'] = {}
db['Ne+'].atomicConstituents = {Ne=1,}
db['Ne+'].charge = 1
db['Ne+'].M = {
   value = 20.1797000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2::thermo.inp'
}
db['Ne+'].gamma = {
   value = 5/3,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'monatomic gas'
}
db['Ne+'].ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower = 298.150,
      T_upper = 1000.0007,
      coeffs = {
	 7.281551480E+04,
        -8.695697990E+02, 
         6.108646970E+00,
        -5.841356930E-03, 
         5.041044170E-06,
        -2.293759207E-09, 
         4.339065680E-13,                 
         2.545996890E+05,
        -1.673449355E+01
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = {
        -1.112742658E+05, 
         4.765697970E+02, 
         2.196650531E+00, 
         1.102593151E-04,
        -2.287564425E-08,
         2.510218183E-12,
        -1.126646096E-16,                
         2.472536944E+05, 
         7.466140540E+00
      }
   },
   segment2 = {
      T_lower = 6000.0,
      T_upper = 20000.0,
      coeffs = {
        -5.615474110E+04, 
         1.418980160E+02,
         2.475716842E+00, 
         1.944430992E-06,
        -6.323099200E-11,
        -1.313313446E-16, 
         3.534699010E-20,               
         2.494452217E+05,
         5.366882220E+00
      }
   },
   reference = 'cea2::thermo.inp'
}
-- No CEA transport data for Ne+, just use Ne
db['Ne+'].ceaViscosity = db.Ne.ceaViscosity 
db['Ne+'].ceaThermCond = db.Ne.ceaThermCond 
