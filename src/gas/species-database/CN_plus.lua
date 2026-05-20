-- Ionized CN data compiled by NNG and RW
-- In addition to the usual sources, data has been taken from: 
-- "Tables of Internal Partition Functions and Thermodynamic  Properties of
--  High-Temperature Mars-Atmosphere Species from 50K to 50000K"
-- ESA, STR-256, October 2005
-- M. Capitelli et al.

db["CN+"] = {}
db["CN+"].type = "molecule"
db["CN+"].molecule_type = "linear"
db["CN+"].theta_v = {
   value = 2273.267,
   units = 'K',
   description = 'characteristic vibrational temperature',
   reference = 'Capitelli et al (2005), Table 9. omega_e in ground state converted to K'
}
db["CN+"].theta_D = {
   value = 90051.66,
   units = 'K',
   description = 'characteristic dissociation temperature',
   reference = 'No Data Available (Borrowed from CN)'
}
db["CN+"].atomicConstituents = {C=1,N=1,}
db["CN+"].charge = 1
db["CN+"].M = {
   value = 26.0168514e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db["CN+"].electronic_levels = {
   g = {
      value = {6,1,1,1,1},
      units = "N/A",
      description = "Degeneracy of electronic energy levels",
      reference = "Table 9 Capitelli (2005)",
   },
   Te = {
      value = {0.0, 1500.0, 6000.0, 8000.0, 15000.0},
      units = "cm^(-1)",
      description = "Electronic excitation energy",
      reference = "Table 9 Capitelli (2005)",
   }
}
db["CN+"].gamma = {
   value = 1.399,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db["CN+"].sigma = {
   value = 3.856,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance.',
   reference = 'Borrowed from CN'
}
db["CN+"].epsilon = {
   value = 75.0,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'Borrowed from CN'
}
db["CN+"].SSH_mass_factor = {
   value = 0.9883,
   units = 'unitless',
   description = 'Mass factor = ( M ( Ma^2 + Mb^2 ) / ( 2 Ma Mb ( Ma + Mb ) )',
   reference = 'Borrowed from CN'
}
db["CN+"].r_eq = {
   value = 1.29e-10,
   units = 'm',
   description = 'Equilibrium intermolecular distance',
   reference = 'Capitelli 2005, Table 9. r_e in ground state'
}
db["CN+"].ceaThermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
     -8.302909570e+05,
      8.775687500e+03,
     -2.977443560e+01,
      4.976897060e-02,
     -1.302225951e-05,
     -2.058325353e-08,
      1.126843895e-11,
      1.703860539e+05,
      2.039918818e+02
   },
   segment1 = {
     -7.153463080e+06,
      1.857250421e+04,
     -1.084534159e+01,
      6.106681430e-03,
     -1.191208566e-06,
      1.184848778e-10,
     -4.799838730e-15,
      9.242644960e+04,
      1.135340573e+02
   },
   segment2 = {
     -2.354919695e+08,
      1.433776703e+05,
     -2.975360271e+01,
      4.280545600e-03,
     -2.707260413e-07,
      8.178340660e-12,
     -9.629506200e-17,
     -9.229047140e+05,
      2.964624987e+02
   }
}
db["CN+"].Hf = {
   value = 1798890.904,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
