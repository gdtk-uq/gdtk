-- Author: Rowan J. Gollan
-- Date: 14-Jan-2015
--
-- Sample input data for CEA viscosity class.

cea = {
   nsegments = 3,
   segment0 = {T_lower = 200.0, T_upper = 1000.0,
	       A=0.60916180e+00, B=-0.52244847e+02,
	       C=-0.59974009e+03, D=0.20410801e+01 },
   segment1 = {T_lower = 1000.0, T_upper = 5000.0,
	       A=0.72216486e+00, B=0.17550839e+03,
	       C=-0.57974816e+05, D=0.10901044e+01 },
   segment2 = {T_lower = 5000.0, T_upper = 15000.0,
	       A=0.73981127e+00, B=0.39194906e+03,
	       C=-0.37833168e+06, D=0.90931780e+00 }
}

Sutherland = {
   T_ref = 273.0,
   mu_ref = 1.716e-5;
   S = 111.0;
}
