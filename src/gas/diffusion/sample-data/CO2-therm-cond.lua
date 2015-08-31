-- Author: Rowan J. Gollan
-- Date: 19-Jan-2015
--
-- Sample input thermal conductivity classes.

cea = {
   nsegments = 3,
   segment0 = {T_lower = 200.0, T_upper = 1000.0,
	       A=0.48056568e+00, B=-0.50786720e+03,
	       C=0.35088811e+05, D=0.36747794e+01},
   segment1 = {T_lower = 1000.0, T_upper = 5000.0,
	       A=0.69857277E+00, B=-0.11830477E+03,
	       C=-0.50688859E+05, D=0.18650551E+01 },
   segment2 = {T_lower = 5000.0, T_upper = 10000.0,
	       A=0.10518358E+01, B=-0.42555944E+04,
	       C=0.14288688E+08, D=-0.88950473E+00 }
}


