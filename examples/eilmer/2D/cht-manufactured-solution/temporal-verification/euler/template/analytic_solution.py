# analytic_solution.py
#

from sympy import *
# constants: read from file
execfile('constants.txt')

# Analytic function coefficients and constants
TA = 300.0; TB = 0.2;
pConst = 1.0e5
rhoA = 1.0; rhoB = 1.2; rhoD = 1.5; rhoE = 0.9
uA = 0.1; uB = 0.0; uD = 2.8; uE = 0.0
vA = 0.1; vB = 0.0; vD = 0.5; vE = 0.0
Tperiod = 1.0e-3; tA = 1.0; tB = 0.5

x, y, t, p, u, v, T, T_s = symbols('x y t p u v T T_s')

p = (tA + tB * sin(2.0 * pi * t/Tperiod)) * pConst # (tA + tB * sin(2.0 * pi * t/Tperiod)) *
u = (tA + tB * sin(2.0 * pi * t/Tperiod)) *uA*(H-y)
v = (tA + tB * sin(2.0 * pi * t/Tperiod)) *vA*(H-y) 
T = (tA + tB * sin(2.0 * pi * t/Tperiod)) *(TA + TB*(H-y))
T_s = (tA + tB * sin(2.0 * pi * t/Tperiod)) *(TA + (k_g/k_s)*TB*(H-y))


