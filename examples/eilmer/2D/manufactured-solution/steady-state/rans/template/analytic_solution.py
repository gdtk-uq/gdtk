# analytic_solution.py
# C.J. Roy, E. Tendean, S. Veluri, R. Rifki, S. Herbet and E. Luke
# Verification of RANS turbulence models in Loci-CHEM using the method
# of manufactured solutions
# AIAA Computational Fluid Dynamics Conference, 2007, 2007-4203
#
# KD, 2018
# It is essentially Rowan's original code formulated for a RANS solution

from __future__ import print_function
from sympy import *
exec(open('./constants.txt').read())

rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
u0=70.0; ux=7.0; uy=-8.0; uxy=5.5; aux=1.5; auy=1.5; auxy=0.6;
v0=90.0; vx=-5.0; vy=10.0; vxy=-11.0; avx=1.5; avy=1.0; avxy=0.9;
p0=1.0e5; px=0.2e5; py=0.175e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75;
tke0=780.0; tkex=160.0; tkey=-120.0; tkexy=80.0; atkex=0.65; atkey=0.7; atkexy=0.8;
omega0=150.0; omegax=-30.0; omegay=22.5; omegaxy=40.0; a_omegax=0.75; a_omegay=0.875; a_omegaxy=0.6

x, y, rho, u, v, p, tke, omega = symbols('x y rho u v p tke omega')

rho = rho0 + rhox*cos(arhox*pi*x/L) + rhoy*sin(arhoy*pi*y/L) + \
   rhoxy*cos(arhoxy*pi*x*y/(L*L));
u =  u0 + ux*sin(aux*pi*x/L) + uy*cos(auy*pi*y/L) + uxy*cos(auxy*pi*x*y/(L*L));
v =  v0 + vx*sin(avx*pi*x/L) + vy*cos(avy*pi*y/L) + vxy*cos(avxy*pi*x*y/(L*L));
p =  p0 + px*cos(apx*pi*x/L) + py*sin(apy*pi*y/L) + pxy*sin(apxy*pi*x*y/(L*L));
tke =  tke0 + tkex*cos(atkex*pi*x/L) + tkey*sin(atkey*pi*y/L) + tkexy*cos(atkexy*pi*x*y/(L*L));
omega = omega0 + omegax*cos(a_omegax*pi*x/L) + omegay*sin(a_omegay*pi*y/L) + omegaxy*cos(a_omegaxy*pi*x*y/(L*L));
