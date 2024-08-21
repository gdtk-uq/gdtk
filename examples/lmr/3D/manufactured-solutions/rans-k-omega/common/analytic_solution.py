# analytic_solution.py
# Python version of the analytic solution described in Appendix A of
# S.P. Veluri, C.J. Roy and E.A. Luke
# Comprehensive code verification techniques for finite volume CFD codes
# Computers & Fluids, 2012;
#
# KD, 2018
# It is essentially Rowan's original code formulated for three-dimensional solution

from __future__ import print_function
from sympy import *
exec(open('./constants.txt').read())

rho0=1.0; rhox=0.15; rhoy=-0.1; rhoz=0.1; rhoxy=0.08; rhoyz=0.05; rhozx=0.12; arhox=0.75; arhoy=0.45; arhoz=0.8; arhoxy=0.65; arhoyz=0.75; arhozx=0.5;
u0=70.0; ux=7.0; uy=-15.0; uz=-10.0; uxy=7.0; uyz=4.0; uzx=-4.0; aux=0.5; auy=0.85; auz=0.4; auxy=0.6; auyz=0.8; auzx=0.9;
v0=90.0; vx=-5.0; vy=10.0; vz=5.0; vxy=-11.0; vyz=-5.0; vzx=5.0; avx=0.8; avy=0.8; avz=0.5; avxy=0.9; avyz=0.4; avzx=0.6;
w0=80.0; wx=-10.0; wy=10.0; wz=12.0; wxy=-12.0; wyz=11.0; wzx=5.0; awx=0.85; awy=0.9; awz=0.5; awxy=0.4; awyz=0.8; awzx=0.75;
p0=1.0e5; px=0.2e5; py=0.5e5; pz=0.2e5; pxy=-0.25e5; pyz=-0.1e5; pzx=0.1e5; apx=0.4; apy=0.45; apz=0.85; apxy=0.75; apyz=0.7; apzx=0.8;
tke0=780.0; tkex=160.0; tkey=-120.0; tkez=80.0; tkexy=80.0; tkeyz=60.0; tkezx=-70.0; atkex=0.65; atkey=0.7; atkez=0.8; atkexy=0.8; atkeyz=0.85; atkezx=0.6;
omega0=150.0; omegax=-30.0; omegay=22.5; omegaz=20.0; omegaxy=40.0; omegayz=-15.0; omegazx=25.0; a_omegax=0.75; a_omegay=0.875; a_omegaz=0.65; a_omegaxy=0.6; a_omegayz=0.75; a_omegazx=0.8;

x, y, z, rho, u, v, w, p, tke, omega = symbols('x y z rho u v w p tke omega')

S = 1.0

rho = rho0 + S*rhox*cos(arhox*pi*x/L) + S*rhoy*sin(arhoy*pi*y/L) + \
   S*rhoz*sin(arhoz*pi*z/L) + S*rhoxy*cos(arhoxy*pi*x*y/(L*L)) + \
    S*rhozx*sin(arhozx*pi*x*z/(L*L)) + S*rhoyz*cos(arhoyz*pi*y*z/(L*L));
u =  u0 + S*ux*sin(aux*pi*x/L) + S*uy*cos(auy*pi*y/L) + S*uz*cos(auz*pi*z/L) + \
     S*uxy*cos(auxy*pi*x*y/(L*L)) + S*uzx*sin(auzx*pi*x*z/(L*L)) + \
     S*uyz*cos(auyz*pi*y*z/(L*L));
v =  v0 + S*vx*sin(avx*pi*x/L) + S*vy*cos(avy*pi*y/L) + S*vz*cos(avz*pi*z/L) + \
      S*vxy*cos(avxy*pi*x*y/(L*L)) + S*vzx*sin(avzx*pi*x*z/(L*L)) + \
     S*vyz*cos(avyz*pi*y*z/(L*L));
w =  w0 + S*wx*cos(awx*pi*x/L) + S*wy*sin(awy*pi*y/L) + S*wz*cos(awz*pi*z/L) + \
      S*wxy*sin(awxy*pi*x*y/(L*L)) + S*wzx*sin(awzx*pi*x*z/(L*L)) + \
     S*wyz*cos(awyz*pi*y*z/(L*L));
p =  p0 + S*px*cos(apx*pi*x/L) + S*py*cos(apy*pi*y/L) +  S*pz*sin(apz*pi*z/L) + \
      S*pxy*cos(apxy*pi*x*y/(L*L)) + S*pzx*sin(apzx*pi*x*z/(L*L)) + \
     S*pyz*cos(apyz*pi*y*z/(L*L));
tke =  tke0 + tkex*cos(atkex*pi*x/L) + tkey*cos(atkey*pi*y/L) + tkez*sin(atkez*pi*z/L) + tkexy*cos(atkexy*pi*x*y/(L*L)) + \
       + tkeyz*cos(atkeyz*pi*y*z/(L*L)) + tkezx*sin(atkezx*pi*x*z/(L*L));
omega = omega0 + omegax*cos(a_omegax*pi*x/L) + omegay*cos(a_omegay*pi*y/L) + omegaz*sin(a_omegay*pi*z/L) + omegaxy*cos(a_omegaxy*pi*x*y/(L*L)) + \
        + omegayz*cos(a_omegayz*pi*y*z/(L*L)) + omegazx*sin(a_omegazx*pi*x*z/(L*L));
