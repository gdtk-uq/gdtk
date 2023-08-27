# analytic_solution.py
# Python version of the analytic solution described in Appendix A of
# C.J. Roy, C.C. Nelson, T.M. Smith and C.C. Ober
# Verification of Euler/Navier-Stokes codes using the method
# of manufactured solutions.
# Int J for Numerical Methods in Fluids 2004; 44:599-620
#
# PJ, 28-May-2011
# It essentially Rowan's code with more and renamed variables
# to bring it closer to the original paper.
# PJ, 30-June-2012
# Scale the disturbance to reduce its magnitude away from the centre.
# RJG, 06-June-2014
# Re-worked completely to use sympy

from __future__ import print_function
from sympy import *
exec(open('constants.txt').read())
L = 1.0
# Supersonic flow
rho0=1.0; rhox=0.15; rhoy=-0.1; rhoz=-0.12; rhoxy=0.0; rhoxz=0.0; rhoyz=0.0; rhoxyz=0.0;
arhox=1.0; arhoy=0.5; arhoz=1.5; arhoxy=0.0; arhoxz=0.0; arhoyz=0.0; arhoxyz=0.0;
u0=800.0; ux=50.0; uy=-30.0; uz=-18.0; uxy=0.0; uxz=0.0; uyz=0.0; uxyz=0.0;
aux=1.5; auy=0.6; auz=0.5; auxy=0.0; auxz=0.0; auyz=0.0; auxyz=0.0;
v0=800.0; vx=-75.0; vy=40.0; vz=-30.0; vxy=0.0; vxz=0.0; vyz=0.0; vxyz=0.0;
avx=0.5; avy=2.0/3.0; avz=1.25; avxy=0.0; avxz=0.0; avyz=0.0; avxyz=0.0;
w0=800.0; wx=15.0; wy=-25.0; wz=35.0; wxy=0.0; wxz=0.0; wyz=0.0; wxyz=0.0;
awx=1.0/3.0; awy=1.5; awz=1.0; awxy=0.0; awxz=0.0; awyz=0.0; awxyz=0.0;
p0=1.0e5; px=0.2e5; py=0.5e5; pz=-0.35e5; pxy=0.0; pxz=0.0; pyz=0.0; pxyz=0.0;
apx=2.0; apy=1.0; apz=1.0/3.0; apxy=0.0; apxz=0.0; apyz=0.0; apxyz=0.0;
    
x, y, z, rho, u, v, w, p, S = symbols('x y z rho u v w p S')

S = 1.0

rho = rho0 + S*rhox*sin(arhox*pi*x/L) + S*rhoy*cos(arhoy*pi*y/L) + \
    S*rhoz*cos(arhoz*pi*z/L) + S*rhoxy*cos(arhoxy*pi*x*y/(L*L)) + \
    S*rhoxz*sin(arhoxz*pi*x*z/(L*L)) + S*rhoyz*cos(arhoyz*pi*y*z/(L*L)) + \
    S*rhoxyz*sin(arhoxyz*pi*x*y*z/(L*L*L));
u =  u0 + S*ux*sin(aux*pi*x/L) + S*uy*cos(auy*pi*y/L) + S*uz*cos(auz*pi*z/L) + \
    S*uxy*cos(auxy*pi*x*y/(L*L)) + S*uxz*sin(auxz*pi*x*z/(L*L)) + \
    S*uyz*cos(auyz*pi*y*z/(L*L)) + S*uxyz*sin(auxyz*pi*x*y*z/(L*L*L));
v =  v0 + S*vx*sin(avx*pi*x/L) + S*vy*cos(avy*pi*y/L) + S*vz*cos(avz*pi*z/L) + \
    S*vxy*cos(avxy*pi*x*y/(L*L)) + S*vxz*sin(avxz*pi*x*z/(L*L)) + \
    S*vyz*sin(avyz*pi*y*z/(L*L)) + S*vxyz*cos(avxyz*pi*x*y*z/(L*L*L));
w =  w0 + S*wx*sin(awx*pi*x/L) + S*wy*cos(awy*pi*y/L) + S*wz*cos(awz*pi*z/L) + \
    S*wxy*cos(awxy*pi*x*y/(L*L)) + S*wxz*cos(awxz*pi*x*z/(L*L)) + \
    S*wyz*sin(awyz*pi*y*z/(L*L)) + S*wxyz*sin(awxyz*pi*x*y*z/(L*L*L));
p =  p0 + S*px*cos(apx*pi*x/L) + S*py*sin(apy*pi*y/L) +  S*pz*sin(apz*pi*z/L) + \
    S*pxy*sin(apxy*pi*x*y/(L*L)) + S*pxz*sin(apxz*pi*x*z/(L*L)) + \
    S*pyz*cos(apyz*pi*y*z/(L*L)) + S*pxyz*cos(apxyz*pi*x*y*z/(L*L*L));

def ref_function(x1, y1, z1, t):
    inp = {x:x1, y:y1, z:z1}
    rho1 = rho.subs(inp).evalf()
    p1 = p.subs(inp).evalf()
    T1 = p1/(rho1*R_air)
    u1 = u.subs(inp).evalf()
    v1 = v.subs(inp).evalf()
    w1 = w.subs(inp).evalf()
    return {"rho":rho1, "p":p1, "T":T1, "vel.x":u1, "vel.y":v1, "vel.z":w1}

if __name__ == "__main__":
    pt = {x:0.5, y:0.5, z:0.5}
    print('rho=', rho.subs(pt).evalf(), \
        'u=', u.subs(pt).evalf(), \
        'v=', v.subs(pt).evalf(), \
        'w=', w.subs(pt).evalf(), \
        'p=', p.subs(pt).evalf())
        
        
