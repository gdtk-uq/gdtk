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
exec(open('./constants.txt').read())
# Read case no.
fp = open('case.txt', 'r');
scale_str = fp.readline().strip()
scale = int(scale_str)
fp.close()
# constants
L = 1.0

# Supersonic flow
rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.0; arhox=1.0; arhoy=0.5; arhoxy=0.0;
u0=800.0; ux=50.0; uy=-30.0; uxy=0.0; aux=1.5; auy=0.6; auxy=0.0;
v0=800.0; vx=-75.0; vy=40.0; vxy=0.0; avx=0.5; avy=2.0/3; avxy=0.0;
p0=1.0e5; px=0.2e5; py=0.5e5; pxy=0.0; apx=2.0; apy=1.0; apxy=0.0

x, y, rho, u, v, p, S = symbols('x y rho u v p S')

# Scaling
if scale == 0:
    S = 1.0
else:
    S = exp(-16.0*((x-L/2)*(x-L/2) + (y-L/2)*(y-L/2))/(L*L))

rho = rho0 + S*rhox*sin(arhox*pi*x/L) + S*rhoy*cos(arhoy*pi*y/L) + \
   S*rhoxy*cos(arhoxy*pi*x*y/(L*L));
u =  u0 + S*ux*sin(aux*pi*x/L) + S*uy*cos(auy*pi*y/L) + S*uxy*cos(auxy*pi*x*y/(L*L));
v =  v0 + S*vx*cos(avx*pi*x/L) + S*vy*sin(avy*pi*y/L) + S*vxy*cos(avxy*pi*x*y/(L*L));
p =  p0 + S*px*cos(apx*pi*x/L) + S*py*sin(apy*pi*y/L) + S*pxy*sin(apxy*pi*x*y/(L*L));

def ref_function(x1, y1, z1, t):
    inp = {x:x1, y:y1}
    rho1 = rho.subs(inp).evalf()
    p1 = p.subs(inp).evalf()
    T1 = p1/(rho1*R_air)
    u1 = u.subs(inp).evalf()
    v1 = v.subs(inp).evalf()
    return {"rho":rho1, "p":p1, "T":T1, "vel.x":u1, "vel.y":v1}

if __name__ == "__main__":
    pt = {x:0.5, y:0.5}
    print 'rho=', rho.subs(pt).evalf(), \
        'u=', u.subs(pt).evalf(), \
        'v=', v.subs(pt).evalf(), \
        'p=', p.subs(pt).evalf()
        
