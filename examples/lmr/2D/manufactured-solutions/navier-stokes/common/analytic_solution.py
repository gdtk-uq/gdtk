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
# RJG, 2023-09-21
# Removed scaling so as to simplify this file.
# We do not seem to need the scaling function now
# that we have well-behaved boundary implementations.

from sympy import *
exec(open('./constants.txt').read())

# Subsonic flow
rho0=1.0; rhox=0.1; rhoy=0.15; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
u0=70.0; ux=4.0; uy=-12.0; uxy=7.0; aux=5.0/3; auy=1.5; auxy=0.6;
v0=90.0; vx=-20.0; vy=4.0; vxy=-11.0; avx=1.5; avy=1.0; avxy=0.9;
p0=1.0e5; px=-0.3e5; py=0.2e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75

x, y, rho, u, v, p = symbols('x y rho u v p')

rho = rho0 + rhox*sin(arhox*pi*x/L) + rhoy*cos(arhoy*pi*y/L) + \
   rhoxy*cos(arhoxy*pi*x*y/(L*L));
u =  u0 + ux*sin(aux*pi*x/L) + uy*cos(auy*pi*y/L) + uxy*cos(auxy*pi*x*y/(L*L));
v =  v0 + vx*cos(avx*pi*x/L) + vy*sin(avy*pi*y/L) + vxy*cos(avxy*pi*x*y/(L*L));
p =  p0 + px*cos(apx*pi*x/L) + py*sin(apy*pi*y/L) + pxy*sin(apxy*pi*x*y/(L*L));

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
    print('rho=', rho.subs(pt).evalf(), \
        'u=', u.subs(pt).evalf(), \
        'v=', v.subs(pt).evalf(), \
        'p=', p.subs(pt).evalf())
        
