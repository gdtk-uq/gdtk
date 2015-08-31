# analytic_solution.py
# Conjugate heat transfer analytic solution first formulated by Rowan 
# and subsequently modified by Anand based on 
# A. Veeraragavan and C.P. Cadou
# Theoretical study of conjugate heat transfer effects on temperature profiles
# in parallel flow with embedded heat sources.
# Int J Heat and Mass Transfer 2010;  53 (9), 1699-1711.
#

from sympy import *
R_air = 287.0


# constants
L = 1.0
k_g = 10000.0
k_s = 10*k_g

# Analytic function coefficients and constants
rho0=1.0; rhox=0.1; rhoy=0.15; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
u0 = 1.0; ux = 0.1; uy = u0; uxy = -ux; aux = 5.0/3; auy = -1.0; auxy = aux;
v0 = 0.9; vx = -0.02; vy = -v0; vxy = -vx; avx = 1.5; avy = 0.5; avxy = avx;	
p0=1.0e5; px=-0.3e5; py=0.2e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75
T0 = 350; Tx = -10.0; Ty = 25.0; aTx = 1.5; aTy = 1.0; Ti = 350.0; aTx2 = 0.75;
x, y, rho, u, v, T, T_s = symbols('x y rho u v T T_s')


rho = rho0 + rhox*sin(arhox*pi*x/L) + rhoy*cos(arhoy*pi*y/L) + \
      rhoxy*cos(arhoxy*pi*x*y/(L*L));
u =  u0 + ux*cos(aux*pi*x/L) + uy*cos(auy*pi*y/L) + uxy*cos(auxy*pi*x*y/(L*L));
v =  v0 + vx*cos(avx*pi*x/L) + vy*sin(avy*pi*y/L) + vxy*cos(avxy*pi*x*y/(L*L));
T = T0 + Tx*cos(aTx*pi*x/L) + Ty*cos(aTx2*pi*x/L)*sin(aTy*pi*y/L);
T_s = T0 + Tx*cos(aTx*pi*x/L) +  Ty*(k_g/k_s)*cos(aTx2*pi*x/L)*sin(aTy*pi*y/L);

def ref_function(x1, y1, z1, t):
    inp = {x:x1, y:y1}
    rho1 = rho.subs(inp).evalf()
    T1 = T.subs(inp).evalf()
    p1 = rho1*R_air*T1	
    u1 = u.subs(inp).evalf()
    v1 = v.subs(inp).evalf()
    return {"rho":rho1, "p":p1, "T[0]":T1, "vel.x":u1, "vel.y":v1}        
