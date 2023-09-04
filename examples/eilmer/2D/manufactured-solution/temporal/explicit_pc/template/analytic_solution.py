from sympy import *
# constants: read from file
exec(open('constants.txt').read())

# Analytic function coefficients and constants
p0   = 1.0e5
u0   = 60.0
v0   = 30.0
rho0 = 1.0
T0   = 1.0e-03

x, y, t, p, u, v, rho = symbols('x y t p u v rho')

p   =   p0*(1.0 + 0.5*sin(2*pi*t/T0))
u   =   u0*(1.0 + 0.5*cos(2*pi*t/T0))**2
v   =   v0*(1.0 + 0.5*sin(2*pi*t/T0))**2
rho = rho0*(1.0 + 0.5*cos(2*pi*t/T0))

