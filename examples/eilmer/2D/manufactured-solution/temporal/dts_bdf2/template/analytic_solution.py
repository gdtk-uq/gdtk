#
# Python version of the analytic solution described in Table 2 from
#     M.R. Nived, Sai Saketha Chandra Athkuri, Vinayak Eswaran
#     On the application of higher-order Backward Diﬀerence (BDF) methods for computing turbulent ﬂows
#     Computers and Mathematics with Applications, vol. 117 (2022) pg. 299-311
#
# author: Kyle Damm
# date:   22-08-2023
#
from sympy import *
exec(open('constants.txt').read())

# Analytic function coefficients and constants
p0 = 100000.0
u0 = 60.0
v0 = 60.0
r0 = 1.0

x, y, t, p, u, v, rho = symbols('x y t p u v rho')

# Analytic function definition
p = p0*(sin(t)+cos(t))
u = u0*cos(t)
v = v0*cos(t)
rho = r0+0.5*sin(t)

