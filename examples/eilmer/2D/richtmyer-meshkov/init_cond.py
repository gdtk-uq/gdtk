# init_cond.py
# Mach 2 shock -- PJ 2015-05-16
from __future__ import print_function
from cfpylib.gasdyn.ideal_gas_flow import m2_shock, p2_p1, u2_u1, T2_T1, r2_r1
from math import sqrt
Rgas = 287.1
g = 1.4
M1 = 2.1
p1 = 100.0e3 # Pa
T1 = 300.0 # degrees K
a1 = sqrt(Rgas*g*T1)
u1 = M1 * a1
rho1 = p1/(Rgas*T1)
print("rho1=", rho1, "kg/m**3")
print("a1=", a1, "m/s")
print("u1=", u1, "m/s")
print("M2=", m2_shock(M1))
print("p2/p1=", p2_p1(M1))
print("u2/u1=", u2_u1(M1))
print("T2/T1=", T2_T1(M1))
p2 = p2_p1(M1)*p1
print("p2=", p2, "Pa")
T2 = T2_T1(M1)*T1
print("T2=", T2, "K")
ug = (1.0 - u2_u1(M1))*u1
a2 = sqrt(T2*Rgas*g)
print("ug=", ug, "m/s")
print("Mg=", ug/a2)
rho2 = r2_r1(M1) * rho1
print("rho2=", rho2, "kg/m**3")
p3 = p1
rho3 = 3.0 * rho1
T3 = p3 / (Rgas * rho3)
print("T3=", T3, "K")


