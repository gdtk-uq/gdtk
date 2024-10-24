-- parameters.lua
-- Extracted from Riggin and Nelson's paper

h = 0.015 -- body thickness, metres
R = h/2 -- easier to define geometry on nose radius
xmax = 0.0165 -- downstream extent of simulation domain

M_inf = 6.5
p_inf = 1185.5 -- Pa
V_inf = 1981.0 -- m/s
rho_inf = 0.01786 -- kg/m^3

Rgas = 287.1 -- J/kg/K for ideal air
T_inf = p_inf/Rgas/rho_inf
