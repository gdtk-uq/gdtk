#! /usr/bin/env python
# surface_properties.py
#
# Pick up the simulation data at the last simulated time
# compute an estimate of the shear-stress coefficient and
# output both shear and pressure along the cylinder and flare.
#
# KD, June-2016
from numpy import *
# inflow conditions ----------------------------------
rho_inf = 8.81e-4 # kg/m**3
p_inf = 31.88 # Pa
u_inf =  2304.0 # m/s
T_inf = 120.4 # K
T_wall = 295.2 # K
#mu_inf = sutherland.mu(T_inf, 'N2')

# geometry -------------------------------------------
mm = 0.001  # metres
R = 32.5*mm
xcorner = 101.7*mm
theta = 30.0 * (pi/180.0)
# corner = is at point (xcorner,R)

# wall conditions ------------------------------------
posX = []; posY = []; pArray = []; tempArray = []
rhoArray = []; velxArray = []; velyArray = []
eArray = []; kArray = []; muArray = []
with open('wallData', 'r') as f:
    data = f.readlines()
    count = 0
    for line in data:
        dat = line.split()
        if (count>0):
            posX.append(float(dat[0]))
            posY.append(float(dat[1]))
            pArray.append(float(dat[8]))
            tempArray.append(float(dat[19]))
            rhoArray.append(float(dat[4]))
            velxArray.append(float(dat[5]))
            velyArray.append(float(dat[6]))
            eArray.append(float(dat[18]))
            kArray.append(float(dat[11]))
            muArray.append(float(dat[10]))
        count += 1

outfile = open("eilmer4-surface.data", "w")
outfile.write("# x(m) Cf p(Pa) Cp q(W/m**2) Ch\n")

for i in range(len(posX)):    
    # Get vertices on surface, for this cell.
    x = posX[i] 
    y = posY[i]
    if x < xcorner:
        # Surface to cell-centre distance.
        dy = y - R
    else:
        # Surface to cell-centre distance.
        y_surface = (x - xcorner) * tan(theta) + R
        y_vertical = y - y_surface
        dy = y_vertical*cos(theta)
    # Cell-centre flow data.
    rho = rhoArray[i] 
    ux = velxArray[i]
    uy = velyArray[i]
    if x < xcorner:
        vt = ux
    else:
        vt = -ux * cos(theta) + uy * sin(theta) # velocity component tangent to surface
    mu = muArray[i]
    kgas = kArray[i]
    p = pArray[i]
    Cp = (p-p_inf)/(0.5*rho_inf*u_inf*u_inf)
    T = tempArray[i]
    # Shear stress
    dudy = (vt - 0.0) / dy # no-slip wall
    tau_w = mu * dudy    # wall shear stress
    Cf = tau_w / (0.5*rho_inf*u_inf*u_inf)
    u_tau = sqrt(abs(tau_w) / rho) # friction velocity
    y_plus = u_tau * dy * rho / mu
    # Heat flux
    dTdy = (T - T_wall) / dy # conductive heat flux at the wall
    q = kgas * dTdy
    Ch = q / (0.5*rho_inf*u_inf*u_inf*u_inf)
    #
    outfile.write("%f %f %f %f %f %f \n" % 
                  (x, Cf, p, Cp, q, Ch))
outfile.close()

