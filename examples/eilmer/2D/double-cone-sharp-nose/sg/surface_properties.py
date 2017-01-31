#! /usr/bin/env python
# surface_properties.py
#
# Pick up the simulation data at the last simulated time
# compute an estimate of the shear-stress coefficient and
# output both shear and pressure along the model surface.
#
# KD, June-2016
from numpy import *
# inflow conditions ----------------------------------
rho_inf = 6.081e-4 # kg/m**3
p_inf = 18.55 # Pa
u_inf =  2576.0 # m/s
T_inf = 102.2 # K
T_wall = 295.8 # K
#mu_inf = sutherland.mu(T_inf, 'N2')

# geometry -------------------------------------------
mm = 0.001  # metres
corner1_x = 92.08*mm
corner1_y = 42.94*mm
corner2_x = 153.69*mm
corner2_y = 130.925*mm
theta1 = 25.0 * (pi/180.0)
theta2 = 55.0 * (pi/180.0)
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
    if x < corner1_x: # first cone
        # Surface to cell-centre distance.
        y_surface = (x - 0.0) * tan(theta1)
        y_vertical = y - y_surface
        dy = y_vertical*cos(theta1)
    elif x > corner1_x and x < corner2_x: # second cone
        # Surface to cell-centre distance.
        y_surface = (x - corner1_x) * tan(theta2) + corner1_y
        y_vertical = y - y_surface
        dy = y_vertical*cos(theta2)
    else: # flat base
        # Surface to cell-centre distance.
        dy = y - corner2_y
    # Cell-centre flow data.
    rho = rhoArray[i] 
    ux = velxArray[i]
    uy = velyArray[i]
    if x < corner1_x: # first cone
        # velocity component tangent to surface.
        vt = -ux * cos(theta1) + uy * sin(theta1)
    if x > corner1_x and x < corner2_x: # second cone
        # velocity component tangent to surface.
        vt = -ux * cos(theta2) + uy * sin(theta2)
    else: # flat base
        # velocity component tangent to surface.
        vt = ux
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

