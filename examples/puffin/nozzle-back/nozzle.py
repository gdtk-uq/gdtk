# Axisymmetric conical nozzle from Back et al AIAA J 3(9):1606-1614.
# PJ 2022-03-04

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
# Start with slightly supersonic condition at throat to fudge the transonic there.
gas1.p = 500.0e3 * 0.4684
gas1.T = 300.0 * 0.8052
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 1.1
V1 = M1 * gas1.a
print("V1=", V1)

config.axisymmetric = True
config.max_step_relax = 80

# Define geometry of the supersonic part of the nozzle.
# The original paper by Back etal (for a lab experiment on supersonic nozzles)
# specifies sizes in inches, Puffin works in metres.
inch = 0.0254 # metres
L_subsonic = 3.0 * inch
L_nozzle = 6.0 * inch
R_tube = 1.5955 * inch
R_throat = 0.775 * inch
R_curve = 1.55 * inch # radius of curvature of throat profile
height = R_throat + R_curve
import math
theta = math.radians(15.0) # nominal straight-cone angle

# Derived points at start and end of spline.
p1x = R_curve*math.sin(theta)
p1y = height - R_curve*math.cos(theta)
L_cone = L_nozzle - p1x
p2x = L_nozzle
p2y = p1y + L_cone*math.tan(theta)

# Circular arc profile from throat to start of straight conical profile.
def arc(x):
    theta = math.asin(x/R_curve)
    return height - R_curve*math.cos(theta)

def straight_line(x):
    return p1y + (x-p1x)*math.tan(theta)

from gdtk.geom.xpath import XPath
upper_y = XPath([0.0, p1x, p2x], [arc, straight_line])
def lower_y(x): return 0.0
def lower_bc(x): return 0
def upper_bc(x): return 0

config.max_x = p2x
config.dx = config.max_x/1000
config.plot_dx = config.max_x/100

st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=lower_y, y1=upper_y, bc0=lower_bc, bc1=upper_bc,
                 ncells=60)
