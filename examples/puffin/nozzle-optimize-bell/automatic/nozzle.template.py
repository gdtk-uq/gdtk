# Axisymmetric nozzle for optimization example.
# Corresponds to the rocket-nozzle exercise in the Eilmer examples.
# PJ 2022-02-08

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
gas1.p = 100.0e3
gas1.T = 300.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 1.01
V1 = M1 * gas1.a
print("V1=", V1)

config.axisymmetric = True
config.max_step_relax = 40

# Define geometry of our rocket motor and nozzle.
# The original paper by Back etal (for a lab experiment on supersonic nozzles)
# specifies sizes in inches, Puffin works in metres.
inch = 0.0254 # metres
L_subsonic = 3.0 * inch
L_nozzle = 6.0 * inch
R_tube = 1.5955 * inch
R_throat = 0.775 * inch
R_curve = 1.55 * inch # radius of curvature of throat profile
height = R_throat + R_curve

# Adjustable parameters.
import math
theta_cone = math.radians($theta_cone) # nominal straight-cone angle
theta_init = math.radians($theta_init) # starting angle for thrust nozzle
alpha = math.radians($alpha) # angle for setting b2(p42) in Bezier curve
beta = math.radians($beta) # angle for setting b3(p43) in Bezier curve

# Derived points at start and end of spline.
p4x = R_curve*math.sin(theta_init)
p4y = height - R_curve*math.cos(theta_init)
L_cone = L_nozzle - p4x
p5x = L_nozzle
p5y = p4y + L_cone*math.tan(theta_cone)
# Interior points of the Bezier curve.
p41x = p4x + 0.2*L_cone; p41y = p4y + 0.2*L_cone*math.tan(theta_init)
p42x = p4x + 0.4*L_cone; p42y = p4y + 0.4*L_cone*math.tan(theta_init+alpha)
p43x = p5x - 0.3*L_cone; p43y = p5y - 0.3*L_cone*math.tan(theta_init-beta)
#
from eilmer.geom.xpath import XBezier
bell = XBezier([p4x,p41x,p42x,p43x,p5x], [p4y,p41y,p42y,p43y,p5y])

# Circular arc profile from throat to start of spline.
def arc(x):
    theta = math.asin(x/R_curve)
    return height - R_curve*math.cos(theta)

from eilmer.geom.xpath import XPath
upper_y = XPath([0.0, p4x, p5x], [arc, bell])
def lower_y(x): return 0.0
def lower_bc(x): return 0
def upper_bc(x): return 0

config.max_x = p5x
config.dx = config.max_x/1000
config.plot_dx = config.max_x/100

st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=lower_y, y1=upper_y, bc0=lower_bc, bc1=upper_bc,
                 ncells=40)
