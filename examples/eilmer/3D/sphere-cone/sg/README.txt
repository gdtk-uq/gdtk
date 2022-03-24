Example for building a structured mesh around a sphere cone.
The example is based on one of the sphere cone geometries simulated in MECH4480.
The test case used the following parameters:
GEOMETRY:
R_nose = 15 mm
Angle = 20 degree (defined as 90-20=70 in simualtion)
Cone height = 25 mm
FLUID:
Ideal gas air
INFLOW:
T = 384.698 K
p = 21434.397 Pa
v = 1509.033 m/s

During experiments a shock stand-off of ~2.2 mm was measured. 
In 2D, simulations shock stand-off of 2.47 mm at the nose (using thermally
perfect air) or 2.69mm (with ideal air), have been obtained, which is close
enough without conducting an in-depth investigations.
In 3D, with structured grid (limited refinement) shock stand off of
approximately 2.48mm (with ideal air) have been predicted. 

The mesh is for the sphere section is constructed by first considering a square
being mapped to a unit circle using the relationship:
x' = x * sqrt(1 - y^2/2)
y' = y * sqrt(1 - x^2/2)
The inner portion of this mesh, defined by variable nose_fraction is used for
the square forming the nose of the sphere cone. 
The outer part of the unit circle is meshed by adding 4 patches to the outer
edges of the resulting square, which interpolate between the edges and the unit
circle. 
To create the sphere mesh, the resulting 2-D mesh is mapped to the surface of a
sphere using polar coordinates. 

Ingo Jahn 2022-03-25