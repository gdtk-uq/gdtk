-- A config file to store some parameters of the simulation that are used in several places.

-- Size of domain
L_x = 0.001-- m
L_y = 0.001 -- m 

-- Mesh Refinement
N_refine = 5 -- set integer value (1 --> 1 x 5 cell mesh)

-- Initial conditions
startP = 100000 -- Pa
startT = 300.0 -- K

-- Properties of wall movement
u_max = 50.  -- m/s
omega = 200.*2.*math.pi -- rad/s

