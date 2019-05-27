-- A config file to store some parameters of the simulation that are used in several places.

-- set piston mass and height (adjusted from diameter D to height)
D = 0.01 -- m Diameter
H = math.sqrt( D*D /4 * math.pi)-- m  Height to get same area per mass
massPiston = 0.001/H -- kg/m                Adjusted to get equivalent mass per unit width    


-- Initial conditions in reservoir
startP = 100000 -- Pa
startT = 348.4 -- K
startL = 4-0.005 -- m

-- Set filename to store movement data
outfile_name = "output.dat"
