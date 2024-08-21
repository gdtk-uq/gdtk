-- Author: Rowan J. Gollan & Peter J.
-- Date: 2017-01-04 -- 2019-05-21, 2024-08-21 lmr(5) port
--
-- This simulation mimics the evolution of the gas in
-- the right-half of the classic Sod shock tube problem,
-- and exercises the moving-grid capability of Eilmer.
-- The piston has a constant velocity, high enough to drive
-- a shock wave through the gas in front of it.
--
print("Fast-moving piston with constant velocity.")
config.dimensions = 2
-- Geometry, grid and block setup.
L = 0.5; H = 0.1
-- Gas region that drives piston.
patch0 = CoonsPatch:new{p00=Vector3:new{x=0, y=0},
                        p10=Vector3:new{x=L, y=0},
                        p11=Vector3:new{x=L, y=H},
                        p01=Vector3:new{x=0, y=H}}
grid0 = StructuredGrid:new{psurface=patch0, niv=51, njv=3}
registerFluidGrid{grid=grid0, fsTag='initial', bcTags={west='piston_face'}}
