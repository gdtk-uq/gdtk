print("Piston in tube, 2 FluidBlocks.")
-- Authors: Rowan G., FabZ, Ingo J., Peter J. 
-- Date: 2017-01-04 -- 2019-05-18, 2024-08-22 lmr(5) port
--
-- Dimensions match verification case 2 in the 16th AFMC paper (2007)
-- "Development of Casbar: a Two-phase Flow Code for the Interior Ballistics
-- Problem", R.J Gollan, I.A. Johnston, B.T. O'Flaherty and P.A. Jacobs
--
config.dimensions = 2
dofile('sim-config.lua')
-- Geometry, grid and block setup.
-- Gas region that drives piston.
driver_patch = CoonsPatch:new{p00=Vector3:new{x=0,  y=0},
                              p10=Vector3:new{x=L1, y=0},
                              p11=Vector3:new{x=L1, y=H},
                              p01=Vector3:new{x=0,  y=H}}
-- Gas region that is compressed by piston.
driven_patch = CoonsPatch:new{p00=Vector3:new{x=L2, y=0},
                              p10=Vector3:new{x=L3, y=0},
                              p11=Vector3:new{x=L3, y=H},
                              p01=Vector3:new{x=L2, y=H}}
nxc = 100; nyc = 2
grid0 = StructuredGrid:new{psurface=driver_patch, niv=nxc+1, njv=nyc+1}
grid1 = StructuredGrid:new{psurface=driven_patch, niv=nxc+1, njv=nyc+1}
registerFluidGrid{grid=grid0, fsTag='initialLeft', bcTags={east='piston_upstream_face'}}
registerFluidGrid{grid=grid1, fsTag='initialRight', bcTags={west='piston_downstream_face'}}
--
-- Note that there are no direct connections between the grids
-- because the piston is in between them.
-- Ignore that warning on finding no inter-block connections.
