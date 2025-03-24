-- =====================================================================
--
-- Example of the superellipsoid structured grid patch used to make a varietry of shapes.
-- 
-- Author: Flynn M. Hack
-- Date: 2023-06-22
--
-- The patch is formed by scaling points on a cube face so they 
-- lie on the superellipsoid surface. The general surface is defined by
-- the surface function:
--     1 = ((x / ax) ** (2/e2) + (y / ay) ** (2/e2)) ** (e2/e1) + (z / az) ** (2/e1)
--
-- Currently, the following simplified superellipsoid surface function has been implemented:
--     1 = (x) ** (2/e) + (y) ** (2/e) + z ** (2/e)
--
-- As e --> 0, the surface tends to a perfect cube with sharp edges.
--
-- To run these examples and produce VTK files:
--
--   $ e4shared --custom-script --script-file=super_ellipsoid-examples.lua
--
-- =====================================================================

-- Shape side length 2a

a_c = 1.0 -- m
e_shape1 = 0.1
e_shape2 = 1.0
e_shape3 = 2.0

niv = 51
njv = 51

require "super_ellipsoid_patch"

faces = {"east", "west", "north", "south", "top", "bottom"}

-- Shape 1 - Cube with rounded edges

for _,f in ipairs(faces) do
    
    super_ellipsoid = super_ellipsoid_patch:new{face = f, a = a_c, e = e_shape1}
   
    patch = super_ellipsoid:create()

    grid = StructuredGrid:new{psurface=patch, niv=niv, njv=njv}
    grid:write_to_vtk_file(string.format('rounded_cube_%s.vtk', f))  
end

-- Shape 2 - Sphere (Note: a sphere patch with a nicer distribution of cells already exists in eilmer)

for _,f in ipairs(faces) do
    
    super_ellipsoid = super_ellipsoid_patch:new{face = f, a = a_c, e = e_shape2}
   
    patch = super_ellipsoid:create()

    grid = StructuredGrid:new{psurface=patch, niv=niv, njv=njv}
    grid:write_to_vtk_file(string.format('sphere_%s.vtk', f))  
end

-- Shape 3 - Diamond

for _,f in ipairs(faces) do
    
    super_ellipsoid = super_ellipsoid_patch:new{face = f, a = a_c, e = e_shape3}
   
    patch = super_ellipsoid:create()

    grid = StructuredGrid:new{psurface=patch, niv=niv, njv=njv}
    grid:write_to_vtk_file(string.format('diamond_%s.vtk', f))  
end
