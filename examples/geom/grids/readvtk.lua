-- Read in a structured grid using the Lua interface. Unfortunately the unstructured vtk
-- reader is not yet working.
-- 
-- @author: Nick Gibbons


H = 0.5
L = 1.0

a = Vector3:new{x=0.0, y=H};   b = Vector3:new{x=L/2, y=H};
d = Vector3:new{x=0.0, y=0.0}; e = Vector3:new{x=L/2, y=0.0};

patch  = CoonsPatch:new{p00=d, p10=e, p11=b, p01=a}

sgrid = StructuredGrid:new{psurface=patch, niv=64, njv=64}
sgrid:write_to_vtk_file('sgrid.vtk')

ugrid = UnstructuredGrid:new{sgrid=sgrid}
ugrid:write_to_vtk_file('ugrid.vtk')

read_sgrid = StructuredGrid:new{filename="sgrid.vtk", fmt="vtk"}
print(read_sgrid)

--read_ugrid = UnstructuredGrid:new{filename="ugrid.vtk", fmt="vtktext"}
--print(read_ugrid)
