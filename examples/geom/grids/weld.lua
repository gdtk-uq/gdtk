-- Join initially structured grids and output to su2. This can be handy if you want
-- to make a grid using Eilmer's native grid generator and then partition it using metis.
-- 
-- @author: Nick Gibbons


H = 0.5
L = 1.0

a = Vector3:new{x=0.0, y=H};   b = Vector3:new{x=L/2, y=H};   c = Vector3:new{x=L, y=H};
d = Vector3:new{x=0.0, y=0.0}; e = Vector3:new{x=L/2, y=0.0}; f = Vector3:new{x=L, y=0.0};


left_patch  = CoonsPatch:new{p00=d, p10=e, p11=b, p01=a}
right_patch = CoonsPatch:new{p00=e, p10=f, p11=c, p01=b}

left_grid = StructuredGrid:new{psurface=left_patch, niv=64, njv=64}
right_grid = StructuredGrid:new{psurface=right_patch, niv=64, njv=64}

left_ugrid = UnstructuredGrid:new{sgrid=left_grid}
right_ugrid = UnstructuredGrid:new{sgrid=right_grid}

left_ugrid:joinGrid(right_ugrid)
left_ugrid:write_to_vtk_file('welded.vtk')
