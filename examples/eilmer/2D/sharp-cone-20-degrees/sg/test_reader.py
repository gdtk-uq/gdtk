# test_reader.py
# Try out the Python service functions for reading grids and flow fields.
# $ python3 test_reader.py
# PJ 2023-04-07

from gdtk.geom.sgrid import StructuredGrid
from gdtk.flow.field import Field
from gdtk.flow.vtk_writer import write_vtk_structured_grid_file

g = StructuredGrid(gzfile="grid/t0000/cone20.grid.b0000.t0000.gz")

print("g.label=", g.label)
print("g.niv=", g.niv)
print("g.njv=", g.njv)
print("g.nkv=", g.nkv)
# print("g=", g)
g.write_to_gzip_file("test-grid.gz", format_version="1.1")

f = Field(gzfile="flow/t0004/cone20.flow.b0000.t0004.gz")
print("f=", f)
for i in range(f.nic):
    x = f.data['pos.x'][i,0]
    y = f.data['pos.y'][i,0]
    p = f.data['p'][i,0]
    print(f"i={i} pos=({x},{y}) pressure={p}")

write_vtk_structured_grid_file("b0.vts", g, f)

g1 = StructuredGrid(gzfile="grid/t0000/cone20.grid.b0001.t0000.gz")
f1 = Field(gzfile="flow/t0004/cone20.flow.b0001.t0004.gz")
write_vtk_structured_grid_file("b1.vts", g1, f1)
