# test_reader.py
# Try out the Python service functions for reading grids and flow fields.
# $ python3 test_reader.py
# PJ 2023-04-07

from gdtk.geom.sgrid import StructuredGrid
from gdtk.flow.field import Field

g = StructuredGrid(gzfile="grid/t0000/cone20-nff.grid.b0000.t0000.gz")

print("g.label=", g.label)
print("g.niv=", g.niv)
print("g.njv=", g.njv)
print("g.nkv=", g.nkv)
# print("g=", g)
g.write_to_gzip_file("test-grid.gz", format_version="1.1")

f = Field(ziparchive="CellData/field/t0000/cone20-nff.field.b0000.t0000.zip")
print("f=", f)
for i in range(f.nic):
    x = f.data['pos.x'][i,0]
    y = f.data['pos.y'][i,0]
    p = f.data['p'][i,0]
    print(f"i={i} pos=({x},{y}) pressure={p}")

