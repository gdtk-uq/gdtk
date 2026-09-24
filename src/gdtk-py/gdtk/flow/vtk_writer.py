# vtk_writer.py
# Write VTK-format files using the data from grid and flow-field objects.
# Peter J.
# 2023-04-09 adapted from the lorikeet-postprocessing code.
#

import math
from gdtk.geom.sgrid import StructuredGrid
from gdtk.flow.field import Field


def write_vtk_structured_grid_file(vtkFile, grid, flow):
    """
    Combine the finite-volume grid and flow data for a structured-grid block
    into a VTK StructuredGrid file.

    Input:
      vtkFile: name of the resulting VTK file, typically xxxx.vts
      grid:    StructuredGrid object
      flow:    flow/Field object

    It is the caller's responsibility to have matching grid and flow-field objects.
    """
    assert isinstance(grid, StructuredGrid), "grid is expected to be a StructuredGrid object"
    assert isinstance(flow, Field), "flowfield is expected to a flow/Field object"
    assert grid.niv-1 == flow.nic, "i-direction mismatch in cell and vertex numbers"
    assert grid.njv-1 == flow.njc, "j-direction mismatch in cell and vertex numbers"
    if grid.dimensions == 2:
        assert grid.nkv == 1 and flow.nkc == 1, "2D grid needs to have niv==1 and nkc==1"
    else:
        assert grid.nkv-1 == flow.nkc, "3D grid mismatch in cell and vertex numbers"
    #
    fp = open(vtkFile, mode='w')
    fp.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian">\n')
    if grid.dimensions == 2:
        fp.write('<StructuredGrid WholeExtent="%d %d %d %d 0 0">\n' % (0, grid.niv-1, 0, grid.njv-1))
        fp.write('<Piece Extent="%d %d %d %d 0 0">\n' % (0, grid.niv-1, 0, grid.njv-1))
    else:
        fp.write('<StructuredGrid WholeExtent="%d %d %d %d %d %d">\n' % (0, grid.niv-1, 0, grid.njv-1, 0, grid.nkv-1))
        fp.write('<Piece Extent="%d %d %d %d %d %d">\n' % (0, grid.niv-1, 0, grid.njv-1, 0, grid.nkv-1))
    #
    fp.write('<CellData>\n')
    for var in flow.variables:
        fp.write('<DataArray Name="%s" type="Float64" NumberOfComponents="1" format="ascii">\n' % (var,))
        # Our arrays of data are indexed as [i, j, k]
        # VTK format has k as outer loop and i as inner loop.
        data = flow.data[var].transpose().flatten()
        for item in data: fp.write('%g\n' % (item,))
        fp.write('</DataArray>\n')
    if "vel.x" in flow.variables:
        vxdata = flow.data['vel.x'].transpose().flatten()
        vydata = flow.data['vel.y'].transpose().flatten()
        vzdata = flow.data['vel.z'].transpose().flatten()
        if "a" in flow.variables:
            adata = flow.data['a'].transpose().flatten()
            fp.write('<DataArray Name="Mach" type="Float64" NumberOfComponents="1" format="ascii">\n');
            for a,vx,vy,vz in zip(adata, vxdata, vydata, vzdata):
                fp.write('%g\n' % (math.sqrt(vx*vx+vy*vy+vz*vz)/a))
            fp.write('</DataArray>\n')
        # Write the velocity vector.
        fp.write('<DataArray Name="vel.vector" type="Float64" NumberOfComponents="3" format="ascii">\n');
        for velxyz in zip(vxdata, vydata, vzdata): fp.write('%g %g %g\n' % velxyz)
        fp.write('</DataArray>\n')
    if "B.x" in flow.variables:
        Bxdata = flow.data['B.x'].transpose().flatten()
        Bydata = flow.data['B.y'].transpose().flatten()
        Bzdata = flow.data['B.z'].transpose().flatten()
        fp.write('<DataArray Name="B.vector" type="Float64" NumberOfComponents="3" format="ascii">\n');
        for Bxyz in zip(Bxdata, Bydata, Bzdata): fp.write('%g %g %g\n' % Bxyz)
        fp.write('</DataArray>\n')
    fp.write('</CellData>\n')
    #
    fp.write('<Points>\n')
    fp.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
    xs = grid.vertices.x.transpose().flatten()
    ys = grid.vertices.y.transpose().flatten()
    zs = grid.vertices.z.transpose().flatten()
    for xyz in zip(xs, ys, zs):
        fp.write('%g %g %g\n' % xyz)
    fp.write('</DataArray>\n')
    fp.write('</Points>\n')
    #
    fp.write('</Piece>\n')
    fp.write('</StructuredGrid>\n')
    fp.write('</VTKFile>\n')
    fp.close()
    return

