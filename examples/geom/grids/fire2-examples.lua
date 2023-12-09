-- Example structured grids on a pre-shock-fitted FireII  domain.
-- Author: Rowan J. Gollan
-- Date: 2022-01-13
--
-- To generate this domain, I first ran the simulation built by Nick G.
-- found in: examples/eilmer/2D/fireII/gupta1989.
-- Then I extracted the shock shape as a series of (x, y) coordinates.
-- The other boundaries --- body, symmetry and outflow --- are taken
-- directly from Nick's example.
--
-- To run these examples and produce VTK files:
--
-- $ e4shared --custom-script --script-file=fire2-examples.lua
--


niv = 15
njv = 50

Ri = 0.9347             
ri = 0.0102
A  = 0.3358
thetai = math.asin((A-ri)/(Ri-ri))
L = 0.05
diffo = 0.07
Ro = Ri + diffo
ro = ri + diffo
psio = math.rad(5.0)

oi = Vector3:new{x=Ri, y=0.0}
ai = oi + Ri*Vector3:new{x=-1.0, y=0.0}
bi = oi + Ri*Vector3:new{x=-math.cos(thetai), y=math.sin(thetai)}
pi = oi + (Ri-ri)*Vector3:new{x=-math.cos(thetai), y=math.sin(thetai)}
ci = pi + ri*Vector3:new{x=0.0, y=1.0}
di = ci + L*Vector3:new{x=1.0, y=0.0}

aibi = Arc:new{p0=ai, p1=bi, centre=oi}
bici = Arc:new{p0=bi, p1=ci, centre=pi}
cidi = Line:new{p0=ci, p1=di}
aidi = Polyline:new{segments={aibi, bici, cidi}}

aodo = Spline2:new{filename='shock-shape.dat'}
aoai = Line:new{p0=aodo(0.0), p1=ai}
dodi = Line:new{p0=aodo(1.0), p1=di}

-- Example 1: Coons-patch grid.
coonsPatch = CoonsPatch:new{north=dodi, south=aoai, east=aidi, west=aodo}
coonsGrid = StructuredGrid:new{psurface=coonsPatch, niv=niv, njv=njv}
coonsGrid:write_to_vtk_file('fire2-coons-patch-grid.vtk')

-- Example 2: AO-patch grid
AOPatch = AOPatch:new{north=dodi, south=aoai, east=aidi, west=aodo}
AOGrid = StructuredGrid:new{psurface=AOPatch, niv=niv, njv=njv}
AOGrid:write_to_vtk_file('fire2-AO-patch-grid.vtk')

-- Example 3: Channel patch
channelPatch = ChannelPatch:new{south=aodo, north=aidi}
channelGrid = StructuredGrid:new{psurface=channelPatch, niv=njv, njv=niv}
channelGrid:write_to_vtk_file('fire2-channel-patch-grid.vtk')

-- Example 3: Ruled surface
ruledSurface = RuledSurface:new{edge0=aodo, edge1=aidi, ruled_direction='r'}
ruledGrid = StructuredGrid:new{psurface=ruledSurface, niv=niv, njv=njv}
ruledGrid:write_to_vtk_file('fire2-ruled-surface-grid.vtk')

-- Example 4: Control-point patch
-- Use ChannelPatch as means to generate control pointe because these will be orthogonal to surface.
ncpi = 4
ncpj = 13

ctrlPtPatch = ControlPointPatch:new{north=dodi, south=aoai, east=aidi, west=aodo,
                                    ncpi=ncpi, ncpj=ncpj, guide_patch="channel-e2w"}
ctrlPtPatch:writeCtrlPtsAsVtkXml('fire2-ctrl-pts')
ctrlPtGrid = StructuredGrid:new{psurface=ctrlPtPatch, niv=niv, njv=njv}
ctrlPtGrid:write_to_vtk_file('fire2-ctrl-pt-patch-grid.vtk')
