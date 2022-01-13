-- Example structured grids on a duct-like shape
-- Author: Rowan J. Gollan
-- Date: 2022-01-13
--
-- This domain is a facsimile of the duct-like domain
-- in Eiseman's paper "A control point form of algebraic grid generation."
-- The particular boundaries have been taken from Carrie Xie's
-- example in the technical note on Control Point Form grid generation.
--
-- To run these examples and produce VTK files:
--
--   $ e4shared --custom-script --script-file=duct-examples.lua
--

niv = 41
njv = 41

--x coordinates of control points on east boundary (5x5 control points)
L0 = 12
L1 = 11.8
L2 = 10
L3 = 8.2
L4 = 8

--y coordinates of control points on north boundary (5x5 control points)
H0 = 3.5
H1 = 3.7
H2 = 5.25
H3 = 6.8
H4 = 7

N = 5
M = 5

L={L0,L1,L2,L3,L4}
H={H0,H1,H2,H3,H4}

-- To ensure uniform distribution of coordinate curves,
-- the control points adjacent to boundaries have increments of half unit spacing.
-- x = C*(L/(N-2); y = D*(H/(M-2))
-- The function coeff(index,N) is used to compute coefficients C and D
function coeff(index,N)
   if index == 1 then return 0
   elseif index == 2 then return 0.5
   elseif index == N then return N-2
   else return index-1.5
   end
end

-- unit spacing of each horizontal line
usx={}
for j=1,M do
   usx[j] = L[j]/(N-2)
end

-- unit spacing of each vertical line
usy={}
for i=1,N do
   usy[i] = H[i]/(M-2)
end

-- Compute and store coefficients in vertical direction
cj={}
for j = 1,M do
   cj[j]=coeff(j,M)
end

-- Locate each control point
ctrl_pts = {}
for i=1,N do
   ctrl_pts[i] = {}
   -- Compute coefficient in horizontal direction
   ci = coeff(i,N)
   for j=1,M do
      ctrl_pts[i][j]=Vector3:new{x=ci*usx[j],y=cj[j]*usy[i]}
   end
end

-- west straight line boundary
west = Line:new{p0=ctrl_pts[1][1],p1=ctrl_pts[1][M]}

-- north Bezier boundary
n0 = ctrl_pts[1][M]
n1 = Vector3:new{x=ctrl_pts[2][M-1].x,y=(M-2.25)*usy[2]}
n2 = ctrl_pts[math.ceil(N/2)][M]
n3 = Vector3:new{x=ctrl_pts[N-1][M].x,y=(M-1.8)*usy[N-1]}
n4 = ctrl_pts[N][M]
north = Bezier:new{points={n0,n1,n2,n3,n4}}

-- east Bezier boundary
e0 = ctrl_pts[N][1]
e1 = Vector3:new{x=(N-1.95)*usx[2],y=0.75*usy[N]}
e2 = ctrl_pts[N][math.ceil(M/2)]
e3 = Vector3:new{x=(N-2.1)*usx[M-1],y=(N-2.75)*usy[N]}
e4 = ctrl_pts[N][M]
east = Bezier:new{points={e0,e1,e2,e3,e4}}

--south straight line boundary
south =  Line:new{p0=ctrl_pts[1][1],p1=ctrl_pts[N][1]}

-- Example 1: Coons patch
coonsPatch = CoonsPatch:new{north=north, south=south, east=east, west=west}
coonsGrid = StructuredGrid:new{psurface=coonsPatch, niv=niv, njv=njv}
coonsGrid:write_to_vtk_file('duct-coons-patch-grid.vtk')

-- Example 2: AO patch
aoPatch = AOPatch:new{north=north, south=south, east=east, west=west}
aoGrid = StructuredGrid:new{psurface=aoPatch, niv=niv, njv=njv}
aoGrid:write_to_vtk_file('duct-AO-patch-grid.vtk')

-- Example 3: Control point patch grid
ctrlPtPatch = ControlPointPatch:new{north=north, south=south, east=east, west=west,
                                    control_points=ctrl_pts}
ctrlPtGrid = StructuredGrid:new{psurface=ctrlPtPatch, niv=niv, njv=njv}
ctrlPtGrid:write_to_vtk_file('duct-ctrl-pt-patch-grid.vtk')

dofile('write_ctrl_points.lua')
writeCtrlPts2VTK(ctrl_pts, 'duct-ctrl-pts.vtk')



