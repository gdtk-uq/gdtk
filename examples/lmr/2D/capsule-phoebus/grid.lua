-- Grid generation for the phoebus capsule geometry from ref. [1].
--
-- References:
-- [1] D. Bianchi et al. “Numerical Analysis and Wind Tunnel Validation of 
--     Low-Temperature Ablators undergoing Shape Change”. In: International 
--     Journal of Heat and Mass Transfer 177 (2021), p. 121430. issn: 0017-9310.
--     doi: 10.1016/j.ijheatmasstransfer.2021.121430.
--
-- @author: Reece B. Otto (2024-11-19)

config.dimensions = 2

-- geometric parameters
Rn = 20.0e-3          -- nose radius [m]
beta = math.rad(45.0) -- cone angle [rad]
Rs = 1.57e-3          -- shoulder radius [m]
Rb = 20.0e-3          -- base radius [m]
L_aft = Rs            -- length of artificial wall downstream of shoulder [m]
Ri = Rn*1.2           -- initial radius of inflow boundary [m]
Rf = Rn*1.5           -- final radius of inflow boundary  [m]

-- points
A = Vector3:new{x=0.0, y=0.0}
B = Vector3:new{x=Rn*(1.0-math.cos(beta)), y=Rn*math.sin(beta)}
Dy = Rb
Cy = Dy + Rs*(math.cos(beta) - 1.0)
C = Vector3:new{x=B.x + (Cy-B.y)/math.tan(beta), y=Cy}
D = Vector3:new{x=C.x + Rs*math.sin(beta), y=Rb}
E = Vector3:new{x=D.x + L_aft, y=D.y}
theta = math.atan(E.y/(Rn - E.x)) -- angle of outflow boundary
F = Vector3:new{x=Rn - Rf*math.cos(theta), y=Rf*math.sin(theta)}
G = Vector3:new{x=Rn-Ri, y=0.0}
centre_n = Vector3:new{x=Rn, y=0.0}
centre_s = Vector3:new{x=D.x, y=D.y-Rs}

-- paths
AB = Arc:new{p0=A, p1=B, centre=centre_n}
BC = Line:new{p0=B, p1=C}
CD = Arc:new{p0=C, p1=D, centre=centre_s}
DE = Line:new{p0=D, p1=E}
AE = Polyline:new{segments={AB,BC,CD,DE}}
AE = ArcLengthParameterizedPath:new{underlying_path=AE}
FE = Line:new{p0=F, p1=E}
GA = Line:new{p0=G, p1=A}
GF1 = Vector3:new{x=G.x, y=B.y}
GF2 = Vector3:new{x=B.x, y=0.8*F.y}
GF = Bezier:new{points={G, GF1, GF2, F}}
GF = ArcLengthParameterizedPath:new{underlying_path=GF}

-- surfaces
ncpi = 4
ncpj = 13
surf = ControlPointPatch:new{south=GA, north=FE, west=GF, east=AE,
                             ncpi=ncpi, ncpj=ncpj, guide_patch="channel_e2w"}
surf:setCtrlPt(1, 11, {x=0.011, y=0.0255})
surf:setCtrlPt(2, 11, {x=0.012, y=0.022})
--surf:writeCtrlPtsAsVtkXml("cntrl_pts")

-- grid vertex counts
grid_level = 1 -- grid refinement factor (set as 4-5 for accurate heat flux)
n_refine = 2^(grid_level/2)
niv = math.ceil(100*n_refine + 1)
njv = math.ceil(20*n_refine + 1)
print("Total cells = ", (niv-1)*(njv-1))

-- clustering
r = 1.2
h_wall = 1.0e-4 -- cell height at wall (set to ~ 1.0e-6 for adequate BL resolution)
d_south = (A-G):abs()
cluster_south = GeometricFunction:new{a=h_wall/d_south, r=r, N=njv, reverse=true}
d_north = (F-E):abs()
cluster_north = GeometricFunction:new{a=h_wall/d_north, r=r, N=njv, reverse=true}
cf_list = {south=cluster_south, north=cluster_north}

-- make grid
sgrid = StructuredGrid:new{psurface=surf, njv=niv, niv=njv, cfList=cf_list}
grid = registerFluidGridArray{
  grid=sgrid,
  fsTag="initial",
  shock_fitting = true,
  bcTags={east="wall", west="inflow", south="symmetry", north="outflow"},
  nib=2,
  njb=4
}
identifyGridConnections()
