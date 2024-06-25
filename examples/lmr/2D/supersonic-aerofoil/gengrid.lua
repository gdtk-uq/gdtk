
function area_of_a_triangle(A, B, C)
    -- The "shoelace formula" from https://en.wikipedia.org/wiki/Area_of_a_triangle
    thing = A.x*B.y - A.x*C.y + B.x*C.y - B.x*A.y + C.x*A.y - C.x*B.y
    area = 0.5*math.abs(thing)
    return area
end

function centroid_of_a_triangle(A, B, C)
    -- https://en.wikipedia.org/wiki/List_of_centroids
    ox = (A.x + B.x + C.x)/3.0
    oy = (A.y + B.y + C.y)/3.0
    return Vector3:new{x=ox, y=oy}
end

f = 5e-3
t = 10e-3

H = 100e-3
l = 100e-3
h = 5e-3
tailfac = 1.0/5.0
xtail = l*(1-tailfac)
ytail = h/2.0*tailfac/(1-0.5) -- the half here is a in your derivation

base_angle = math.tan((h/2.0)/(l*(1-0.5)))
trim_angle = math.rad(0.0)
r = math.sqrt(ytail*ytail + tailfac*l*tailfac*l)
xtip = xtail + r*math.cos(base_angle-trim_angle)
ytip = ytail + r*math.sin(base_angle-trim_angle)

a = Vector3:new{x=-f, y=0.0}
b = Vector3:new{x=0.0, y=0.0}
c = Vector3:new{x=l/2.0, y=h/2.0}
q = Vector3:new{x=xtail, y=ytail}
d = Vector3:new{x=xtip, y=ytip}
e = Vector3:new{x=xtip+r+t, y=ytip}

A = Vector3:new{x=-f, y=H/2.0}
B = Vector3:new{x=0.0, y=H/2.0}
C = Vector3:new{x=l/2.0, y=H/2.0}
Q = Vector3:new{x=xtail, y=H/2.0}
D = Vector3:new{x=l, y=H/2.0}
E = Vector3:new{x=l+r+t, y=H/2.0}

ab = Line:new{p0=a,p1=b}
bc = Line:new{p0=b,p1=c}
cq = Line:new{p0=c,p1=q}
qd = Line:new{p0=q,p1=d}
de = Line:new{p0=d,p1=e}

AB = Line:new{p0=A,p1=B}
BC = Line:new{p0=B,p1=C}
CQ = Line:new{p0=C,p1=Q}
QD = Line:new{p0=Q,p1=D}
DE = Line:new{p0=D,p1=E}

aA = Line:new{p0=a,p1=A}
bB = Line:new{p0=b,p1=B}
cC = Line:new{p0=c,p1=C}
dD = Line:new{p0=d,p1=D}
qQ = Line:new{p0=q,p1=Q}
eE = Line:new{p0=e,p1=E}

patch0 = CoonsPatch:new{north=AB, east=bB, south=ab, west=aA}
patch1 = CoonsPatch:new{north=BC, east=cC, south=bc, west=bB}
patch2 = CoonsPatch:new{north=CQ, east=qQ, south=cq, west=cC}
patch3 = CoonsPatch:new{north=QD, east=dD, south=qd, west=qQ}
patch4 = CoonsPatch:new{north=DE, east=eE, south=de, west=dD}

-- Mirror grid
h = -h
H = -H
ytail = -ytail
a_ = Vector3:new{x=-f, y=0.0}
b_ = Vector3:new{x=0.0, y=0.0}
c_ = Vector3:new{x=l/2.0, y=h/2.0}
q_ = Vector3:new{x=xtail, y=ytail}
d_ = Vector3:new{x=xtip, y=ytip}
e_ = Vector3:new{x=xtip+r+t, y=ytip}

A_ = Vector3:new{x=-f, y=H/2.0}
B_ = Vector3:new{x=0.0, y=H/2.0}
C_ = Vector3:new{x=l/2.0, y=H/2.0}
Q_ = Vector3:new{x=xtail, y=H/2.0}
D_ = Vector3:new{x=l, y=H/2.0}
E_ = Vector3:new{x=l+r+t, y=H/2.0}

ab_ = Line:new{p0=a_,p1=b_}
bc_ = Line:new{p0=b_,p1=c_}
cq_ = Line:new{p0=c_,p1=q_}
qd_ = Line:new{p0=q_,p1=d_}
de_ = Line:new{p0=d_,p1=e_}

AB_ = Line:new{p0=A_,p1=B_}
BC_ = Line:new{p0=B_,p1=C_}
CQ_ = Line:new{p0=C_,p1=Q_}
QD_ = Line:new{p0=Q_,p1=D_}
DE_ = Line:new{p0=D_,p1=E_}

Aa_ = Line:new{p0=A_,p1=a_}
Bb_ = Line:new{p0=B_,p1=b_}
Cc_ = Line:new{p0=C_,p1=c_}
Qq_ = Line:new{p0=Q_,p1=q_}
Dd_ = Line:new{p0=D_,p1=d_}
Ee_ = Line:new{p0=E_,p1=e_}

-- Interlude, compute the centroid of the aerofoil shape
centrepoint = Vector3:new{x=l/2.0, y=0.0}
alpha   = {A=b, B=c, C=c_}
beta    = {A=c, B=q, C=centrepoint}
charlie = {A=centrepoint, B=q, C=q_}
delta   = {A=centrepoint, B=q_, C=c_}
echo    = {A=q, B=d, C=q_}

alpha_area   = area_of_a_triangle(alpha.A, alpha.B, alpha.C)
beta_area    = area_of_a_triangle(beta.A, beta.B, beta.C)
charlie_area = area_of_a_triangle(charlie.A, charlie.B, charlie.C)
delta_area   = area_of_a_triangle(delta.A, delta.B, delta.C)
echo_area    = area_of_a_triangle(echo.A, echo.B, echo.C)

alpha_centroid   = centroid_of_a_triangle(alpha.A,   alpha.B,   alpha.C)
beta_centroid    = centroid_of_a_triangle(beta.A,    beta.B,    beta.C)
charlie_centroid = centroid_of_a_triangle(charlie.A, charlie.B, charlie.C)
delta_centroid   = centroid_of_a_triangle(delta.A,   delta.B,   delta.C)
echo_centroid    = centroid_of_a_triangle(echo.A,    echo.B,    echo.C)

aerofoil_area = alpha_area + beta_area + charlie_area + delta_area + echo_area

aerofoil_centroid = (alpha_centroid   * alpha_area   +
                     beta_centroid    * beta_area    +
                     charlie_centroid * charlie_area +
                     delta_centroid   * delta_area   +
                     echo_centroid    * echo_area  )/aerofoil_area

print("Centre of Mass: ", aerofoil_centroid)

patch5 = CoonsPatch:new{north=ab_, east=Bb_, south=AB_, west=Aa_}
patch6 = CoonsPatch:new{north=bc_, east=Cc_, south=BC_, west=Bb_}
patch7 = CoonsPatch:new{north=cq_, east=Qq_, south=CQ_, west=Cc_}
patch8 = CoonsPatch:new{north=qd_, east=Dd_, south=QD_, west=Qq_}
patch9 = CoonsPatch:new{north=de_, east=Ee_, south=DE_, west=Dd_}

factor = 0.5
niv = 128*factor; njv=64*factor
a = 0.0002/factor
if factor==0.5 then
    a = 0.0040
end

xlength = f + l + r
--cfx = RobertsFunction:new{end0=false,end1=true,beta=1.5}
cflist = {north=cfx, east=GeometricFunction:new{a=a, r=1.2, N=njv},
          south=cfx, west=GeometricFunction:new{a=a, r=1.2, N=njv}}
grid0 = StructuredGrid:new{psurface=patch0, niv=math.floor(niv*f/xlength), njv=njv, cfList=cflist}
grid1 = StructuredGrid:new{psurface=patch1, niv=math.floor(niv*l/2.0/xlength), njv=njv, cfList=cflist}
grid2 = StructuredGrid:new{psurface=patch2, niv=math.floor(niv*l*(1.0-tailfac-0.5)/xlength), njv=njv, cfList=cflist}
grid3 = StructuredGrid:new{psurface=patch3, niv=math.floor(niv*l*tailfac/xlength), njv=njv, cfList=cflist}
grid4 = StructuredGrid:new{psurface=patch4, niv=math.floor(niv*(r+t)/xlength), njv=njv, cfList=cflist}

cflist2 = {north=cfx, east=GeometricFunction:new{a=a, r=1.2, N=njv, reverse=true},
           south=cfx, west=GeometricFunction:new{a=a, r=1.2, N=njv, reverse=true}}

grid5 = StructuredGrid:new{psurface=patch5, niv=math.floor(niv*f/xlength), njv=njv, cfList=cflist2}
grid6 = StructuredGrid:new{psurface=patch6, niv=math.floor(niv*l/2.0/xlength), njv=njv, cfList=cflist2}
grid7 = StructuredGrid:new{psurface=patch7, niv=math.floor(niv*l*(1.0-tailfac-0.5)/xlength), njv=njv, cfList=cflist2}
grid8 = StructuredGrid:new{psurface=patch8, niv=math.floor(niv*l*tailfac/xlength), njv=njv, cfList=cflist2}
grid9 = StructuredGrid:new{psurface=patch9, niv=math.floor(niv*(r+t)/xlength), njv=njv, cfList=cflist2}

ugrid0 = UnstructuredGrid:new{sgrid=grid0}
ugrid1 = UnstructuredGrid:new{sgrid=grid1}
ugrid2 = UnstructuredGrid:new{sgrid=grid2}
ugrid3 = UnstructuredGrid:new{sgrid=grid3}
ugrid4 = UnstructuredGrid:new{sgrid=grid4}

ugrid5 = UnstructuredGrid:new{sgrid=grid5}
ugrid6 = UnstructuredGrid:new{sgrid=grid6}
ugrid7 = UnstructuredGrid:new{sgrid=grid7}
ugrid8 = UnstructuredGrid:new{sgrid=grid8}
ugrid9 = UnstructuredGrid:new{sgrid=grid9}

-- New Style
WEST = 0
EAST = 1
SOUTH = 2
NORTH = 3

-- Old Style
--NORTH = 0
--EAST = 1
--SOUTH = 2
--WEST = 3

ugrid0:set_boundaryset_tag(WEST,  "inflow");  ugrid0:set_boundaryset_tag(NORTH, "upper")
ugrid1:set_boundaryset_tag(SOUTH, "wall");    ugrid1:set_boundaryset_tag(NORTH, "upper")
ugrid2:set_boundaryset_tag(SOUTH, "wall");    ugrid2:set_boundaryset_tag(NORTH, "upper")
ugrid3:set_boundaryset_tag(SOUTH, "wall");    ugrid3:set_boundaryset_tag(NORTH, "upper")
ugrid4:set_boundaryset_tag(EAST,  "outflow"); ugrid0:set_boundaryset_tag(NORTH, "upper")

ugrid5:set_boundaryset_tag(WEST, "inflow");  ugrid5:set_boundaryset_tag(SOUTH, "lower")
ugrid6:set_boundaryset_tag(NORTH, "wall");    ugrid6:set_boundaryset_tag(SOUTH, "lower")
ugrid7:set_boundaryset_tag(NORTH, "wall");    ugrid7:set_boundaryset_tag(SOUTH, "lower")
ugrid8:set_boundaryset_tag(NORTH, "wall");    ugrid8:set_boundaryset_tag(SOUTH, "lower")
ugrid9:set_boundaryset_tag(EAST, "outflow"); ugrid9:set_boundaryset_tag(SOUTH, "lower")

ugrid0:joinGrid(ugrid1)
ugrid0:joinGrid(ugrid2)
ugrid0:joinGrid(ugrid3)
ugrid0:joinGrid(ugrid4)
ugrid0:joinGrid(ugrid5)
ugrid0:joinGrid(ugrid6)
ugrid0:joinGrid(ugrid7)
ugrid0:joinGrid(ugrid8)
ugrid0:joinGrid(ugrid9)

ugrid0:write_to_su2_file('grid.su2')
