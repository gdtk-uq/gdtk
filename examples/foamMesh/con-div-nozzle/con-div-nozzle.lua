-- Mesh for optimisation of convergent-divergent nozzle
-- Author: Ingo Jahn
-- last modified: 08/04/2017

axisymmetric = true
--#########################################
--# Create Geometry                     ###
--#########################################
--     +-----------+--\   rt  /--+--------+  R
--     |           |   \-+---/   |        |
--     |   A0      | A1  |  A2   |   A3   |
--     |           |     |       |        |
--     +-----------+-----+-------+--------+  r=0
--
--
--     a0---------a1--\        /--a3------a4
--     |           |   \-a2---/   |        |
--     |   A0      | A1  |   A2   |   A3   |  nr
--     |           |     |        |        |
--     b0---------b1----b2--------b3------b4
--     x0         x1    x2        x3      x4
--          n0        n1     n3      n3     
--
-- Bezier Curve Control Points (define as fractions of corner points)
--
--     a1---c0                      (fixed at R) 
-- 
--                c1                c1_rf (free to move in r and x)
--                                  
--                        c2---a2   (fixed at rt)
--       c0_xf   c1_xf    c2_xf

-- ###################
-- Input variables to parmetrically define nozzle
-- ###################
-- Geometric parameters
x0 = -0.2
x1 = -0.10
x3 = 0.1
x4 = 0.2
R = 0.1
Rc = 1e-6
x2 = 0.0  -- = (1-OP[0])*x1 + OP[0]*x3
Rt = 0.05  -- = (1-OP[1])*Rc + OP[1]*R
-- Bezier Curve Control points for a1a2 and a2a3
c0_xf = 0.2 --= OP[2]
c1_xf = 0.5 --= OP[3]
c1_rf = 0.5 --= OP[4]
c2_xf = 0.8 --= OP[5]
-- Bezier Curve Control points for a2a3
d0_xf = 0.2 --= OP[6]
d1_xf = 0.5 --= OP[7]
d1_rf = 0.5 --= OP[8]
d2_xf = 0.8 --= OP[9]

-- Define fixed points
a0 = Vector3:new{x=x0, y=R}
a1 = Vector3:new{x=x1, y=R}
a2 = Vector3:new{x=x2, y=Rt}
a3 = Vector3:new{x=x3, y=R}
a4 = Vector3:new{x=x4, y=R}

b0 = Vector3:new{x=x0, y=Rc}
b1 = Vector3:new{x=x1, y=Rc}
b2 = Vector3:new{x=x2, y=Rc}
b3 = Vector3:new{x=x3, y=Rc}
b4 = Vector3:new{x=x4, y=Rc}

-- define Bezier control points
c0 = Vector3:new{x=((1-c0_xf)*x1+c0_xf*x2), y=R}
c1 = Vector3:new{x=((1-c1_xf)*x1+c1_xf*x2), y=((1-c1_rf)*Rt+c1_rf*R)}
c2 = Vector3:new{x=((1-c2_xf)*x1+c2_xf*x2), y=Rt}
d0 = Vector3:new{x=((1-d0_xf)*x2+d0_xf*x3), y=Rt}
d1 = Vector3:new{x=((1-d1_xf)*x2+d1_xf*x3), y=((1-d1_rf)*Rt+c1_rf*R)}
d2 = Vector3:new{x=((1-d2_xf)*x2+d2_xf*x3), y=R}

-- create Bezier Curves
a1a2 = Bezier:new{points={a1,c0,c1,c2,a2}}
a2a3 = Bezier:new{points={a2,d0,d1,d2,a3}}

-- Define patch (which are parametric surfaces, no discretisation at this point.)
surf = {}
surf[0] = CoonsPatch:new{p00=b0, p10=b1, p11=a1, p01=a0}
surf[1] = CoonsPatch:new{north=a1a2, south=Line:new{p0=b1, p1=b2}, west=Line:new{p0=b1, p1=a1}, east=Line:new{p0=b2, p1=a2} }
surf[2] = CoonsPatch:new{north=a2a3, south=Line:new{p0=b2, p1=b3}, west=Line:new{p0=b2, p1=a2}, east=Line:new{p0=b3, p1=a3} }
surf[3] = CoonsPatch:new{p00=b3, p10=b4, p11=a4, p01=a3}

--     a0---------a1--\        /--a3------a4
--     |           |   \-a2---/   |        |
--     |   A0      | A1  |   A2   |   A3   |  nr
--     |           |     |        |        |
--     b0---------b1----b2--------b3------b4
--     x0         x1    x2        x3      x4
--          n0        n1     n3      n3     

-- Define 2D grid on patch, clustering can be added if desired
n0=20; n1=30; n2=30; n3=20
nr=20

cfr = RobertsFunction:new{end0=false, end1=true, beta=1.05}
cf0 = RobertsFunction:new{end0=false, end1=true, beta=1.12}
cf1 = RobertsFunction:new{end0=true, end1=false, beta=1.12}
 
--cfr = None --RobertsFunction:new{end0=true, end1=true, beta=1.05}
grid = {}
grid[0] = StructuredGrid:new{psurface=surf[0], niv=n0, njv=nr, 
                cfList={east=cfr,west=cfr,north=cf0,south=cf0} }
grid[1] = StructuredGrid:new{psurface=surf[1], niv=n1, njv=nr, 
              cfList={east=cfr,west=cfr} }
grid[2] = StructuredGrid:new{psurface=surf[2], niv=n2, njv=nr, 
                cfList={east=cfr,west=cfr} }
grid[3] = StructuredGrid:new{psurface=surf[3], niv=n3, njv=nr, 
                cfList={east=cfr,west=cfr,north=cf1,south=cf1} }

-- Define OpenFoam block (a "grid" with labels)
block = {}
block[0] = FoamBlock:new{grid=grid[0],
		     bndry_labels={west="i-00", north="w-00", south="s-00"}}
block[1] = FoamBlock:new{grid=grid[1], bndry_labels={north="w-00", south="s-00"}}
block[2] = FoamBlock:new{grid=grid[2], bndry_labels={north="w-00", south="s-00"}}
block[3] = FoamBlock:new{grid=grid[3], bndry_labels={north="w-00", south="s-00", east="o-00"}}


