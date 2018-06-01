-- Script for the generation a NACA00xx aerofoil mesh in 3-D. 
-- This example uses a mix of advanced mesh generation tools, 
-- so don't be to alarmed if it initially looks intimidating.
-- 
-- Make sure you attempt the other meshing tutorials first.
--
-- Author: Ingo Jahn
-- last modified: 18/05/2018

--#########################################
--# Create Geometry                     ###
--#########################################
--
--            ------f3--------N-------t0---N--t1
--          /       |                 |        |
--         N        |        blk1     |  blk0  |
--       /  blk2   -a3----\           /        E
--      /         /XXXXXXXX--------- /     S   |
--     f2-------a2XXX NACA FOIL XXXXa0--------a4
--      \         \XXXXXXXX--------- \     N   |
--       \  blk3   -a1----/           \        E
--        N         |        blk4     |  blk5  |
--          \       |                 |        |
--            ------f1--------N-------b0---S--b1


-- ###################
-- Input variables to parmetrically define the mesh
-- ###################
-- Geometric parameters
turbulenceModel = 'S-A'
thickness = 0.24  -- thickness (in percent) of NACA00XX aerofoil
c = 0.2 -- (m) chord length of the wing
S = 1.0 -- (m) span of the wing
L = 3*c -- (m) distance to far-field boundary
frac = 0.4 -- fraction to define position of point a1 and a3 along foil

-- ###################
-- Define Lines that create block boundaries
-- ###################

-- we use a function to decribe the outline of the aerofoil
function foil(t)
   -- function that return x/y positon along profiles as a function
   -- of paramter t, which start at t=0 at bottom rear, t=0.5 at leading
   -- edge and and finishes as t=1 at top rear.
   -- calculate y from NACA polynominal with last value adjusted to close trailing edge.
   LE=0.005
   if t < 0.5-LE then -- do lower edge
      x = 1.-2*t
      y = -5*thickness * ( 0.2969*x^0.5 - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1036*x^4)
   elseif t > 0.5+LE then -- do upper edge
      x = (t-0.5)*2
      y = 5*thickness * ( 0.2969*x^0.5 - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1036*x^4)
   else
      -- to get good clustering at the leading edge we can approximate the tip by a circle 
      r = 1.1019 *thickness^2
      xL = 1.-2*(0.5-LE)
      yL = -5*thickness * ( 0.2969*xL^0.5 - 0.1260*xL - 0.3516*xL^2 + 0.2843*xL^3 - 0.1036*xL^4)
      theta_start = math.atan(yL/(r-xL))
      t = (t-(0.5-LE))/LE -- re-dscretise t to suit circle
      y = r *     math.sin(theta_start + t* math.abs(theta_start))
      x = r *(1 - math.cos(theta_start + t*math.abs(theta_start)))
   end
   return {x=x, y=y}
end
--and finally we can create a paths and normalise them with respct to path length
foil_path = LuaFnPath:new{luaFnName="foil"}
foil_path_norm = ArcLengthParameterizedPath:new{underlying_path=foil_path}

-- to get the far-field bounding boundary, we are going to do some construction
function farField(t)
   -- find a position along the wing surface
   --F = foil(t) -- we can do this by simply evaluating the foil() function
   F = foil_path_norm(t)  -- here F is a table with entries F.x and F.y
   --print("F,",F.x,F.y)
   -- find the local gradient by numerical differentiation 
   delta_t = 0.001
   if t < delta_t then
      tm = t; tp = t+delta_t
   elseif t > 1-delta_t then
      tm = t-delta_t; tp = t      
   else
      tm = t-delta_t; tp = t+delta_t
   end
   Fp = foil_path_norm(tp)
   Fm = foil_path_norm(tm)
   delta_x = Fp.x-Fm.x
   delta_y = Fp.y-Fm.y
   -- for best mesh quality we want grid lines to be perpendicular to the 
   -- wall. This can be achieved by placing points on the farField line
   -- along a line perpendicular to the surface. Lets find the unit vector
   -- for this perpendicular line.
   X = -delta_y; Y = delta_x 
   X = X / (delta_x^2 + delta_y^2)^0.5; Y = Y / (delta_x^2 + delta_y^2)^0.5;
   -- Now find the corresponding far field point by vector addition.
   xval = F.x + X*L
   yval = F.y + Y*L
   return {x=xval, y=yval}
end
-- and to cretae the corresponding path
farField_path_norm = LuaFnPath:new{luaFnName="farField"}

-- Having these two paths defined, we can subivide these to create lines used
-- for mesh construction and also to find the points along the lines. 
b0 = farField_path_norm(0.)  -- we can simply evaluate the path to get coordinates
f1 = farField_path_norm(frac)
f2 = farField_path_norm(0.5)
f3 = farField_path_norm(1.-frac)
t0 = farField_path_norm(1.)

a0 = foil_path_norm(0.)
a1 = foil_path_norm(frac)
a2 = foil_path_norm(0.5)
a3 = foil_path_norm(1.-frac)

b0f1 = SubRangedPath:new{underlying_path=farField_path_norm, t0=0., t1=frac}
f1f2 = SubRangedPath:new{underlying_path=farField_path_norm, t0=frac, t1=0.5}
f2f3 = SubRangedPath:new{underlying_path=farField_path_norm, t0=0.5, t1=1.-frac}
f3t0 = SubRangedPath:new{underlying_path=farField_path_norm, t0=1.-frac, t1=1.}
a0a1 = SubRangedPath:new{underlying_path=foil_path_norm, t0=0., t1=frac}
a1a2 = SubRangedPath:new{underlying_path=foil_path_norm, t0=frac, t1=0.5}
a2a3 = SubRangedPath:new{underlying_path=foil_path_norm, t0=0.5, t1=1.-frac}
a3a0 = SubRangedPath:new{underlying_path=foil_path_norm, t0=1.-frac, t1=1.}

-- the remaining points we can define relative to the others
t1 = Vector3:new{x=t0.x+L, y=t0.y}
b1 = Vector3:new{x=b0.x+L, y=b0.y}
a4 = Vector3:new{x=a0.x+L, y=a0.y}

-- Define patch (which are parametric surfaces, no discretisation at this point.)
surf = {}
surf[0] = CoonsPatch:new{p00=a0, p10=a4, p11=t1, p01=t0}
surf[1] = CoonsPatch:new{north=f3t0, south=a3a0, 
           west=Line:new{p0=a3, p1=f3}, east=Line:new{p0=a0, p1=t0} }
surf[2] = CoonsPatch:new{north=f2f3, south=a2a3, 
           west=Line:new{p0=a2, p1=f2}, east=Line:new{p0=a3, p1=f3} }
surf[3] = CoonsPatch:new{north=f1f2, south=a1a2, 
           west=Line:new{p0=a1, p1=f1}, east=Line:new{p0=a2, p1=f2} }
surf[4] = CoonsPatch:new{north=b0f1, south=a0a1, 
           west=Line:new{p0=a0, p1=b0}, east=Line:new{p0=a1, p1=f1} } -- error ???
surf[5] = CoonsPatch:new{p00=b0, p10=b1, p11=a4, p01=a0}

-- We could stop here if we were to build a 2-D mesh. But for a 3-D mesh
-- an extra surface is needed to go on the end of the wing.
-- Note all mesh sections beyond the wing tips have numbers starting with 1x
surf[16] = AOPatch:new{north=a2a3, south=ReversedPath:new{underlying_path=a0a1}, 
           west=a1a2, east=ReversedPath:new{underlying_path=a3a0} }

-- To create the 3-D mesh we will extrude the surfaces along a vector. 
-- This vector will go from the respective p00 points and extrude in the 
-- +z direction.
-- Volumes that sit around the wing
volume = {} 
volume[0] = SweptSurfaceVolume:new{face0123=surf[0], edge04=Line:new{p0=a0, 
                p1=Vector3:new{x=a0.x, y=a0.y, z=a0.z+S}} }
volume[1] = SweptSurfaceVolume:new{face0123=surf[1], edge04=Line:new{p0=a3,  
                p1=Vector3:new{x=a3.x, y=a3.y, z=a3.z+S}} }
volume[2] = SweptSurfaceVolume:new{face0123=surf[2], edge04=Line:new{p0=a2,  
                p1=Vector3:new{x=a2.x, y=a2.y, z=a2.z+S}} }
volume[3] = SweptSurfaceVolume:new{face0123=surf[3], edge04=Line:new{p0=a1,  
                p1=Vector3:new{x=a1.x, y=a1.y, z=a1.z+S}} }
volume[4] = SweptSurfaceVolume:new{face0123=surf[4], edge04=Line:new{p0=a0,  
                p1=Vector3:new{x=a0.x, y=a0.y, z=a0.z+S}} }
volume[5] = SweptSurfaceVolume:new{face0123=surf[5], edge04=Line:new{p0=b0,  
                p1=Vector3:new{x=b0.x, y=b0.y, z=b0.z+S}} }

-- Volumes that sit beyond wing tip. Here we set the edge04 to start at y = L
volume[10] = SweptSurfaceVolume:new{face0123=surf[0], edge04=Line:new{p0=Vector3:new{x=a0.x, y=a0.y, z=a0.z+S}, p1=Vector3:new{x=a0.x, y=a0.y, z=a0.z+S+L}} }
volume[11] = SweptSurfaceVolume:new{face0123=surf[1], edge04=Line:new{p0=Vector3:new{x=a3.x, y=a3.y, z=a3.z+S}, p1=Vector3:new{x=a3.x, y=a3.y, z=a3.z+S+L}} }
volume[12] = SweptSurfaceVolume:new{face0123=surf[2], edge04=Line:new{p0=Vector3:new{x=a2.x, y=a2.y, z=a2.z+S}, p1=Vector3:new{x=a2.x, y=a2.y, z=a2.z+S+L}} }
volume[13] = SweptSurfaceVolume:new{face0123=surf[3], edge04=Line:new{p0=Vector3:new{x=a1.x, y=a1.y, z=a1.z+S}, p1=Vector3:new{x=a1.x, y=a1.y, z=a1.z+S+L}} }
volume[14] = SweptSurfaceVolume:new{face0123=surf[4], edge04=Line:new{p0=Vector3:new{x=a0.x, y=a0.y, z=a0.z+S}, p1=Vector3:new{x=a0.x, y=a0.y, z=a0.z+S+L}} }
volume[15] = SweptSurfaceVolume:new{face0123=surf[5], edge04=Line:new{p0=Vector3:new{x=b0.x, y=b0.y, z=b0.z+S}, p1=Vector3:new{x=b0.x, y=b0.y, z=b0.z+S+L}} }
volume[16] = SweptSurfaceVolume:new{face0123=surf[16], edge04=Line:new{p0=Vector3:new{x=a1.x, y=a1.y, z=a1.z+S}, p1=Vector3:new{x=a1.x, y=a1.y, z=a1.z+S+L}} }


--            ------f3--------N-------t0---N--t1
--          /       |                 |        |
--         N        |        blk1     |  blk0  |  ny0
--       /  blk2   -a3----\           /        E
--      /         /XXXXXXXX--------- /     S   |
--     f2-------a2XXX NACA FOIL XXXXa0--------a4
--      \         \XXXXXXXX--------- \     N   |
--       \  blk3   -a1----/           \        E
--         N        |        blk4     |  blk5  |  ny0
--     nx0  \       |                 |        |
--            ------f1--------N-------b0---S--b1                                        
--                           nx0          nx1  
--
-- Define number of cells in each block
nx0=61; nx1=30 
ny0=40
nz0=50; nz1=20

-- set up refining function
N_refine = 1
nx0 = math.ceil(N_refine*nx0); nx1 = math.ceil(N_refine*nx1)
ny0 = math.ceil(N_refine*ny0)
nz0 = math.ceil(N_refine*nz0); nz1 = math.ceil(N_refine*nz1)

-- Define Custer Functions.
cfr0 = RobertsFunction:new{end0=true, end1=false, beta=1.02}
cfr1 = RobertsFunction:new{end0=false, end1=true, beta=1.02}
cfz0 = RobertsFunction:new{end0=false, end1=true, beta=1.08} -- z-direction on wing
cfz1 = RobertsFunction:new{end0=true, end1=false, beta=1.05} -- z-direction in far-field
cfx0 = RobertsFunction:new{end0=true, end1=false, beta=1.03}

-- Now we can define the grid!!!  (Hope you have debugging cube)
-- To get an optimum grid you should add some clustering in the wing tangential direction too!
grid = {}
grid[0] = StructuredGrid:new{pvolume=volume[0], niv=nx1, njv=ny0, nkv = nz0, 
            cfList={edge04=cfz0, edge15=cfz0, edge26=cfz0, edge37=cfz0, 
                    edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0,
                    edge01=cfx0, edge32=cfx0, edge76=cfx0, edge45=cfx0} }
grid[1] = StructuredGrid:new{pvolume=volume[1], niv=nx0, njv=ny0, nkv = nz0, 
            cfList={edge04=cfz0, edge15=cfz0, edge26=cfz0, edge37=cfz0, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[2] = StructuredGrid:new{pvolume=volume[2], niv=nx0, njv=ny0, nkv = nz0,
            cfList={edge04=cfz0, edge15=cfz0, edge26=cfz0, edge37=cfz0, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[3] = StructuredGrid:new{pvolume=volume[3], niv=nx0, njv=ny0, nkv = nz0, 
            cfList={edge04=cfz0, edge15=cfz0, edge26=cfz0, edge37=cfz0, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[4] = StructuredGrid:new{pvolume=volume[4], niv=nx0, njv=ny0, nkv = nz0,
            cfList={edge04=cfz0, edge15=cfz0, edge26=cfz0, edge37=cfz0, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[5] = StructuredGrid:new{pvolume=volume[5], niv=nx1, njv=ny0, nkv = nz0, 
            cfList={edge04=cfz0, edge15=cfz0, edge26=cfz0, edge37=cfz0, 
                    edge56=cfr1, edge12=cfr1, edge03=cfr1, edge47=cfr1,
                    edge01=cfx0, edge32=cfx0, edge76=cfx0, edge45=cfx0} }
grid[10] = StructuredGrid:new{pvolume=volume[10], niv=nx1, njv=ny0, nkv = nz1, 
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1, 
                    edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0,
                    edge01=cfx0, edge32=cfx0, edge76=cfx0, edge45=cfx0} }
grid[11] = StructuredGrid:new{pvolume=volume[11], niv=nx0, njv=ny0, nkv = nz1, 
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[12] = StructuredGrid:new{pvolume=volume[12], niv=nx0, njv=ny0, nkv = nz1, 
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[13] = StructuredGrid:new{pvolume=volume[13], niv=nx0, njv=ny0, nkv = nz1, 
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[14] = StructuredGrid:new{pvolume=volume[14], niv=nx0, njv=ny0, nkv = nz1, 
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1, edge56=cfr0, edge12=cfr0, edge03=cfr0, edge47=cfr0} }
grid[15] = StructuredGrid:new{pvolume=volume[15], niv=nx1, njv=ny0, nkv = nz1, 
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1, 
                    edge56=cfr1, edge12=cfr1, edge03=cfr1, edge47=cfr1,
                    edge01=cfx0, edge32=cfx0, edge76=cfx0, edge45=cfx0} }
grid[16] = StructuredGrid:new{pvolume=volume[16], niv=nx0, njv=nx0, nkv = nz1,
            cfList={edge04=cfz1, edge15=cfz1, edge26=cfz1, edge37=cfz1} }
            
-- Define OpenFoam block (a "grid" with labels)
block = {}
block[0] = FoamBlock:new{grid=grid[0],
		     bndry_labels={north="i-00", east="o-00", bottom="s-00"}}
block[1] = FoamBlock:new{grid=grid[1],
		     bndry_labels={north="i-00", south="w-00", bottom="s-00"}}
block[2] = FoamBlock:new{grid=grid[2],
		     bndry_labels={north="i-00", south="w-00", bottom="s-00"}}
block[3] = FoamBlock:new{grid=grid[3],
		     bndry_labels={north="i-00", south="w-00", bottom="s-00"}}
block[4] = FoamBlock:new{grid=grid[4],
		     bndry_labels={north="i-00", south="w-00", bottom="s-00"}}
block[5] = FoamBlock:new{grid=grid[5],
		     bndry_labels={south="i-00", east="o-00", bottom="s-00"}}
block[10] = FoamBlock:new{grid=grid[10],
		     bndry_labels={north="i-00", east="o-00", top="i-00"}}
block[11] = FoamBlock:new{grid=grid[11],
		     bndry_labels={north="i-00", top="i-00"}}
block[12] = FoamBlock:new{grid=grid[12],
		     bndry_labels={north="i-00", top="i-00"}}
block[13] = FoamBlock:new{grid=grid[13],
		     bndry_labels={north="i-00", top="i-00"}}
block[14] = FoamBlock:new{grid=grid[14],
		     bndry_labels={north="i-00", top="i-00"}}
block[15] = FoamBlock:new{grid=grid[15],
		     bndry_labels={south="i-00", east="o-00", top="i-00"}}
block[16] = FoamBlock:new{grid=grid[16],
		     bndry_labels={top="i-00", bottom="w-01"}}
