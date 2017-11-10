-- read the case and find the vertex pairs
config.title = "3D Connection tester"
config.dimensions = 3

setGasModel('ideal-air.lua')
config.max_time = 0.6e-3
config.max_step = 1
config.dt_init = 1.0e-6
config.flux_calculator = "ausmdv"

nn = 4
p_base = 1.0e5
T_base = 348.3

-- Set up the physical aspects of the test
function lftCubeFill(x, y, z)
   p_inc_factor = x * (5*y^2) * (10*z^3)
   p = p_base*(1.0 + p_inc_factor)
   return FlowState:new{p=p, T=T_base}
end

function rghtCubeFill(x, y, z)
   p_inc_factor = x * (5*y^2) * (10*z^3)
   p = p_base*(0.8 + p_inc_factor)
   return FlowState:new{p=p, T=T_base}
end

-- Set up a function to give vertex positions of a cube
--[[
          H------------G
         /|           /|
        / |          / |
       D--+---------C  |
       |  |         |  |
       |  |         |  |
       |  E---------+--F
       | /          | /
       |/           |/
       A------------B

--]]
function cube(origin, length)
   x0 = origin.x; y0 = origin.y; z0 = origin.z
   A = Vector3:new{x=x0, y=y0, z=z0}
   B = Vector3:new{x=x0+length, y=y0, z=z0}
   C = Vector3:new{x=x0+length, y=y0+length, z=z0}
   D = Vector3:new{x=x0, y=y0+length, z=z0}
   E = Vector3:new{x=x0, y=y0, z=z0+length}
   F = Vector3:new{x=x0+length, y=y0, z=z0+length}
   G = Vector3:new{x=x0+length, y=y0+length, z=z0+length}
   H = Vector3:new{x=x0, y=y0+length, z=z0+length}
   return {A=A, B=B, C=C, D=D, E=E, F=F, G=G, H=H}
end

-- build cubes (left and right)
cL = cube(Vector3:new{x=0.0, y=0.0, z=0.0}, 1.0)
lftVtxs = {cL.A, cL.B, cL.C, cL.D, cL.E, cL.F, cL.G, cL.H}
cR = cube(Vector3:new{x=1.0, y=0.0, z=0.0}, 1.0)
rghtVtxs = {cR.A, cR.B, cR.C, cR.D, cR.E, cR.F, cR.G, cR.H}

lftGrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices=lftVtxs},
			     niv=nn+1, njv=nn+1, nkv=nn+1}
lftBlk = FluidBlock:new{grid=lftGrid, initialState=lftCubeFill}

rghtGrid = StructuredGrid:new{pvolume=TFIVolume:new{vertices=rghtVtxs},
			      niv=nn+1, njv=nn+1, nkv=nn+1}
rghtBlk = FluidBlock:new{grid=rghtGrid, initialState=rghtCubeFill}

identifyBlockConnections()


   
   




