-- rmi.lua
-- Redesign by NNG, February 2021
config.title = "Richtmyer-Meshkov Instability"
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)

-- The initial flow is travelling leftwards at some velocity, so that the instability 
-- stays relatively still after the impact. This gives us more time to watch it evolve. 
-- 300 m/s seemed to be about the right number, after some trial and error
axes_velocity = 300
state1 = FlowState:new{p=100.0e3, T=300.0, velx=0.0-axes_velocity, vely=0.0} -- pre-shock
state2 = FlowState:new{p=497.8e3, T=531.1, velx=469.9-axes_velocity, vely=0.0} -- post-shock
state3 = FlowState:new{p=100.0e3, T=75.0, velx=0.0-axes_velocity, vely=0.0} -- dense gas
config.flux_calculator = 'ausmdv'
config.flow_format='rawbinary'


-- Geometry of flow domain.
L = 4.0  -- length of downstream part of duct
H = 1.0  -- height of duct
d = 1.5  -- interface location
delta = 0.1 -- magnitude of interface perturbation
E1= 15.0    -- Length of sponge block 1
E2= 15.0    -- Length of sponge block 2

-- Interface shape is controlled by a custom fill condition.
function fillFunction(x,y,z)
   -- Parametric path with 0<=t<=1.
   interfacex = delta*math.cos(2*math.pi*(y/H)) + d
   if (x < interfacex) then
       return state1
   else
       return state3
   end
end

spongepatch1= CoonsPatch:new{p00=Vector3:new{x=-E1,  y=0.0},
                             p10=Vector3:new{x=0.0, y=0.0},
                             p11=Vector3:new{x=0.0, y=H  },
                             p01=Vector3:new{x=-E1,  y=H  }}

mainpatch = CoonsPatch:new{p00=Vector3:new{x=0.0, y=0.0},
                           p10=Vector3:new{x=L,   y=0.0},
                           p11=Vector3:new{x=L,   y=H  },
                           p01=Vector3:new{x=0.0, y=H  }}

spongepatch2= CoonsPatch:new{p00=Vector3:new{x=L,   y=0.0},
                             p10=Vector3:new{x=L+E2, y=0.0},
                             p11=Vector3:new{x=L+E2, y=H  },
                             p01=Vector3:new{x=L,   y=H  }}

nx = 200
nsponge = math.floor((nx+1)/4)
ny = math.floor((nx+1)*(H/L))
dx_main = L/nx
a_sponge1 = dx_main/E1
a_sponge2 = dx_main/E2

cluster1 = GeometricFunction:new{a=a_sponge1, r=1.14, N=nsponge, reverse=true}
cluster2 = GeometricFunction:new{a=a_sponge2, r=1.14, N=nsponge, reverse=false}
cflist1 = {north=cluster1, east=none, south=cluster1, west=none}
cflist2 = {north=cluster2, east=none, south=cluster2, west=none}
sponge1grid = StructuredGrid:new{psurface=spongepatch1, niv=nsponge, njv=ny, cfList=cflist1}
maingrid   = StructuredGrid:new{psurface=mainpatch, niv=nx+1, njv=ny}
sponge2grid = StructuredGrid:new{psurface=spongepatch2, niv=nsponge, njv=ny, cfList=cflist2}

-- Define the boundary conditions that we care about 
-- and the flow-solution blocks.
sponge1bcs = {north=nil,
             east=nil,
             south=nil,
             west=InFlowBC_Supersonic:new{flowState=state2}}

sponge1blk = FBArray:new{grid=sponge1grid, nib=1, njb=1,
		                initialState=state2, bcList=sponge1bcs}

mainbcs = {north=nil,
           east=nil,
           south=nil,
           west=nil}
mainblk = FBArray:new{grid=maingrid, nib=4, njb=1,
                      initialState=fillFunction, bcList=mainbcs}

sponge2bcs = {north=nil,
              east=OutFlowBC_SimpleExtrapolate:new{},
              south=nil,
              west=nil}
spong2blk = FBArray:new{grid=sponge2grid, nib=1, njb=1,
                      initialState=state3, bcList=sponge2bcs}

identifyBlockConnections()

config.max_time = 200.0e-3  -- seconds
config.dt_plot = 5.0e-4
config.max_step = 10000000
config.dt_init = 5.0e-6
