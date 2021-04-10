-- cubic-ramp.lua
-- PJ, 2017-05-16 adapted from Eilmer3 example
-- Model of Mohammadian's concave surface experiment. 

config.title = "Mohammadian cubic ramp."
print(config.title)
-- Conditions to match those reported in JFM paper.
p_inf = 66.43 -- Pa
u_inf = 1589.8 -- m/s
T_inf = 41.92 -- degree K
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=p_inf/5, T=T_inf, velx=0.0}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
T_wall = 296.0 -- degree K, assumed cold-wall temperature

function ramp(t)
   -- Parametric definition of ramp profile in xy-plane.
   local m_per_inch = 25.4e-3 -- metres per inch
   local alpha = 28.0*math.pi/180.0 -- angle of final straight section
   local tan_alpha = math.tan(alpha)
   local x_join_inch = math.sqrt(50.0*tan_alpha)
   local y_join_inch = x_join_inch^3 / 150.0
   local x_inch = 6.4 * t
   local y_inch = y_join_inch + (x_inch - x_join_inch) * tan_alpha
   if x_inch < x_join_inch then
      y_inch = x_inch^3 / 150.0
   end
   return {x=x_inch*m_per_inch, y=y_inch*m_per_inch, z=0.0}
end

mm = 1.0e-3 -- metres per mm
a0 = Vector3:new(ramp(0.0)); a1 = a0+Vector3:new{y=5*mm} -- leading edge
b0 = Vector3:new(ramp(1.0)); b1 = b0+Vector3:new{x=-5*mm, y=20*mm} -- downstream end
c0 = b0+Vector3:new{x=10*mm, y=0.0}; c1 = Vector3:new{x=c0.x, y=b1.y} -- end of model
patch0 = CoonsPatch:new{north=Line:new{p0=a1,p1=b1},
			south=LuaFnPath:new{luaFnName="ramp"},
			west=Line:new{p0=a0, p1=a1},
			east=Line:new{p0=b0, p1=b1}}
patch1 = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}

rcfx = RobertsFunction:new{end0=true, end1=false, beta=1.2}
rcfy = RobertsFunction:new{end0=true, end1=false, beta=1.1}
ni0 = 200; nj0 = 40 -- We'll scale discretization off these values
factor = 1.0
ni0 = math.floor(ni0*factor); nj0 = math.floor(nj0*factor)
grid0 = StructuredGrid:new{psurface=patch0, niv=ni0+1, njv=nj0+1,
			   cfList={north=rcfx,east=rcfy,south=rcfx,west=rcfy}}
grid1 = StructuredGrid:new{psurface=patch1, niv=math.floor(ni0/10)+1, njv=nj0+1,
			   cfList={east=rcfy, west=rcfy}}

wedge = FBArray:new{grid=grid0, nib=10, njb=2,
			initialState=inflow, label="wedge",
			bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
				south=WallBC_NoSlip_FixedT:new{Twall=T_wall,group="loads"},
				west=InFlowBC_Supersonic:new{flowState=inflow}}}
tail = FBArray:new{grid=grid1, nib=1, njb=2,
		       initialState=initial, label="tail",
		       bcList={north=InFlowBC_Supersonic:new{flowState=inflow},
			       south=WallBC_NoSlip_FixedT:new{Twall=T_wall},
			       east=OutFlowBC_FixedP:new{p_outside=p_inf/5}}}
identifyBlockConnections()

-- Do a little more setting of global data.
config.viscous = true
config.flux_calculator = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 1.0
config.max_time = 1.0e-3 -- long enough for several flow lengths
config.max_step = 2000000
config.dt_init = 1.0e-9
config.dt_plot = 0.1e-3
config.dt_loads = 0.1e-3
