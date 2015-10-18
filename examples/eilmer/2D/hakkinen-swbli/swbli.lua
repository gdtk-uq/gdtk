-- swbli.lua
-- Anand V, 10-October-2015
-- Adapted PJ's script of the model of Hakkinen et al's 1959 experiment

config.title = "Shock Wave Boundary Layer Interaction"
print(config.title)
config.dimensions = 2

-- Conditions to match those of Figure 6: pf/p0=1.4, Re_shock=2.96e5
p_inf = 6205.0 -- Pa
u_inf = 514.0 -- m/s
T_inf = 164.4 -- degree K

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)

inflow = FlowState:new{p=p_inf, velx=u_inf, T=T_inf, vely=0.0}
mm = 1.0e-3 -- metres per mm
L1 = 10.0*mm; L2 = 90.0*mm; L3 = 67*mm
H1 = 37.36*mm
alpha = 3.09*math.pi/180.0 -- angle of inviscid shock generator
tan_alpha = math.tan(alpha)
a0 = Vector3:new{-L1, 0.0}; a1 = a0+Vector3:new{0.0,H1} -- leading edge of shock generator
b0 = Vector3:new{0.0, 0.0}; b1 = b0+Vector3:new{0.0,H1-L1*tan_alpha} -- start plate
c0 = Vector3:new{L3, 0.0}; c1 = c0+Vector3:new{0.0,H1-(L1+L3)*tan_alpha} -- end shock generator
d0 = Vector3:new{L2, 0.0}; d1 = d0+Vector3:new{0.0,H1} -- end plate

-- Number of cells, blocks and clustering
rcf = RobertsFunction:new{end0=true,end1=true,beta=1.1}
ni0 = 20; nj0 = 80 -- We'll scale discretization off these values
factor = 4
ni0 = math.floor(ni0*factor); nj0 = math.floor(nj0*factor)

gridin = StructuredGrid:new{psurface=CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1},
			    cfList={nil,rcf,nil,rcf}, niv=ni0+1, njv=nj0+1}

gridp1 = StructuredGrid:new{psurface=CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1},
			    cfList={nil,rcf,nil,rcf}, niv=7*ni0+1, njv=nj0+1}

gridp2 = StructuredGrid:new{psurface=CoonsPatch:new{p00=c0, p10=d0, p11=d1, p01=c1},
			    cfList={nil,rcf,nil,rcf}, niv=2*ni0+1, njv=nj0+1}

blkin = SBlockArray{grid=gridin, fillCondition=inflow, nib=1, njb=2,
		    bcList={west=InFlowBC_Supersonic:new{flowCondition=inflow},
			    north=WallBC_WithSlip:new{},
			    south=WallBC_WithSlip:new{}}}

blk1 = SBlockArray{grid=gridp1, fillCondition=inflow, nib=7, njb=2,
		   bcList={south=WallBC_NoSlip_Adiabatic:new{},
			   north=WallBC_WithSlip:new{}}}

blk2 = SBlockArray{grid=gridp2, fillCondition=inflow, nib=2, njb=2,
		   bcList={south=WallBC_NoSlip_Adiabatic:new{},
			   north=WallBC_WithSlip:new{},
			   east=OutFlowBC_FixedPT:new{p_outside=p_inf,T_outside=T_inf}}}
identifyBlockConnections()

config.gasdynamic_update_scheme = "classic-rk3"
config.flux_calculator = 'adaptive'
config.viscous = true
config.cfl_value = 1.0
config.max_time = 5.0*L2/u_inf -- time in flow lengths
config.max_step = 200000
config.dt_init = 1.0e-8
config.dt_plot = config.max_time/10
