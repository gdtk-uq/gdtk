-- bcp.lua
-- PJ, 2016-05-25

config.title = "X3 blunted conical probe."
print(config.title)
config.axisymmetric = true
config.viscous = true

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
inflow = FlowState:new{p=1000.0, T=622.0, velx=5000.0, vely=0.0}
initial = FlowState:new{p=100.0, T=300.0, velx=0.0, vely=0.0}

theta = 15.0*math.pi/180.0
tanth = math.tan(theta)
mm = 1.0e-3 -- mm in metres
-- Put virtual apex at (0,0)
Rnose = 0.3*mm; xnose = Rnose/tanth
Rbase = 5.0*mm; xbase = Rbase/tanth
a0 = Vector3:new{x=xnose/2, y=0.0}
a1 = Vector3:new{x=xnose/2, y=Rnose}
a2 = Vector3:new{x=xnose/2, y=2.5*Rnose}
b0 = Vector3:new{x=xnose, y=0.0}
b1 = Vector3:new{x=xnose, y=Rnose}
b2 = Vector3:new{x=xnose-0.5*tanth*1.5*Rnose, y=2.5*Rnose}
patch0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
patch1 = CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2}
c1 = Vector3:new{x=xbase, y=Rbase}
c2 = c1 + Vector3:new{x=-tanth*3.0*mm, y=3.0*mm}
patch2 = CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2}
cfx0 = RobertsFunction:new{end0=false,end1=true,beta=1.2}
cfy0 = RobertsFunction:new{end0=false,end1=true,beta=1.3}
cfy1 = RobertsFunction:new{end0=true,end1=false,beta=1.2}
cfx2 = RobertsFunction:new{end0=true,end1=false,beta=1.2}
factor = 2
ni = math.floor(12*factor); nj = math.floor(16*factor)
grd0 = StructuredGrid:new{psurface=patch0, niv=ni+1, njv=nj+1,
			  cfList={west=cfy0,east=cfy0,north=cfx0,south=cfx0}}
grd1 = StructuredGrid:new{psurface=patch1, niv=ni+1, njv=nj+1,
			  cfList={west=cfy1,east=cfy1,north=cfx0,south=cfx0}}
ni2 = math.floor(200*factor)
grd2 = StructuredGrid:new{psurface=patch2, niv=ni2+1, njv=nj+1,
			  cfList={west=cfy1,east=cfy1,north=cfx2,south=cfx2}}

-- Assemble the block from the grid and boundary data.
blk0 = SBlock:new{grid=grd0, fillCondition=inflow,
		  bcList={west=InFlowBC_Supersonic:new{flowCondition=inflow}}}
blk1 = SBlock:new{grid=grd1, fillCondition=inflow,
		  bcList={west=InFlowBC_Supersonic:new{flowCondition=inflow},
			  north=InFlowBC_Supersonic:new{flowCondition=inflow}}}
blks = SBlockArray{grid=grd2, nib=10, njb=1, 
		   fillCondition=initial,
		   bcList={north=InFlowBC_Supersonic:new{flowCondition=inflow},
			   east=OutFlowBC_Simple:new{},
			   south=WallBC_NoSlip_FixedT:new{Twall=300.0}}}
identifyBlockConnections()

config.flux_calculator = "adaptive"
config.stringent_cfl = true
config.max_time = 10.0e-6  -- seconds
config.max_step = 500000
config.dt_init = 1.0e-11
config.dt_plot = 0.1e-6
