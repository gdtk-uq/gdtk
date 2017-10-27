-- arc.lua
-- This test appeared in papers on Euler solvers for transonic flows.
-- For example:
--   TL Holst and WF Ballhaus (1978)
--   Fast conservative schemes for the full potential equation
--   applied to transonic flows.
--   NASA TM-78469
-- PJ, 2017-10-27, adapted from channel-with-bump example

config.title = "Circular-arc foil in transonic flow."
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
p_stag = 101.3e3 -- Pa
T_stag = 300.0 -- degree K
stagGas = FlowState:new{p=p_stag, T=T_stag}
M_inf = 0.84 -- for supercritical flow over foil
p_inf = p_stag/idealgasflow.p0_p(M_inf)
T_inf = T_stag/idealgasflow.T0_T(M_inf)
a_inf = math.sqrt(1.4*287.1*T_inf)
velx_inf = M_inf * a_inf
print("Free-stream: p=", p_inf, "T=", T_inf, "velx=", velx_inf)
initGas = FlowState:new{p=p_inf, T=T_inf, velx=velx_inf}

-- Geometry of flow domain.
L = 1.0 -- will use for both length and height of domain
h = 0.10 * L -- height of bump
a0 = Vector3:new{x=0.0}; a1 = Vector3:new{y=6*L}
b0 = Vector3:new{x=5*L}; b1 = Vector3:new{x=5*L,y=6*L}
c0 = Vector3:new{x=5.5*L,y=h}
d0 = Vector3:new{x=6*L}; d1 = Vector3:new{x=6*L,y=6*L}
e0 = Vector3:new{x=11*L}; e1 = Vector3:new{x=11*L,y=6*L}
patch0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
patch1 = makePatch{north=Line:new{p0=b1,p1=d1},
		   east=Line:new{p0=d0,p1=d1},
		   south=Arc3:new{p0=b0,pmid=c0,p1=d0},
		   west=Line:new{p0=b0,p1=b1}}
patch2 = CoonsPatch:new{p00=d0, p10=e0, p11=e1, p01=d1}

-- Mesh the patches, with particular discretisation.
rcfx0 = RobertsFunction:new{end0=false, end1=true, beta=1.1}
rcfx1 = RobertsFunction:new{end0=true, end1=true, beta=1.3}
rcfx2 = RobertsFunction:new{end0=true, end1=false, beta=1.1}
rcfy = RobertsFunction:new{end0=true, end1=false, beta=1.1}
ni0 = 24; nj0 = 24 -- We'll scale discretization off these values
factor = 2.0
ni0 = math.floor(ni0*factor); nj0 = math.floor(nj0*factor)
grid0 = StructuredGrid:new{psurface=patch0,
			   cfList={north=rcfx0,east=rcfy,south=rcfx0,west=rcfy},
			   niv=3*ni0+1, njv=3*nj0+1}
grid1 = StructuredGrid:new{psurface=patch1,
			   cfList={north=rcfx1,east=rcfy,south=rcfx1,west=rcfy},
			   niv=ni0+1, njv=3*nj0+1}
grid2 = StructuredGrid:new{psurface=patch2,
			   cfList={north=rcfx2,east=rcfy,south=rcfx2,west=rcfy},
			   niv=3*ni0+1, njv=3*nj0+1}
-- Define the flow-solution blocks and set boundary conditions.
blk0 = FluidBlockArray{grid=grid0, initialState=initGas, nib=3, njb=2,
		       bcList={west=InFlowBC_FromStagnation:new{stagnationState=stagGas}}}
blk1 = FluidBlockArray{grid=grid1, initialState=initGas, nib=1, njb=2}
blk2 = FluidBlockArray{grid=grid2, initialState=initGas, nib=3, njb=2,
		       bcList={east=OutFlowBC_FixedP:new{p_outside=p_inf}}}
identifyBlockConnections()

config.flux_calculator = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 0.8
config.max_time = 0.400 -- 100.0*L/velx_inf=0.366
config.max_step = 50000
config.dt_init = 1.0e-6
config.dt_plot = config.max_time/40
