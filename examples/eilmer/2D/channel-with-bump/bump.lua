-- bump.lua
-- Ported for Eilmer4 by PJ, 2015-05-28
-- This test appeared in papers on Euler solvers for supersonic flows.
-- For example:
-- S. Eidelman, P. Colella and R.P. Shreeve (1984)
-- Application of the Godunov method and its second-order extension to
-- cascade flow modelling.
-- AIAA Journal Vol. 22 No. 11 pp  

config.title = "Channel with circular-arc bump."
print(config.title)

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
p_inf = 101.3e3 -- Pa
T_inf = 288.0 -- degree K
a_inf = math.sqrt(1.4*287.0*T_inf)
u_inf = 1.65 * a_inf -- m/s
print("Inflow: p=", p_inf, " T=", T_inf, " velx=", u_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0}

-- Geometry of flow domain.
L = 1.0 -- will use for both length and height of domain
h = 0.04 * L -- height of bump
a0 = Vector3:new{x=0.0}; a1 = Vector3:new{y=L}
b0 = Vector3:new{x=L}; b1 = Vector3:new{x=L,y=L}
c0 = Vector3:new{x=1.5*L,y=h}
d0 = Vector3:new{x=2.0*L}; d1 = Vector3:new{x=2.0*L,y=L}
e0 = Vector3:new{x=3.0*L}; e1 = Vector3:new{x=3.0*L,y=L}
patch0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
patch1 = makePatch{north=Line:new{p0=b1,p1=d1},
		   east=Line:new{p0=d0,p1=d1},
		   south=Arc3:new{p0=b0,pmid=c0,p1=d0},
		   west=Line:new{p0=b0,p1=b1}}
patch2 = CoonsPatch:new{p00=d0, p10=e0, p11=e1, p01=d1}

-- Mesh the patches, with particular discretisation.
rcfx0 = RobertsFunction:new{end0=false, end1=true, beta=1.2}
rcfx1 = RobertsFunction:new{end0=true, end1=true, beta=1.2}
rcfx2 = RobertsFunction:new{end0=true, end1=false, beta=1.2}
rcfy = RobertsFunction:new{end0=true, end1=false, beta=1.2}
ni0 = 64; nj0 = 64 -- We'll scale discretization off these values
factor = 1.0
ni0 = math.floor(ni0*factor); nj0 = math.floor(nj0*factor)
grid0 = StructuredGrid:new{psurface=patch0,
			   cfList={north=rcfx0,east=rcfy,south=rcfx0,west=rcfy},
			   niv=ni0+1, njv=nj0+1}
grid1 = StructuredGrid:new{psurface=patch1,
			   cfList={north=rcfx1,east=rcfy,south=rcfx1,west=rcfy},
			   niv=ni0+1, njv=nj0+1}
grid2 = StructuredGrid:new{psurface=patch2,
			   cfList={north=rcfx2,east=rcfy,south=rcfx2,west=rcfy},
			   niv=ni0+1, njv=nj0+1}
-- Define the flow-solution blocks and set boundary conditions.
fba0 = FBArray:new{grid=grid0, initialState=inflow, nib=16, njb=4,
                   bcList={west=InFlowBC_Supersonic:new{flowState=inflow}}}
fba1 = FBArray:new{grid=grid1, initialState=inflow, nib=16, njb=4}
fba2 = FBArray:new{grid=grid2, initialState=inflow, nib=16, njb=4,
                   bcList={east=OutFlowBC_Simple:new{}}}
identifyBlockConnections()

config.block_marching = true
config.nib = 48
config.njb = 4
config.propagate_inflow_data = true
config.flux_calculator = "adaptive"
config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 0.8
config.max_time = 20.0*L/u_inf -- long enough, tunnel has 15ms steady time
config.max_step = 50000
config.dt_init = 1.0e-6
config.dt_plot = config.max_time
