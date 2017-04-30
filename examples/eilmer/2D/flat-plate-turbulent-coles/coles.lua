-- coles.lua : Turbulent flow over a flat plate
-- Anand V  2016-05-26, Indooroopilly
-- Ported PJ's example from Eilmer3
-- 
config.title = "Coles Mach 3.7 flow over a flat plate (k-omega)"
print(config.title)
-- Some of the global config options affect the way blocks and 
-- boundary conditions are built.  They need to be set first.
config.dimensions = 2
config.turbulence_model = "k_omega"
config.viscous = true
config.flux_calculator = 'adaptive'
 	
-- Gas model and flow conditions
nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
p_inf = 1.358e3  -- Pa
u_inf = 677.4    -- m/s
T_inf = 83.34    -- K
gas_inf = FlowState:new{p=p_inf, T=T_inf}:toTable()
-- Turbulence quantities estimate
tke_inf = 1.5 * (u_inf * 0.05)^2
mu_t_inf = 100.0 * gas_inf.mu
omega_inf = gas_inf.rho * tke_inf / mu_t_inf
-- print("Inflow turbulence: tke=", tke_inf, "omega=", omega_inf)
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, tke=tke_inf, omega=omega_inf}
print("Inflow Check\n", inflow)

-- Geometry of the flow domain
L = 0.60 -- metres
H = 0.4 * L
--
--         wall
--        c---------b
-- flow=> |         |
--        d         |
--          -\-     |
--    flow=>    -\- |
--        0         a ----> x
-- 
a = Vector3:new{x=L, y=0.0}; b = Vector3:new{x=L, y=H}; 
c = Vector3:new{x=0.0, y=H}; d = Vector3:new{x=0.0, y=3.0*H/4.0}
patch = CoonsPatch:new{p00=d, p10=a, p11=b, p01=c}
cfx = RobertsFunction:new{end0=true,end1=false,beta=1.05}
cflist = {north=cfx, east=RobertsFunction:new{end0=false,end1=true,beta=1.01},
	  south=cfx, west=RobertsFunction:new{end0=false,end1=true,beta=1.05}}
grd = StructuredGrid:new{psurface=patch, niv=100, njv=90, cfList=cflist}

blks = FluidBlockArray{grid=grd, nib=4, njb=1, fillCondition=inflow,
		       bcList={north=WallBC_NoSlip_Adiabatic:new{},
			       east=OutFlowBC_Simple:new{},
			       south=InFlowBC_Supersonic:new{flowCondition=inflow},
			       west=InFlowBC_Supersonic:new{flowCondition=inflow}}}
identifyBlockConnections()

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- should allow a few flow lengths   
config.dt_plot =  1.0e-3
config.dt_history = 1.0e-5
config.max_step = 200 

config.gasdynamic_update_scheme = "classic-rk3"
config.cfl_value = 0.5	
config.cfl_count = 3
config.stringent_cfl = false -- true is more robust
config.dt_init = 1.0e-9	-- only an initial guess, the simulation will take this over
